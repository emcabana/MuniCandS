// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double calc_est(NumericMatrix X, IntegerVector H) {
  int n = X.nrow();
  double T = 0.0;
  
  for (int h = 0; h < n; ++h) {
    for (int i = 0; i < n; ++i) {
      double prod = 1.0;
      for (int j = 0; j < H.size(); ++j) {
        int col = H[j] - 1; // Convertir de 1-based (R) a 0-based (C++)
        double val_h = X(h, col);
        double val_i = X(i, col);
        double term = (pow(val_h, 2.0) + pow(val_i, 2.0)) / 2.0 - std::max(val_h, val_i) + 1.0 / 3.0;
        prod *= term;
      }
      T += prod;
    }
  }
  
  return T / n;
}



// [[Rcpp::export]]
double calc_est_arma(const arma::mat& X, const arma::uvec& H) {
  int n = X.n_rows;
  double T = 0.0;

  for (int h = 0; h < n; ++h) {
    for (int i = 0; i < n; ++i) {
      double prod = 1.0;
      for (size_t j = 0; j < H.n_elem; ++j) {
        double val_h = X(h, H[j]);
        double val_i = X(i, H[j]);
        double term = (std::pow(val_h, 2.0) + std::pow(val_i, 2.0)) / 2.0 - std::max(val_h, val_i) + 1.0 / 3.0;
        prod *= term;
      }
      T += prod;
    }
  }

  return T / n;
}

// [[Rcpp::export]]
NumericVector S2Cpi(NumericVector x) {
	int p = x.size() - 1;
	NumericVector phi(p);
	for (int j = 0; j < p; ++j) {
    double denom = 1.0;
    for (int k = 0; k < j; ++k) {
      denom *= sin(phi[k]);
    }
    phi[j] = acos(x[j] / (denom == 0 ? 1.0 : denom));
  }
  if ( x[p] < 0 ) {
        phi[p-1]=2*M_PI-phi[p-1];
    } 
  return phi;
}



// -------------------------------------------------------------
// intlin: interpolación lineal escalar (equivalente al intlin de tu R)
// -------------------------------------------------------------
 // [[Rcpp::export]]
double intlin(double x, NumericVector X, NumericVector Y) {
  int n = X.size();
  if (n == 0) return NA_REAL;

  // extremos
  if (x <= X[0]) return Y[0];
  if (x >= X[n-1]) return Y[n-1];

  // buscar primer índice pos tal que X[pos] >= x
  int pos = std::lower_bound(X.begin(), X.end(), x) - X.begin();

  // si coincide exactamente
  if (pos < n && X[pos] == x) return Y[pos];

  int i1 = std::max(0, pos - 1);
  int i2 = std::min(n - 1, pos);

  double x1 = X[i1], x2 = X[i2];
  double y1 = Y[i1], y2 = Y[i2];

  if (x2 == x1) return y1; // degenerado
  return ((x2 - x) * y1 + (x - x1) * y2) / (x2 - x1);
}

// -----------------------------------------------------------------------------
// ajus: función principal (incluye la lógica de 'preaj' integrada)
// devuelve pchisq(x, df = 1) según tu versión corregida
// -----------------------------------------------------------------------------
 // [[Rcpp::export]]
double ajus(double y, NumericVector Y) {
  int Len = Y.size();
  if (Len < 2) return NA_REAL;

  // indi = (1:Len)/(Len+1)
  NumericVector indi(Len);
  for (int i = 0; i < Len; ++i) indi[i] = (i + 1.0) / (Len + 1.0);

  // X = qchisq(indi, df = 1)
  NumericVector X(Len);
  for (int i = 0; i < Len; ++i) X[i] = R::qchisq(indi[i], 1.0, /*lower_tail*/1, /*log_p*/0);

  // L1 = round(.3*Len), L2 = round(.95*Len)
  int L1 = std::lround(0.3 * Len);
  int L2 = std::lround(0.95 * Len);

  // acotar entre 1 y Len (R es 1-based)
  if (L1 < 1) L1 = 1;
  if (L2 < 1) L2 = 1;
  if (L1 > Len) L1 = Len;
  if (L2 > Len) L2 = Len;
  if (L2 < L1) L2 = L1;

  // convertir a 0-based
  int start = L1 - 1;
  int end = L2 - 1;

  // XM = X[L1:L2], YM = Y[L1:L2]
  int m = end - start + 1;
  NumericVector XM(m), YM(m);
  for (int i = 0; i < m; ++i) {
    XM[i] = X[start + i];
    YM[i] = Y[start + i];
  }

  // medias MX, MY
double MX = std::accumulate(XM.begin(), XM.end(), 0.0) / XM.size();
double MY = std::accumulate(YM.begin(), YM.end(), 0.0) / YM.size();
  // bb = sum((XM-MX)*(YM-MY))/sum((XM-MX)^2)
  double num = 0.0, den = 0.0;
  for (int i = 0; i < m; ++i) {
    double dx = XM[i] - MX;
    double dy = YM[i] - MY;
    num += dx * dy;
    den += dx * dx;
  }
  double bb = (den == 0.0) ? 0.0 : num / den;
  double aa = MY - bb * MX;

  // Construir YY = Y ajustada (misma longitud que Y)
  NumericVector YY = clone(Y);

// límites del bloque de mezcla
double XM_first = XM[0];                  // primer elemento
double XM_last  = XM[XM.size() - 1];      // último elemento
double denom_mix = XM_last - XM_first;

  for (int j = 0; j < Len; ++j) {
    double Xj = X[j];
    if (Xj < XM_first) {
      YY[j] = Y[j];
    } else if (Xj > XM_last) {
      YY[j] = aa + bb * Xj;
    } else {
      // mezcla entre Y[j] y aa+bb*Xj
      double v1 = Y[j];
      double v2 = aa + bb * Xj;
      if (denom_mix == 0.0) {
        // fallback (muy raro: XM_first == XM_last)
        YY[j] = 0.5 * (v1 + v2);
      } else {
        double w1 = (XM_last - Xj);
        double w2 = (Xj - XM_first);
        YY[j] = (w1 * v1 + w2 * v2) / denom_mix;
      }
    }
  }

  // ahora determinar x: si y > YY[Len]  o bien intlin
  double YY_last = YY[Len - 1];
  double xval;
  if (y > YY_last) {
    if (bb == 0.0) {
      // fallback si bb == 0: evitamos division por 0
      xval = y - aa;
    } else {
      xval = (y - aa) / bb;
    }
  } else {
    xval = intlin(y, YY, X);
  }

  // devolver pchisq(x, df = 1)
  double pval = R::pchisq(xval, 1.0, /*lower_tail*/1, /*log_p*/0);
  return 1.0-pval;
  
//  return xval;
}

// [[Rcpp::export]]
NumericVector Cpi2C(NumericVector phi) {
  int p = phi.size();
  NumericVector u(p);

  for (int j = 0; j < p - 1; ++j) {
    // Compute I_j(phi_j) and normalize
    double num = 0.0, denom = 0.0;
    int N = 1000;
    double h = M_PI / N;
    for (int i = 0; i <= N; ++i) {
      double t = i * h;
      double w = (i == 0 || i == N) ? 0.5 : 1.0;
      double f = pow(sin(t), p - j - 1);
      if (t <= phi[j]) num += w * f;
      denom += w * f;
    }
    u[j] = h * num / (h * denom);
  }

  // Final coordinate: φ_p / (2π)
  u[p - 1] = phi[p - 1] / (2.0 * M_PI);

  return u;
}


