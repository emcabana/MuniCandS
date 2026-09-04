// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;


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
// intlin: interpolación lineal escalar
// -----------------------------------------------------------
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

//--------------------------------------------------------------
// ajus: función principal 
// devuelve pchisq(x, df = 1)
//--------------------------------------------------------------
 // [[Rcpp::export]]
// -----------------------------------------------------------------------------
// ajus: versión con parámetro entero h
// devuelve el valor pv definido en tu fórmula
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double ajus(double y, Rcpp::NumericVector Y, int h) {
  int R = Y.size();
  if (R < 1) return NA_REAL;

  // Construir YS = (0, sort(Y), Inf)
  Rcpp::NumericVector YS(R + 2);
  YS[0] = 0.0;
  std::vector<double> Ys_sorted(Y.begin(), Y.end());
  std::sort(Ys_sorted.begin(), Ys_sorted.end());
  for (int i = 0; i < R; ++i) {
    YS[i + 1] = Ys_sorted[i];
  }
  YS[R + 1] = R_PosInf;

  // r = max índice tal que YS[r] <= y
  int r = 0;
  for (int i = 0; i < R + 2; ++i) {
    if (YS[i] <= y) r = i;
  }

  // parámetros de la gamma
  double a  = std::pow(5.0, h) / std::pow(2.0, h + 1);
  double la = std::pow(15.0, h) / 2.0;

  // valores de distribución gamma
  double p1 = R::pgamma(YS[r],   a, 1.0/la, /*lower_tail*/1, /*log_p*/0);
  double p2 = R::pgamma(YS[r+1], a, 1.0/la, /*lower_tail*/1, /*log_p*/0);
  double pp = R::pgamma(y,       a, 1.0/la, /*lower_tail*/1, /*log_p*/0);

  // cálculo del p‑valor
  double pv = (R - r + (p2 - pp) / (p2 - p1)) / (R + 1);

  return pv;
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


