// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// Estimador simultaneo de Hettmansperger-Randles (centro + forma, det V = 1)
// [[Rcpp::export]]
Rcpp::List hr_cpp(const arma::mat& X, double tol = 1e-6, int maxit = 500) {
  const uword n = X.n_rows, p = X.n_cols;
  rowvec mu = median(X, 0);              // arranque: mediana coordenada a coordenada
  mat V = eye(p, p), Vh(p, p), Vmh(p, p), Xc(n, p), Z(n, p), U(n, p), M(p, p);
  vec r(n), rinv(n), eval; mat evec;
  bool conv = false; int k = 0;

  auto standardize = [&]() {
    if (!eig_sym(eval, evec, V) || eval.min() <= 0)
      Rcpp::stop("hr_cpp: V no es definida positiva (datos degenerados?)");
    Vh  = evec * diagmat(sqrt(eval)) * evec.t();
    Vmh = evec * diagmat(1.0 / sqrt(eval)) * evec.t();
    Xc = X; Xc.each_row() -= mu;
    Z = Xc * Vmh;
    r = sqrt(sum(Z % Z, 1));
  };

  for (; k < maxit; ++k) {
    standardize();
    uword neff = 0;
    for (uword i = 0; i < n; ++i) {
      if (r[i] > 1e-12) { rinv[i] = 1.0 / r[i]; ++neff; } else rinv[i] = 0.0;
    }
    U = Z; U.each_col() %= rinv;                       // u_i (0 si z_i = 0)
    rowvec s = sum(U, 0);
    M = (static_cast<double>(p) / neff) * (U.t() * U);

    if (abs(s).max() / neff < tol && abs(M - eye(p, p)).max() < tol) {
      conv = true; break;
    }
    mu += (s / sum(rinv)) * Vh;                        // paso de la mediana espacial
    V = Vh * M * Vh;                                   // paso de Tyler
    V = 0.5 * (V + V.t());
    double ld, sgn; log_det(ld, sgn, V);
    V *= std::exp(-ld / p);                            // det V = 1
  }
  if (!conv) { Rcpp::warning("hr_cpp: no convergio en %d iteraciones", maxit); standardize(); }

  return Rcpp::List::create(
    Rcpp::_["mu"] = mu, Rcpp::_["V"] = V, Rcpp::_["Vh"] = Vh,
    Rcpp::_["Z"] = Z, Rcpp::_["r"] = r,
    Rcpp::_["iter"] = k, Rcpp::_["converged"] = conv);
}