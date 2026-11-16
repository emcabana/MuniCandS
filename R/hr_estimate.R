#' Estimador simultáneo de Hettmansperger-Randles
#' @param X matriz n x p
#' @param tol tolerancia de convergencia
#' @param maxit máximo de iteraciones
#' @export
hr_estimate <- function(X, tol = 1e-6, maxit = 500)
  hr_cpp(as.matrix(X), tol, as.integer(maxit))