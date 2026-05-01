#' @title Multivariate tests of uniformity, normality, spherical and elliptical symmetry, and independence
#'
#' @description Implements uniformity tests on \eqn{C_p = [0,1]^p} and \eqn{S^p},
#' as well as tests of normality, isotropy (spherical symmetry), ellipticity (elliptical symmetry),
#' and independence on \eqn{R^p}.
#'
#' @details The core of the method is a fast implementation of the m- and s-tests
#' introduced in the reference.
#' Given a multivariate sample, the function either:
#' \itemize{
#'   \item directly computes p-values for uniformity tests,
#'   \item transforms the sample from \eqn{S^p} to \eqn{C_p} for uniformity on \eqn{S^p},
#'   \item transforms the sample from \eqn{R^p} to \eqn{C_p} for normality, isotropy,
#'         ellipticity, or independence, and then computes p-values for the m- and s-tests.
#' }
#'
#' The main test statistics are the squared norms \eqn{b_H^2} of the zero-marginal random
#' components of the empirical process associated with subsets
#' \eqn{H \subset J := \{1,2,\dots,p\}}. By default, all nonempty subsets are included,
#' unless restricted by \code{hmin} and \code{hmax}.
#' The \eqn{b_H^2} are asymptotically independent with distribution function \eqn{P_H},
#' so that the vector of p-values \eqn{1 - P_H(b_H^2)} is asymptotically uniform on
#' \eqn{[0,1]^{2^p-1}}.
#'
#' A first Monte Carlo simulation estimates the finite-sample distribution of the \eqn{b_H^2},
#' yielding p-values and the combined statistics
#' \deqn{p_{min} = 1 - (1 - \min(pvals))^{\mathrm{length}(pvals)}}
#' and
#' \deqn{p_{sum} = 1 - \mathrm{pchisq}(\sum(\mathrm{qchisq}(1 - pvals, df = 1)), df = \mathrm{length}(pvals))},
#' both asymptotically uniform on \eqn{[0,1]} under the null hypothesis.
#'
#' A second Monte Carlo simulation of \eqn{p_{\min}} and \eqn{p_{\mathrm{sum}}}
#' corrects their null distribution by accounting for the finite-sample dependence
#' among the p-values.
#' @param Z Numeric \eqn{n \times p} matrix. Each row represents one sample element.
#' @param type Character string specifying the null hypothesis:
#'   \itemize{
#'     \item `"UC"` — Uniformity on \eqn{[0,1]^p}.
#'     \item `"US"` — Uniformity on \eqn{S^{p-1}}.
#'     \item `"N"` — Normality on \eqn{R^p}.
#'     \item `"I"` — Isotropy (spherical symmetry) on \eqn{R^p}.
#'     \item `"E"` — Elliptical symmetry on \eqn{R^p}.
#'     \item `"IN"` — Independence on \eqn{R^p}.
#'   }
#' @param hmin Integer. Minimum subset size to include.
#' @param hmax Integer. Maximum subset size to include (use \code{Inf} for no upper bound).
#' @param n_sim Integer. Number of Monte Carlo repetitions for the first simulation.
#' @param n_mc Integer. Number of Monte Carlo repetitions for the second simulation.
#' @param use_parallel Logical. If \code{TRUE}, parallelization is used.
#' @param cache Optional list containing precomputed Monte Carlo simulations.
#' @param return_cache Logical. If TRUE, the Monte Carlo simulations are returned for reuse.
#' @param graph Logical. If \code{TRUE}, plots the vector of p-values
#'   and prints the list of subsets used in the computation.
#' @param full Logical. If \code{TRUE}, includes all p-values resulting
#'   from the first Monte Carlo simulation in the returned values.
#'
#' @return A list of length 2 containing the p-values of the
#' m- and s-tests. If \code{full = TRUE}, the return also includes the intermediate
#' p-values obtained before the second Monte Carlo simulation.
#'
#' @references
#' Cabana, A. and Cabana, E. M. (2025).
#' Brownian sheet and uniformity tests on the hypercube.
#' arXiv:2509.06134.
#'
#' @examples
#' \donttest{
#' Z <- matrix(runif(150), ncol = 3)
#' MuniCandS(Z, type = "UC", n_sim=100, n_mc = 100)
#' }
##########################################
### FUNCION PRINCIPAL 
##########################################
#' @export
MuniCandS <- function(
  Z, type,
  hmin = 1, hmax = Inf,
  n_sim = 1000,
  n_mc = 1000,
  use_parallel = FALSE,
  cache = NULL,
  return_cache = FALSE,
  full = FALSE,
  graph = FALSE
) {

  stopifnot(is.matrix(Z))
  stopifnot(type %in% c("UC","US","N","I","E","IN"))


  n <- nrow(Z)
  p <- ncol(Z)
  px <- if (type == "US") p - 1 else p

  H_list <- genlist(hmin, hmax, px)

  # --- validar cache ---
  valid_cache <- FALSE
  if (!is.null(cache)) {
    if (!is.null(cache$params)) {
      pars <- cache$params
      valid_cache <- all.equal(
        list(n=n, p=p, hmin=hmin, hmax=hmax, type=type, n_sim=n_sim, n_mc=n_mc),
        pars
      ) == TRUE
    }
  }

  # --- backend ---
  lapply_fun <- lapply
  if (use_parallel) {
    lapply_fun <- function(X, FUN) {
      future.apply::future_lapply(X, FUN, future.seed = TRUE)
    }
  }

  # --- BH ---
  if (valid_cache && !is.null(cache$BH)) {
    BH <- cache$BH
  } else {
    BH <- do.call(rbind,
      lapply_fun(1:n_sim, function(i) {
        X2bH(Z2X(G2Z(n, p, type), type), H_list)
      })
    )
  }

  # --- PV ---
  if (valid_cache && !is.null(cache$PV)) {
    PV <- cache$PV
  } else {
    PV <- do.call(rbind,
      lapply_fun(1:n_mc, function(i) {
        pvals2pv(
          bHBH2pvals(
            X2bH(Z2X(G2Z(n, p, type), type), H_list),
            BH
          ),
          H_list
        )
      })
    )
  }

  # --- datos observados ---
  pvals <- bHBH2pvals(
    X2bH(Z2X(Z, type), H_list),
    BH
  )

  pvmys0 <- pvals2pv(pvals, H_list)
  pvmys  <- pvPV2mys(pvmys0, PV)

  names(pvmys) <- c("m-test", "s-test")

  if (graph) {
    plot(pvals, ylim = c(0, 1), pch = 15)
    print(H_list)
  }

  new_cache <- list(
    BH = BH,
    PV = PV,
    params = list(
      n = n, p = p,
      hmin = hmin, hmax = hmax,
      type = type,
      n_sim = n_sim, n_mc = n_mc
    )
  )

  if (return_cache) {
    return(list(result = pvmys, cache = new_cache))
  }

  if (!full) return(list(pvals=pvmys))
  return(list(pvals.all=pvmys0, pvals=pvmys))
}

##########################################
### FUNCIONES BASICAS 
##########################################

G2Z=function(n,p,type){
	 if(type %in% c("E","N")){
  	ZS <- matrix(stats::rnorm(n*p),n,p)
  }
  if(type %in% c("UC","IN")){
    ZS <- matrix(runif(n * p), n, p)
    }
  if(type%in% c("US","I")){
  	 ZS=matrix(stats::rnorm(n*p),n,p)
  	 if(type=="US")ZS=ZS/sqrt(rowSums(ZS^2))
  }
return(ZS)
}

Z2X=function(Z,type){
	  n <- nrow(Z)
  p <- ncol(Z)
  X <- Z
  variance <- NA

  if (type == "US") {
    X <- S2C(Z) #matrix(NA, n, p-1)
  #  for (i in 1:n) X[i, ] <- S2C(Z[i, ])
  }

  if (type == "E") {
    ZC <- scale(Z, center = TRUE, scale = FALSE)
    variance <- crossprod(ZC) / n
    eig <- eigen(variance)
    Z <- ZC %*% eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)/sqrt(n)
  }

  if (type %in% c("I","E")) {
    rho <- sqrt(rowSums(Z^2))
    Z0 <- Z / rho
    X <- matrix(NA, n, p-1)
    for (i in 1:n) X[i, ] <- S2C(Z0[i, ])
    X <- cbind(X, rank(rho) / (n + 1))
  }

  if (type == "N") {
    ZC <- scale(Z, center = TRUE, scale = FALSE)
    variance <- crossprod(ZC) / n
    eig <- eigen(variance)
    X <- pnorm(ZC %*% eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)/sqrt(n))
  }

 if (type == 'IN') {
  	for(j in 1:p) X[,j]=rank(X[,j])/(n+1)
  }

  return(X)
}

X2bH=function(X,H_list){
	X <- as.matrix(X)
	bH <- sapply(H_list, function(H) calc_est_arma(X, as.integer(H) - 1))
	# Armadillo usa 0-based
	return(bH)
}

bHBH2pvals=function(bH,BH){
  resu <- numeric(ncol(BH))
  for (i in seq_len(ncol(BH))) {
    resu[i] <- ajus(bH[i], sort(BH[, i]))
  }
  resu[resu < 0] <- 0
  return(resu)
}

pvals2pv=function(pvals,H_list){
	pmin=min(pvals)
	suma=sum(qchisq(1 - pvals, df = 1))
  
	m_p_value <- 1 - (1 - pmin)^length(H_list)
	s_p_value <- 1 - pchisq(suma, df = length(H_list))
	return(c(m_p_value, s_p_value))
}

pvPV2mys=function(pv,PV){
	mys=c(evlin(pv[1],sort(PV[,1])),evlin(pv[2],sort(PV[,2])))
	return(mys)
}

##########################################
### FUNCIONES AUXILIARES 
##########################################

evlin=function(u,f){ff=c(0,sort(f),1)
	if (u==1) res <- 1 else{sf=max(which(ff<=u))
	res <- ((sf+(u-ff[sf])/(ff[sf+1]-ff[sf])-1)/(length(f)+1))}
	return(res)
}

   # Generation of the subsets H of J={0,1,...,p}
  genlist=function(hmin,hmax,p){
  H_list <- list()
  for (h in max(1,hmin):min(hmax, p)) {
    H_list <- c(H_list, combn(p, h, simplify = FALSE))
  }
  return(H_list)
}

 S2C <- function(Z) {
  if (is.matrix(Z)) {
    n <- nrow(Z)
    U <- matrix(NA, n, ncol(Z)-1)
    for (i in 1:n) U[i, ] <- Cpi2C(S2Cpi(Z[i, ]))
  } else {
    U <- Cpi2C(S2Cpi(Z))
  }
  return(U)
}

Cpi2S=function(Phi){
	Phi=c(Phi,0)
	Z=cos(Phi[1])
	for(j in 2:(length(Phi)))Z=c(Z,cos(Phi[j])*prod(sin(Phi[1:(j-1)])))
	return(Z)
}


