#' @title Multivariate test of uniformity, normality and isotropy on C and S
#'
#' @description Performs uniformity tests on \eqn{C_p=[0,1]^p} and S^p,
#' and normality and isotropy tests on R^p
#'
#' @details
#' The central piece is a fast implementation of the m- and s-tests
#' introduced by the authors in 
#' \emph{Brownian sheet and uniformity tests on the hyperecube}, arXiv 2509.06134.
#' It receives a multivariate sample and either
#' computes directly the p-values of the uniformity tests, or
#' transforms the sample from S^p to C_p for uniformity on S^p,
#' or from R^p to C_p for normality and for isotropy, and
#' computes the p-values for the m- and s-tests.
#' MuniCandS performs the original tests,
#' MUCS corrects the p-values for finite samples
#' by using the results of 10000 simulations of MuniCandS.
#' @param Z n x p matrix. Each row is an element of the sample.
#' @param type Character string specifying the null hypothesis:
#'   \itemize{
#'     \item `"UC"` — Uniformity on the p-dimensional cube C_p.
#'     \item `"US"` — Uniformity on the (p-1)-dimensional sphere.
#'     \item `"N"` — Normality on R^p.
#'     \item `"I"` — Isotropy on R^p.
#'   }
#' @param hmin Integer, minimum subset size to consider.
#' @param hmax Integer, maximum subset size to consider (use Inf for no upper bound).
#' @param repet Integer, number of Monte Carlo repetitions used to estimate the p-value.
#' @param semilla Integer, random seed for reproducibility.
#'
#' @return A named numeric vector with the p-values of the m- and s-tests.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' Z <- matrix(runif(150), ncol = 3)
#' MuniCandS(Z, type = "UC")
#' }
#'
#' @export
MuniCandS <- function(Z, type, hmin = 1, hmax = Inf, repet = 999, semilla = 123) {
  if (!is.matrix(Z)) stop("Z debe ser una matriz")

  n <- nrow(Z)
  
  datap=dataproc(Z,type)
  X=datap$X
  variance=datap$variance
  p <- ncol(X)
  H_list=genlist(hmin,hmax,p)
  if(type=='I') H_list[[p]] <- NULL
  T_H2=calcest(X,H_list)
  pvals=calcpvals(T_H2,n,p,repet,H_list,variance=datap$variance,type=type)
  
  
  pmin=min(pvals)
  suma=sum(qchisq(1 - pvals, df = 1))
  
  m_p_value <- 1 - (1 - pmin)^length(H_list)
  s_p_value <- 1 - pchisq(suma, df = length(H_list))
  return(c(m_p_value = m_p_value, s_p_value = s_p_value))
}


  # Initial processing of the data  
  
dataproc <- function(Z, type) {
  n <- nrow(Z)
  p <- ncol(Z)
  X <- Z
  variance <- NA

  if (type == "US") {
    X <- matrix(NA, n, p-1)
    for (i in 1:n) X[i, ] <- S2C(Z[i, ])
  }

  if (type == "I") {
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
    X <- pnorm(ZC %*% eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors))
  }

  list(X = X, variance = variance)
}

   # Generation of the subsets H of J={0,1,...,p}
  genlist=function(hmin,hmax,p){
  H_list <- list()
  for (h in max(1,hmin):min(hmax, p)) {
    H_list <- c(H_list, combn(p, h, simplify = FALSE))
  }
  return(H_list)
}

  #Computation of the statistics ||T_H||^2
  
calcest <- function(X, H_list) {
  X <- as.matrix(X)
  sapply(H_list, function(H) calc_est_arma(X, as.integer(H) - 1)) # Armadillo usa 0-based
}

	
  # Simulation of samples under the null hypothesis
   
simsamp <- function(n, p, variance = diag(p), type) {
  if (type == "N") {
    X <- mvtnorm::rmvnorm(n, rep(0, p), variance)
    XC <- scale(X, center = TRUE, scale = FALSE)
    eig <- eigen(variance)
    X <- pnorm(XC %*% eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors))
  } else {
    X <- matrix(runif(n * p), n, p)
    if (type == "I") X[, p] <- rank(X[, p]) / (n + 1)
  }
  X
}

  # Estimation of the p-values by simulation
calcpvals <- function(T_H2, n, p, repet, H_list, variance = diag(p), type) {
  sim_stats <- parallel::mclapply(1:repet, function(rep) {
    samp <- simsamp(n, p, variance, type)
    calcest(samp, H_list)
  }, mc.cores = parallel::detectCores())

  sim_matrix <- do.call(rbind, sim_stats)
  resu=c()
  for(i in 1:ncol(sim_matrix))resu <- c(resu,ajus(T_H2[i],sort(sim_matrix[,i])))
  resu[resu<0]=0
  return(resu)
}



S2C <- function(Z) {
  if (is.matrix(Z)) {
    n <- nrow(Z)
    U <- matrix(NA, n, ncol(Z))
    for (i in 1:n) U[i, ] <- Cpi2C(S2Cpi(Z[i, ]))
  } else {
    U <- Cpi2C(S2Cpi(Z))
  }
  U
}

Cpi2S=function(Phi){
	Phi=c(Phi,0)
	Z=cos(Phi[1])
	for(j in 2:(length(Phi)))Z=c(Z,cos(Phi[j])*prod(sin(Phi[1:(j-1)])))
	return(Z)
}

evlin=function(u,f){ff=c(0,sort(f),1)
	sf=max(which(ff<=u))
	if(sf==length(ff))return(1) else return(((u-ff[sf])*(sf)+(ff[sf+1]-u)*(sf-1))/((ff[sf+1]-ff[sf])*(length(f)+1)))
 }


#'
#' @export
MUCS=function(X, type, hmin = 1, hmax = Inf, repet = 999, semilla = 123, full=FALSE){
	pv=MuniCandS(X, type=type, hmin = hmin, hmax = hmax, repet = repet, semilla = semilla)
	n=nrow(X)
	p=ncol(X)
	if(type=='US')p=p-1
	p=min(p,hmax)
	pvd=pvdata
	if(type=='N')pvd=pvdatan
	if(type=='I')pvd=pvdatai
	pvs=matrix(NA,3,2)
	for(i in 1:3){
pvs[i,1]=evlin(pv[1],pvd[,1,i,p-1])
pvs[i,2]=evlin(pv[2],pvd[,2,i,p-1])
	}

evv=c(intlin(n,c(50,100,200),pvs[,1]),intlin(n,c(50,100,200),pvs[,2]))
	
if(full==TRUE) return(c(pv,evv)) else return(evv)
}
