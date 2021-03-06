% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MedLRT.R
\name{MedLRT}
\alias{MedLRT}
\title{Likelihood ratio test for multiple mediators}
\usage{
MedLRT(
  X,
  M,
  Y,
  Z,
  total_iter = 1000,
  tol = 1e-06,
  total_iter_lambda = 10,
  m = 10000,
  seed = 100
)
}
\arguments{
\item{X}{a vector of \eqn{n} observations for one exposure}

\item{M}{a matrix of \eqn{n} observations for \eqn{q} mediators}

\item{Y}{a vector of \eqn{n} observations for one outcome}

\item{Z}{a matrix of \eqn{n} observations for \eqn{L} confounding covariates
including the first column with all ones for the intercept}

\item{total_iter}{the maximum number of iterations for updating regression
coefficients under the null. The default value is 1000.}

\item{tol}{the convergence tolerance. The default value is 1e-6}

\item{total_iter_lambda}{the maximum number of iteration for updating
Lagrange multiplier \eqn{\lambda}}

\item{m}{the Monte Carlo sample size for simulating the null distribution of
LRT. The default value is 10000.}

\item{seed}{the random seed. The default value is 100.}
}
\value{
a list object of five elements:
\describe{
\item{pvalue}{p value of the likelihood ratio test}
\item{chisq_pvalue}{p value based on the chi-squared distribution}
\item{kappa_pvalue}{p value based on the kappa distribution}
\item{test_stat}{likelihood ratio test statistics}
\item{convergence}{a binary indicator. Value 1 indicates
 the convergence of the algorithm. Value 0 indicates that
 the algorithm does not converge}
}
}
\description{
Likelihood ratio test for multiple mediators
}
\examples{
### A simulated example with no significant mediation effects
n = 100
q = 10
L = 3
X = rnorm(n)
Z = cbind(1,matrix(rnorm(n*L),nrow=n,ncol=L-1))
M = X\%*\%t(rep(c(1,0,1),length=q)) + 0.1*matrix(rnorm(n*q),nrow=n,ncol=q)
Y = X + M\%*\%rep(c(0,1,0),length=q) + Z\%*\%matrix(1,nrow=L,ncol=1) + 0.1*rnorm(n)
res <- MedLRT(X=X,M=M,Y=Y,Z=Z)
print(as.data.frame(res),digits=3)
}
\author{
Wei Hao <weihao@umich.edu>
}
