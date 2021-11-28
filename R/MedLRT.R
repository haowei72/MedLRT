
#'@import mvtnorm


#OBTAIN MLES AND LOG LIKELIHOOD FN USING LEAST SQUARES W/O CONSTRAINT
HM_MLE_H1=function(dat)
{
  n=length(dat$Y)
  q=ncol(dat$M)
  L=ncol(dat$Z)

  Y_fit=lm(dat$Y ~ dat$W-1) #no intercept, the intercept is included in Z
  M_fit=lm(dat$M ~ dat$B-1) #no intercept, the intercept is included in Z

  alpha_tilde_un=coefficients(M_fit)#matrix: (L+1)*q
  beta_tilde_un=coefficients(Y_fit)#col vector:1*(1+q+L)

  beta_star=rbind(t(beta_tilde_un[1:q]),matrix(0,nrow=L,ncol=q))
  alpha_star=matrix(c(alpha_tilde_un[1,],rep(0,(L+1))),ncol=1)

  eps_Y=dat$Y-dat$W%*%beta_tilde_un
  sigma2_Y=mean(eps_Y*eps_Y)

  eps_M=dat$M - dat$B%*%alpha_tilde_un
  sigma2_M=t(eps_M)%*%eps_M/n

  logdet=as.numeric(determinant(sigma2_M,logarithm=TRUE)$modulus)
  log_lik=-0.5*n*(log(sigma2_Y)+1+logdet+q)

  theta = list(alpha_tilde_un = alpha_tilde_un, beta_tilde_un= beta_tilde_un,
               sigma2_Y = sigma2_Y, sigma2_M = sigma2_M,logdet=logdet,
               beta_star=beta_star,alpha_star=alpha_star)

  return(list(theta=theta,log_lik=log_lik))
}

########CALCULATE TWO LAMBDAS##################
cal_lambda=function(theta,dat)#theta is initial theta
{
  q=ncol(dat$M)

  a1_T=theta$alpha_tilde_un[1,]
  a2=theta$beta_tilde_un[(1:q)]

  b1_T=(dat$B_inv%*%theta$beta_star%*%theta$sigma2_M)[1,]
  b2=as.matrix((theta$sigma2_Y*dat$W_inv%*%theta$alpha_star)[(1:q),1])

  b=a1_T%*%b2+b1_T%*%a2

  lambda_1=(b-sqrt(b^2-4*b1_T%*%b2%*%a1_T%*%a2))/(2*b1_T%*%b2)
  lambda_2=(b+sqrt(b^2-4*b1_T%*%b2%*%a1_T%*%a2))/(2*b1_T%*%b2)


  return(list(lambda_1=lambda_1,lambda_2=lambda_2))
}


########GIVEN LAMBDA,CALCULATE ALPHA_TILDE, BETA_TILDE###########
alpha_tilde_beta_tilde=function(theta,dat,lambda)
{

  lambda=as.numeric(lambda)
  alpha_tilde=theta$alpha_tilde_un-lambda*dat$B_inv%*%theta$beta_star%*%theta$sigma2_M
  beta_tilde=theta$beta_tilde_un-lambda*theta$sigma2_Y*dat$W_inv%*%theta$alpha_star

  return(list(alpha_tilde=alpha_tilde,beta_tilde=beta_tilde))
}

# ########GIVEN ALPHA_TILDE AND BETA_TILDE,CALCULATE LOG-LIKELIHOOD###########
HM_2_log_likeli=function(alpha_tilde,beta_tilde,sigma2_Y,logdet,inv_sigma2_M,dat)
{
  n =  nrow(dat$Y)

  eps_Y = dat$Y-dat$W%*%beta_tilde

  eps_M = dat$M - dat$B%*%alpha_tilde
  sigma2_M = t(eps_M)%*%eps_M

  log_lik=-0.5*(n*log(sigma2_Y)+n*logdet+sum(eps_Y*eps_Y)/sigma2_Y+sum(inv_sigma2_M*sigma2_M))
  return(log_lik)
}

##########SELECT LAMBDA THAT GIVES THE HIGHER LOG-LIKELIHOOD############
sel_lambda=function(theta,lambda,dat)
{

  inv_sigma2_M=solve(theta$sigma2_M)
  logdet = as.numeric(determinant(theta$sigma2_M,logarithm=TRUE)$modulus)

  alpha_tilde_beta_tilde_est_1=alpha_tilde_beta_tilde(theta=theta,dat=dat,
                                                      lambda=lambda$lambda_1)


  log_like_lambda_1=HM_2_log_likeli(alpha_tilde=alpha_tilde_beta_tilde_est_1$alpha_tilde,
                                    beta_tilde=alpha_tilde_beta_tilde_est_1$beta_tilde,
                                    sigma2_Y=theta$sigma2_Y,logdet=logdet,inv_sigma2_M=inv_sigma2_M,
                                    dat=dat)


  alpha_tilde_beta_tilde_est_2=alpha_tilde_beta_tilde(theta=theta,dat=dat,
                                                      lambda=lambda$lambda_2)

  log_like_lambda_2=HM_2_log_likeli(alpha_tilde=alpha_tilde_beta_tilde_est_2$alpha_tilde,
                                    beta_tilde=alpha_tilde_beta_tilde_est_2$beta_tilde,
                                    sigma2_Y=theta$sigma2_Y,logdet=logdet,inv_sigma2_M=inv_sigma2_M,
                                    dat=dat)


  if (log_like_lambda_1>log_like_lambda_2){
    lambda=lambda$lambda_1
  }
  else{
    lambda=lambda$lambda_2
  }

  return(lambda=as.numeric(lambda))
}

#############GIVEN THE SELECTED LAMBDA, SIGMA2_Y AND SIGMA2_M, UPDATE ALPHA_TILDE, BETA_TILDE#######
cal_para_fix_lambda_null=function(lambda,theta,dat,total_iter,tol)#intial theta
{
  n =  nrow(dat$Y)
  q =  ncol(dat$M)
  L =  ncol(dat$Z)

  alpha_tilde_un=theta$alpha_tilde_un
  beta_tilde_un=theta$beta_tilde_un
  sigma2_Y=theta$sigma2_Y
  sigma2_M=theta$sigma2_M
  alpha_star0=theta$alpha_star
  beta_star0=theta$beta_star
  alpha_tilde0=alpha_tilde_un
  beta_tilde0=beta_tilde_un

  for(iter in 1:total_iter){

    alpha_tilde=alpha_tilde_un-lambda*dat$B_inv%*%beta_star0%*%sigma2_M
    beta_tilde=beta_tilde_un-lambda*sigma2_Y*dat$W_inv%*%alpha_star0

    beta_star=rbind(t(beta_tilde[1:q]),matrix(0,nrow=L,ncol=q))
    alpha_star=matrix(c(alpha_tilde[1,],rep(0,(L+1))),ncol=1)

    err = mean((alpha_tilde-alpha_tilde0)^2)+mean((beta_tilde-beta_tilde0)^2)


    if(err<tol){
      break
    }#once break, we will not calculate the four lines below

    beta_star0=beta_star
    alpha_star0=alpha_star
    alpha_tilde0=alpha_tilde
    beta_tilde0=beta_tilde
  }
  return(theta=list(alpha_tilde_un=alpha_tilde_un,beta_tilde_un=beta_tilde_un,
                    alpha_tilde=alpha_tilde,beta_tilde=beta_tilde,
                    alpha_star=alpha_star,beta_star=beta_star,
                    sigma2_Y=sigma2_Y,sigma2_M=sigma2_M))
}


################################FIND PARAS UNDER H0###################
HM_MLE_H0=function(theta,dat,total_iter,tol,total_iter_lambda)
{

  q=ncol(dat$M)
  n=nrow(dat$Y)

  theta0=theta

  convergence = 0

  for(iter in 1:total_iter_lambda)
  {

    lambda=cal_lambda(theta=theta0,dat=dat)

    lambda_sel=sel_lambda(theta=theta0,dat=dat,lambda=lambda)

    theta1=cal_para_fix_lambda_null(lambda=lambda_sel,
                                    theta=theta0,dat=dat,
                                    total_iter=total_iter,tol=tol)#update alpha and beta

    eps_Y = dat$Y-dat$W%*%theta1$beta_tilde#update variance
    theta1$sigma2_Y = mean(eps_Y*eps_Y)

    eps_M = dat$M - dat$B%*%theta1$alpha_tilde
    theta1$sigma2_M = t(eps_M)%*%eps_M/n

    err = mean((theta1$alpha_tilde-theta0$alpha_tilde)^2)+mean((theta1$beta_tilde-theta0$beta_tilde)^2)

    theta0=theta1

    if(abs(err)<tol)
    {
      convergence=1
      break
    }

  }

  theta=list(alpha_tilde=theta0$alpha_tilde,
             beta_tilde=theta0$beta_tilde,
             alpha_star=theta0$alpha_star,
             beta_star=theta0$beta_star,
             sigma2_Y=theta0$sigma2_Y,
             sigma2_M=theta0$sigma2_M,
             lambda=lambda_sel)

  logdet = as.numeric(determinant(theta$sigma2_M,logarithm=TRUE)$modulus)
  log_lik=-0.5*n*(log(theta$sigma2_Y)+1+logdet+q)

  return(list(theta=theta,log_lik=log_lik,convergence=convergence,
              iter_lambda = iter))

}

#given unconstrained sigma2_M and sigma2_Y (mle) to find i_theta_0
find_eigen_val=function(dat,sigma2_M,sigma2_Y)
{
  n = length(dat$Y)
  q = ncol(dat$M)
  L = ncol(dat$Z)
  V=cbind(dat$Z, dat$X)

  B_tilde=dat$B[,-1]
  M_fit_1 = lm(dat$M ~ B_tilde-1) #no intercept, the intercept is included in Z
  alpha_tilde=coefficients(M_fit_1)# matrix: L*q
  alpha_tilde=rbind(0,alpha_tilde)

  B_prod=t(dat$B)%*%dat$B
  sigma2_M_inv=solve(sigma2_M)
  i_theta_0_1=cbind(kronecker(sigma2_M_inv,B_prod)/n,matrix(0,nrow=(L+1)*q,ncol=q+L+1))

  alpha_t_B_t=t(alpha_tilde)%*%t(dat$B)
  B_alpha=dat$B%*%alpha_tilde
  E_WW_1=alpha_t_B_t%*%B_alpha/n+sigma2_M
  E_WW_2=alpha_t_B_t%*%V/n
  E_WW_3=t(V)%*%B_alpha/n
  E_WW_4=t(V)%*%V/n
  E_WW=rbind(cbind(E_WW_1,E_WW_2),cbind(E_WW_3,E_WW_4))/sigma2_Y
  i_theta_0_2=cbind(matrix(0,nrow=q+L+1,ncol=(L+1)*q),E_WW)
  I_theta_0=rbind(i_theta_0_1,i_theta_0_2)

  I_tilde=rbind(cbind(diag(q),matrix(0,nrow=q,ncol=L+1)),matrix(0,nrow=L*q,ncol=q+L+1))
  H_theta_0=rbind(cbind(matrix(0,nrow=(L+1)*q,ncol=(L+1)*q),I_tilde),cbind(t(I_tilde),matrix(0,nrow=L+1+q,ncol=L+1+q)))

  eig_res = eigen(I_theta_0)
  I_theta_0_inv_sqrt = tcrossprod(eig_res$vectors%*%diag(1/sqrt(eig_res$values)),eig_res$vectors)
  A_theta=I_theta_0_inv_sqrt%*%H_theta_0%*%I_theta_0_inv_sqrt
  eigen_val=eigen(A_theta)$value
  eigen=eigen_val[c(1:q,(L*q+q+L+2):(L*q+2*q+L+1))]
  return(eigen)

}

p_val_00=function(R,eig,m = 10000)
{
  q2 = length(eig)
  xi_mat = matrix(rchisq(m*q2,df=1),nrow=m,ncol=q2)
  Tstat = (xi_mat%*%eig)^2/(xi_mat%*%eig^2)/4
  return(mean(Tstat>R))
}


############LR METHOD###################
LR_pval=function(dat,total_iter,tol,total_iter_lambda,m)
{
  #H1
  H1_fit=HM_MLE_H1(dat=dat)
  log_lik_H1=H1_fit$log_lik #calculate log likelihood under H1
  theta_0=H1_fit$theta#get initial value for theta in H0

  #H0
  H0_fit=HM_MLE_H0(theta=theta_0,dat=dat,total_iter = total_iter,tol=tol,total_iter_lambda=total_iter_lambda )
  log_lik_H0=H0_fit$log_lik#calculate log likelihood under H0


  #calculate test statistics using chisq(1)
  R=2*(log_lik_H1-log_lik_H0)
  LR.pval=pchisq(R, df=1, lower=F)

  #calculate test statistics using 00 null distribution
  eigen=find_eigen_val(dat=dat,sigma2_M=theta_0$sigma2_M,sigma2_Y=theta_0$sigma2_Y)
  LR00.pval=p_val_00(R=R,eig=eigen,m = m)

  return(list(pvalue = max(LR.pval,LR00.pval),
              chisq_pvalue = LR.pval,
              kappa_pvalue = LR00.pval,
              test_stat = R,
              convergence=H0_fit$convergence))
}


#####################
#' Likelihood ratio test for multiple mediators
#' @param X a vector of \eqn{n} observations for one exposure
#' @param M a matrix of \eqn{n} observations for \eqn{q} mediators
#' @param Y a vector of \eqn{n} observations for one outcome
#' @param Z a matrix of \eqn{n} observations for \eqn{L} confounding covariates
#' including the first column with all ones for the intercept
#' @param total_iter the maximum number of iterations for updating regression
#' coefficients under the null. The default value is 1000.
#' @param tol the convergence tolerance. The default value is 1e-6
#' @param total_iter_lambda the maximum number of iteration for updating
#' Lagrange multiplier \eqn{\lambda}
#' @param m the Monte Carlo sample size for simulating the null distribution of
#' LRT. The default value is 10000.
#' @param seed the random seed. The default value is 100.
#' @return a list object of two elements
#' \describe{
#' \item{pvalue}{p value of the likelihood ratio test}
#' \item{chisq_pvalue}{p value based on the chi-squared distribution}
#' \item{kappa_pvalue}{p value based on the kappa distribution}
#' \item{test_stat}{likelihood ratio test statistics}
#' \item{convergence}{a binary indicator. Value 1 indicates
#'  the convergence of the algorithm. Value 0 indicates that
#'  the algorithm does not converge}
#' }
#' @author Wei Hao <weihao@umich.edu>
#' @examples
#' n = 100
#' q = 10
#' L = 3
#' X = rnorm(n)
#' Z = cbind(1,matrix(rnorm(n*L),nrow=n,ncol=L-1))
#' M = X%*%t(rep(c(1,0,1),length=q)) + 0.1*matrix(rnorm(n*q),nrow=n,ncol=q)
#' Y = X + M%*%rep(c(0,1,0),length=q) + Z%*%matrix(1,nrow=L,ncol=1) + 0.1*rnorm(n)
#' res <- MedLRT(X=X,M=M,Y=Y,Z=Z)
#' print(as.data.frame(res),digits=3)
#' @export
MedLRT<-function(X,M,Y,Z,total_iter=1000,tol=1e-6,
                 total_iter_lambda=10,m=10000,seed=100)
{
  W =cbind(M,Z,X)
  B=cbind(X,Z)
  B_inv=solve(t(B)%*%B)
  W_inv=solve(t(W)%*%W)
  dat = list(X = X, M=M, Y=Y,Z=Z, W = W, B=B, B_inv=B_inv,W_inv=W_inv)

  LRres = LR_pval(dat=dat,total_iter=total_iter,tol=tol,total_iter_lambda=total_iter_lambda,m=m)

  return(LRres)
}
