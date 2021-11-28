# R package `MedLRT`

## Overview 
MedLRT is a R package for implementing the likelihood ratio test for making statistical inferences on the joint mediation effects of multiple mediators

## Install from GitHub
Please make sure you have installed R package `mvtnorm` and `devtools`. If not, you can install it using the following command:

```r
install.packages("mvtnorm")
install.packages("devtools")
```
Install `MedLRT` from GitHub with 

```r
devtools::install_github("haowei72/MedLRT")
```

## Documentaiton

You can check the help document for the main function

```r
?MedLRT
```

## Examples

```r
n = 100
q = 10
L = 3
X = rnorm(n)
Z = cbind(1,matrix(rnorm(n*L),nrow=n,ncol=L-1))
M = X%*%t(rep(c(1,0,1),length=q)) + 0.1*matrix(rnorm(n*q),nrow=n,ncol=q)
Y = X + M%*%rep(c(,1,0),length=q) + Z%*%matrix(1,nrow=L,ncol=1) + 0.1*rnorm(n)
res <- MedLRT(X=X,M=M,Y=Y,Z=Z)
print(as.data.frame(res))
```
