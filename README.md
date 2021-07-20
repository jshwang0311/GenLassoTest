# ESPgenlasso (Exact Solution Path of generalized lasso problem for rank-deficient design matrices)
## Description
This R package provides a solution path of generalized lasso problem for rank-deficient design matrices, which is an exact entrie solution path. And also when the Lasso penalty is contained in the penalty matrix, the specific solutions, which are the (partial) minimal (or/and) maximal $l_1$ norm solution and the (partial) $l_{\infty}$ norm solution, are obtained. Details are described in the following paper:
- Jaesung Hwang, Joongho Won, Yongdai Kim, Characterization of the Solutions Set of the Generalized Lasso for Rank-deficient Design Matrices.

## Install
In R terminal,

```R
install.packages("devtools")
devtools::install_github("jshwang0311/ESPgenlasso")
```

## Example
Import library, generate simulated data and make a penalty matrix using the following commands.
```R
library(ESPgenlasso)
library(mvtnorm)

n <- 100
p <- 500

X <- rmvnorm(n,sigma = diag(1,p))
beta.coef <- rep(0,p)
beta.coef[c((p/10+1):(2*p/10))] <- 1
beta.coef[c((2*p/10+1):(4*p/10))] <- 2
beta.coef <- matrix(beta.coef,ncol=1)
err.val <- matrix(rmvnorm(1,sigma = diag(1,n)),ncol=1)
y <- X%*%beta.coef + err.val

D <- matrix(0,ncol=p,nrow=200)
for(i in c(1:200)){
  D[i,(299+i)] <- -1
  D[i,(i+300)] <- 1
}
eps <- 0.1
temp.D <- diag(eps,p)
D <- rbind(D,temp.D)
```

Find the exact solution path algorithm with the following command.
```R
new.prop.rslt <- ESPgenlasso::ESPgenlasso(y,X,D,iter = 50000)
```

To find the (partial) minimal (or/and) maximal $l_1$ norm solution for the given lambdas, enter the following command.

