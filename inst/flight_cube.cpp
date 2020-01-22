#include <RcppArmadillo.h>
#include "onestep.h"
#include "onestepineq.h"



/*** R

rm(list = ls())
library(sampling)
library(MASS)
EPS=0.0000001
N=1000
n=300
p=2
q=7
z=runif(N)
#z=rep(1,N)
pik=inclusionprobabilities(z,n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
r=c(ceiling(pik%*%B))
r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
s=ineq(X,pik,B,r,EPS=0.000001)
round(s,3)
# c(t(B)%*%pik)
# c(r)
# c(t(B)%*%s)
# colSums(X)
# c(t(X)%*%(s/pik))

piks=as.vector(flightphase_arma2(as.matrix(cbind(X,B*pik)),pik))
pikstar=ineq(as.matrix(X/pik*piks),piks,B,r,EPS=0.000001)
sum(pikstar>EPS & pikstar<1-EPS )

piks=flightphase_arma2(as.matrix(cbind(pik,B*pik)),pik,EPS=0.000001)
pikstar=ineq(piks,piks,B,r,EPS=0.000001)
sum(pikstar>EPS & pikstar<1-EPS )

s=landingcube(cbind(pik,B*pik),pikstar,pik,comment=TRUE)


rm(list = ls())
EPS=0.0000001
N=500
n=50
p=4
q=7
z=runif(N)
#z=rep(1,N)
pik=sampling::inclusionprobabilities(z,n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
r=c(ceiling(pik%*%B))
r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
system.time(s <- ineq(X,pik,B,r,EPS=0.000001))
system.time(s <- fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.000001))



t(B)%*%s <= r
t(A)%*%s
t(A)%*%pik





s <- flightphase_arma2(X,pik)
t(A)%*%s
t(A)%*%pik



round(s,3)
c(t(B)%*%pik)
c(r)
round(c(t(B)%*%pik)-c(t(B)%*%s),3)
c(t(B)%*%pikstar)-c(t(B)%*%s)
colSums(X)
c(t(X)%*%(s/pik))


abs(c(t(B)%*%pik)-r)<=EPS

rm(list = ls())
N=5000
n=200
p=20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))

system.time(test <- flightphase_arma2(X,pik))
system.time(IneqCube::flightphase(pik,X))
system.time(BalancedSampling::flightphase(pik,X))

(X/pik)%%


x <- t(X[1:6,]/pik[1:6])

onestepflightphase_arma(x,pik[1:6])



rm(list = ls())
N=200
n=50
p=4
pik=inclusionprobabilities(runif(N),n)

X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
# X=cbind(pik)
# X = c()



system.time(piks <- fast.flight.cube(X,pik,EPS=0.000001))
system.time(s <- IneqCube:::flightphase_arma(X,pik))
system.time(sBal <- flightphase(pik,X))



rm(list = ls())
library(sampling)
library(MASS)
EPS=0.0000001
N=8000
n=500
p=1
q=7
z=runif(N)
#z=rep(1,N)
pik=inclusionprobabilities(z,n)
# pik[3] <- 0
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
B=cbind(Z,-Z)
r=c(ceiling(pik%*%B))
r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
# s=fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.000001)
# round(s,3)
# c(t(B)%*%pik)
# c(r)
# c(t(B)%*%s)
# colSums(X)
# c(t(X)%*%(s/pik))
#

dim(cbind(X,B*pik))

X <- cbind(X,B*pik)
A <- X/pik
dim(A)
test <- MASS::Null(A)

system.time(piks <- fast.flight.cube(X,pik,deepness=1,EPS=0.000001))

system.time(test <- IneqCube:::flightphase_arma(X,pik))


*/
