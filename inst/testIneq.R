
library(sampling)
library(MASS)
EPS=0.0000001
N=1000
n=300
p=1
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
s=fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.000001)
round(s,3)
c(t(B)%*%pik)
c(r)
c(t(B)%*%s)
colSums(X)
c(t(X)%*%(s/pik))

piks=fast.flight.cube(cbind(X,B*pik),pik,deepness=1,EPS=0.000001)
pikstar=fast.flight.cube.ineq(X/pik*piks,piks,B,r,deepness=1,EPS=0.000001)
sum(pikstar>EPS & pikstar<1-EPS )

piks=fast.flight.cube(cbind(pik,B*pik),pik,deepness=1,EPS=0.000001)
pikstar=fast.flight.cube.ineq(piks,piks,B,r,deepness=1,EPS=0.000001)
sum(pikstar>EPS & pikstar<1-EPS )

s=landingcube(cbind(pik,B*pik),pikstar,pik,comment=TRUE)


round(s,3)
c(t(B)%*%pik)
c(r)
round(c(t(B)%*%pik)-c(t(B)%*%s),3)
c(t(B)%*%pikstar)-c(t(B)%*%s)
colSums(X)
c(t(X)%*%(s/pik))


a=1;j=0
SIM=100
PIK=rep(0,N)
while(a>=-0.00001 & j<SIM)
{
  j=j+1
  print(j)
  piks=fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.0000001)
  PIK=PIK+piks
  a=sum(r-c(t(B)%*%piks))
}
j
round((r-c(t(B)%*%piks))/(c(t(B)%*%pik)),4)
pike=PIK/SIM
sum(abs((pike-pik)/sqrt(pik*(1-pik))*sqrt(SIM))<1.96)/N








