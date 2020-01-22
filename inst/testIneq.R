
rm(list = ls())
comment = TRUE
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
# pikstar <- as.vector(ineq(X,pik,B,r,EPS=0.000001))
# t(B)%*%pikstar
# r
#
#
# s <- landIneq(X,pikstar,pik,B,r)
# s <- sampling::landingcube(X,pikstar,pik)
# t(B)%*%s -r
#
# t(X/pik)%*%s
# t(X/pik)%*%pik
#
# test <- which(t(B)%*%s > r)
# r[test]
# as.vector(t(B)%*%s)[test]



# ON COMMENCE PAR EQUILIBRER TOUTEs
# LES VARIABLES (Y COMPRIS LES INEGALITE)
# A LA FIN DE LA PHASE DE VOL ON A DONC AU MAXIUM
# p+q+1 UNITEES NON-EGALE A 0 OU 1

piks <- as.vector(flightphase_arma(cbind(X,B*pik),pik))
p+q+1 # + 1 pour la colonne 1
length(which(piks > EPS & piks <(1-EPS)))
# CHECK QUE TOUT SOIT SATISFAIT
round(t(B)%*%piks - r,2)
t(X/pik)%*%piks
t(X/pik)%*%pik

# ON RELACHE LES CONTRAINTES D'INEGALITE
# IL RESTE AU MAXIMUM p+1 UNITE
pikstar <- ineq(X/pik*piks,piks,B,r)
# pikstar <- ineq(X/pik,pik,B,r)
# pikstar <- fast.flight.cube.ineq(X/pik*piks,piks,B,r)
p+1 # + 1 pour la colonne 1
length(which(pikstar > EPS & pikstar <(1-EPS)))
# check -> CERTAINE INEGALITE NE SONT PLUS SATISFAITE
round(t(B)%*%pikstar - r,3)
t(X/pik)%*%pikstar
t(X/pik)%*%pik


# LAND DES p+1 UNITE SUR LES ECHANTILLONS QUI SATISFONT LES CONTRAINTES RESTANTES.
s <- landIneq(X,pikstar,pik,B,r)

round(t(B)%*%s - r,3)
t(X/pik)%*%s
t(X/pik)%*%pik


#####################################################################################
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
s=as.vector(fast.flight.cube.ineq(X,pik,B,r,EPS=0.000001))
t(B)%*%pik <= r
t(B)%*%s <= r

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

s=landIneq(cbind(pik,B*pik),pikstar,pik,comment=TRUE)


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









