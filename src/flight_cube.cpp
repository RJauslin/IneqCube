#include <RcppArmadillo.h>
#include "onestep.h"
#include "onestepineq.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param x x
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' func
//'
//' @examples
//'
//' @export
// [[Rcpp::export]]
arma::vec flightphase_arma(arma::mat X,arma::vec pik,double EPS=0.0000001){

  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B

  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat B = (A.rows(i)).t(); // extract B of A

  while(i.size() > 0){
    // std::cout << i.size() << std::endl;
    pik.elem(i) =  onestep(B,pik.elem(i));
    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    B = (A.rows(i)).t();
    // std::cout << B << std::endl;
    if(i.size() < (J+1)){
      arma::mat kern = arma::null(B);
      if(kern.empty()){
        break;
      }
    }
  }
  return(pik);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param x x
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' func
//'
//' @examples
//'
//' @export
// [[Rcpp::export]]
arma::vec ineq(arma::mat X,
               arma::vec pik,
               arma::mat B,
               arma::vec r,
               double EPS=0.0000001){

  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B

  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat Ar = (A.rows(i)).t(); // extract B of A
  arma::mat Br = (B.rows(i)).t();

  arma::vec c = B.t()*pik; // initialize c
  arma::vec num = r-c; // initizalize r-c

  while(i.size() > 0){

    pik(i) = onestepineq(Ar,pik(i),Br,num);

    arma::vec cond = arma::abs(B.t()*pik-r);
    r = r(find(cond >= EPS));
    c = c(find(cond >= EPS));
    num = r-c;
    for(arma::uword k = 0; k < cond.size(); k++){
      if( cond[k] <= EPS){
        A = join_rows( A, B.col(k));
        B.shed_col(k);
        J = A.n_cols;
      }
    }

    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    Ar = (A.rows(i)).t();
    Br = (B.rows(i)).t();

    if(i.size() < (J+1)){
      arma::mat kern = arma::null(Ar);
      if(kern.empty()){
        break;
      }
    }
  }
  return(pik);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

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
