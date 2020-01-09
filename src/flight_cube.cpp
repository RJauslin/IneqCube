#include <RcppArmadillo.h>
#include "onestepflightphase_arma.h"

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
  unsigned int J = X.n_cols;


  // arma::uvec i arma::regspace<arma::uvec>(0,  1,  J); // J+1 first units
  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find index of B
  // std::cout << i << std::endl;

  arma::mat B = A.rows(i);

  B = B.t();


 // arma::mat NN = arma::null(B.t());
 // arma::vec u;
 // int ncolNN = NN.n_cols;

 // while(ncolNN >= 1){
  // std::cout << ncolNN << std::endl;
 while(i.size() >= J+1){

   pik(i) = onestepflightphase_arma(B,pik(i));
   i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
   B = A.rows(i);
   B = B.t();


  //   // find u in kern of B
  //    u = NN.col(0);
  //
  //   // update pikR
  //   double l1 = 1e+200;
  //   double l2 = 1e+200;
  //   double l = 1e-9;
  //   for(unsigned int k = 0; k < i.size(); k++){
  //     if(u[k]> 0){
  //       l1 = std::min(l1,(1.0-pik[i[k]])/u[k]);
  //       l2 = std::min(l2,pik[i[k]]/u[k]);
  //     }
  //     if(u[k]< 0){
  //       l1 = std::min(l1,-pik[i[k]]/u[k]);
  //       l2 = std::min(l2,(pik[i[k]]-1.0)/u[k]);
  //     }
  //   }
  //   if(Rcpp::runif(1)[0]<l2/(l1+l2)){
  //     l = l1;
  //   }else{
  //     l = -l2;
  //   }
  //
  //   for(unsigned int k = 0; k < i.size(); k++){
  //     pik[i[k]] = pik[i[k]] + l*u[k];
  //     if(pik[i[k]] < EPS){
  //       pik[i[k]] = 0;
  //     }
  //     if(pik[i[k]] > (1-EPS)){
  //       pik[i[k]] = 1;
  //     }
  //   }
  //   // std::cout << i << std::endl;
  //
  //   i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
  //   // std::cout << i << std::endl;
  //   J = std::min(i.size(),X.n_cols);
  //   arma::mat B = A.rows(i);
  //   arma::mat NN = arma::null(B.t());
  //   ncolNN = NN.n_cols;
  }
  return(pik);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R


rm(list = ls())
N=5000
n=200
p=50
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))

system.time(test <- flightphase_arma(X,pik))
system.time(IneqCube::flightphase(pik,X))
system.time(BalancedSampling::flightphase(pik,X))

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
