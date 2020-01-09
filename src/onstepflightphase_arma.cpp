#include <RcppArmadillo.h>
#include "rref.h"

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
arma::vec onestepflightphase_arma(arma::mat B,arma::vec pik,double EPS=0.0000001){
  int ncol = B.n_cols;
  int nrow = B.n_rows;
  int i, j;
  arma::vec u(ncol,arma::fill::zeros);
  arma::uvec uset(ncol,arma::fill::zeros);
  // NumericVector u(ncol,0.0);
  // IntegerVector uset(ncol,0);
  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  rref(B);
  // std::cout << B << std::endl;
  for(i=(nrow-1);i>=0;i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(j=0;j<ncol;j++){
      if(B(i,j)==0.0){
        lead++;
      }else{
        break;
      }
    }
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(j=lead+1;j<ncol;j++){
        if( uset[j] == 0 ){
          uset[j] = 1;
          free *= -1.0;
          u[j] = free;
        }
        v -= u[j]*B(i,j);
      }
      u[lead] = v/B(i,lead);
      uset[lead] = 1;
    }
  }
  // unset u[i] are free and are set to 1 or -1, can only exist at beginning
  for(i=0;i<ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{
      break;
    }
  }
  // std::cout << u << std::endl;
  // find lambda1 and lambda2
  for(i=0;i<ncol;i++){
    if(u[i]>0){
      la1 = std::min(la1,(1-pik[i])/u[i]);
      la2 = std::min(la2,pik[i]/u[i]);
    }
    if(u[i]<0){
      la1 = std::min(la1,-pik[i]/u[i]);
      la2 = std::min(la2,(pik[i]-1)/u[i]);
    }
  }
  // random choice of p+lambda1*u and p-lambda2*u
  if(Rcpp::runif(1)[0]<la2/(la1+la2)){
    la = la1;
  }else{
    la = -la2;
  }
  // update pik
  for(i=0;i<ncol;i++){
    pik[i] = pik[i] + la * u[i];
    if(pik[i] < eps){
      pik[i] = 0;
    }
    if(pik[i] > 1-eps){
      pik[i] = 1;
    }
  }
  return pik;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R


rm(list = ls())
N=200
n=50
p=4
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))

x <- t(X[1:6,]/pik[1:6])

onestepflightphase_arma(x,pik[1:6])

IneqCube::flightphase(pik[1:6],X[1:6,])

*/
