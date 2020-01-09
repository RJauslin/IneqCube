#include <RcppArmadillo.h>

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
void rref(arma::mat& M){

  int lead = 0;
  int n = M.n_rows;
  int p = M.n_cols;

  double eps = 1e-11;

  int r,i,k;
  double temp;

  for(r = 0; r < n; r++){// loop on the row
    if(p<=lead){
      return;
    }
    i = r;
    while( std::max(M(i,lead),-M(i,lead)) < eps ){

      M(i,lead) = 0.0;
      i = i + 1;
      if(i == n){
        i = r;
        lead = lead + 1;
        if(p == lead){
          return;
        }
      }
    }
    // swap rows i and r
    for(k = 0; k < p; k++){
      temp = M(i,k);
      M(i,k) = M(r,k);
      M(r,k) = temp;
    }
    // If M(r, lead) is not 0 divide row r by M(r, lead)
    if( M(r,lead) != 0.0 ){
      temp = M(r,lead);
      for(k=0;k<lead;k++){
        M(r,k) = 0.0;
      }
      for(k=lead;k<p;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(i = 0;i < n;i++){
      if( i != r ){
        // Subtract M(i, lead) multiplied by row r from row i
        temp = M(i,lead);
        for(k=0;k<p;k++){
          M(i,k) = M(i,k) - temp * M(r,k);
        }
      }
    }
    lead = lead + 1;
  }
  return;
}


/*** R

A <- matrix(c(1,3,-1,0,1,7),c(2,3),byrow = T)
rref(A)

*/


