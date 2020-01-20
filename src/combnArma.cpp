#include <RcppArmadillo.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <functional>
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
long long int chooseArma(int n, int k) {
  if (k == 0 || k == n)
    return 1;
  return chooseArma(n - 1, k - 1) + chooseArma(n - 1, k);
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
long long int choose(int n, int k)
{
  long long int res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k )
    k = n - k;

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description Rosetta code combn
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
arma::mat comb(int N, int K){

  arma::uword count = choose(N,K);
  arma::mat out(N,count,arma::fill::zeros);
  arma::uvec tmp(K,arma::fill::zeros);
  arma::vec tmp2(N,arma::fill::zeros);

  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's

  // print integers and permute bitmask
  int c = 0;
  do {
    int c1 = 0;
    for(int i = 0; i < N; ++i){
      if(bitmask[i]){
        tmp(c1) = i;
        c1++;
      }
    }
    tmp2(tmp) += 1;
    out.col(c) = tmp2;
    tmp2(tmp) -= 1;
    c++;
  } while(std::prev_permutation(bitmask.begin(), bitmask.end()));
  return(out);
}

/*** R

system.time(test <- IneqCube::choose(10,3))
system.time(test <- base::choose(10,3))


system.time(test <- comb(25,10))
system.time(test2 <- combn(25,10))
system.time(test <- sampling::writesample(10,25))
test
test2

x <- seq(1,5,1)
m = 3
combnArma(x,m)

# system.time(combnArma(50,10))
# choose(50,10)
*/
