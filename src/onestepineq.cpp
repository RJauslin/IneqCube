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
arma::vec onestepineq(arma::mat Ar,
                      arma::vec pik,
                      arma::mat Br,
                      arma::vec num,
                      double EPS=0.0000001){

  arma::mat kern = arma::null(Ar);
  unsigned int N = pik.size();
  arma::vec u(N);
  u = kern.col(0);


  double l1 = 1e+200;
  double l2 = 1e+200;
  double mu1 = 1e+200;
  double mu2 = 1e+200;
  double l = 1e-9;
  for(arma::uword k = 0; k < N; k++){
    if(u[k]> 0){
      mu1 = std::min(mu1,(1.0 - pik[k])/u[k]);
      mu2 = std::min(mu2,pik[k]/u[k]);
    }
    if(u[k]< 0){
      mu1 = std::min(mu1,-pik[k]/u[k]);
      mu2 = std::min(mu2,(pik[k]-1.0)/u[k]);
    }
  }

  arma::vec pet = num;
  arma::vec tmp = Br*u;
  // std::cout << tmp << std::endl;
  // std::cout << pet << std::endl;
  pet = pet/tmp;
  // std::cout << pet << std::endl;
  double nu1 = 1e+200;
  double nu2 = 1e+200;
  for(arma::uword k = 0; k < pet.size(); k++){
    if(pet[k]> 0){
      nu1 = std::min(nu1,pet[k]);
    }
    if(pet[k]< 0){
      nu2 = std::min(nu2,-pet[k]);
    }
  }

  l1 = std::min(mu1,nu1);
  l2 = std::min(mu2,nu2);
  if(Rcpp::runif(1)[0]<l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  for(arma::uword k = 0; k < N; k++){
    pik[k] = pik[k] + l*u[k];
    if(pik[k] < EPS){
      pik[k] = 0;
    }
    if(pik[k] > (1-EPS)){
      pik[k] = 1;
    }
  }
 return(pik);
}
