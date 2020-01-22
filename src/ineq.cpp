#include <RcppArmadillo.h>
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
arma::vec ineq(arma::mat X,
               arma::vec pik,
               arma::mat B,
               arma::vec r,
               double EPS=0.0000001){

  // INITIALIZING A AND J
  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B
  arma::uvec rest = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
  /*
   * MAIN LOOP
   */
  do{

    // UPDATE A AND B IF BOUND IS REACHED
    arma::vec cond = arma::abs(r - B.t()*pik);

    arma::uvec s_cols = find(cond <= EPS);
    A = arma::join_rows(A, B.cols(s_cols));
    B.shed_cols(s_cols);
    J = A.n_cols;


    // DEFINE rest, Arest and Brest
    rest = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    arma::mat Ar = (A.rows(rest)).t();
    arma::mat Br = (B.rows(rest)).t();

    // IF REST IS SMALL THEN NULLSPACE IS CHECKED
    if(rest.size() < (J+1)){
      arma::mat kern = arma::null(Ar);
      if(kern.empty()){
        break;
      }
    }

    // UPDATE c,r and NUMERATOR THAT IS PASSED TO onestepineq
    arma::vec c = B.t()*pik;
    r = r(find(cond >= EPS));
    arma::vec num = r-c;



    // ONE STEP THAT RESPECT BOTH EQUALITY AND INEQUALITY
    pik(rest) = onestepineq(Ar,pik(rest),Br,num);

  } while (rest.size() > 0);

  return(pik);
}
