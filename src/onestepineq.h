#ifndef onestepineq_H
#define onestepineq_H

#include <RcppArmadillo.h>

arma::vec onestepineq(arma::mat Ar,
                      arma::vec pik,
                      arma::mat Br,
                      arma::vec num,
                      double EPS=0.0000001);
#endif
