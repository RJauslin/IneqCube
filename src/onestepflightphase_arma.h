#ifndef onestepflightphase_arma_H
#define onestepflightphase_arma_H

#include <RcppArmadillo.h>

arma::vec onestepflightphase_arma(arma::mat B,arma::vec pik,double EPS=0.0000001);

#endif
