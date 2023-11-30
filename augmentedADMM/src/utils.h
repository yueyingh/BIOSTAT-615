// utils.h
#include <RcppArmadillo.h>
#ifndef UTILS_H
#define UTILS_H

arma::colvec genlasso_shrinkage(arma::colvec a, const double kappa);
double genlasso_objective(const arma::mat &A, const arma::colvec &b, const arma::mat &D, const double lambda, const arma::colvec &x, const arma::colvec &z);
arma::mat genlasso_factor(const arma::mat &A, double rho, const arma::mat &D);

#endif