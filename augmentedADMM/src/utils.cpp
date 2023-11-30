// utils.cpp
#include "utils.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// admm_genlasso ======================================================
arma::colvec genlasso_shrinkage(arma::colvec a, const double kappa) {
  const int n = a.n_elem;
  arma::colvec y(n, fill::zeros);
  for (int i = 0; i < n; i++) {
    // first term : max(0, a-kappa)
    if (a(i) - kappa > 0) {
      y(i) = a(i) - kappa;
    }
    // second term : -max(0, -a-kappa)
    if (-a(i) - kappa > 0) {
      y(i) = y(i) + a(i) + kappa;
    }
  }
  return (y);
}

double genlasso_objective(const arma::mat &A, const arma::colvec &b,
                          const arma::mat &D, const double lambda,
                          const arma::colvec &x, const arma::colvec &z) {
  return (pow(norm(A * x - b, 2), 2) / 2 + lambda * norm(D * x, 1));
}

arma::mat genlasso_factor(const arma::mat &A, double rho, const arma::mat &D) {
  // const int m = A.n_rows;
  // const int n = A.n_cols;
  arma::mat U;
  U = chol(A.t() * A + rho * D.t() * D);
  return (U);
}