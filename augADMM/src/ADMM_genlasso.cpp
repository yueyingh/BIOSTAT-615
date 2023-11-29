#include "ADMM_genlasso.h"
#include "dmatrix.h" // Includes matrix utility functions
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include "utilities.h" 

void ADMM_GenLasso(double *x, double *z, double *y, double *A, double *b,
                   const int length, const double lambda, const double rho,
                   const int max_iter, const double abstol,
                   const double reltol) {

  // Initialize necessary variables
  std::vector<double> AtA(length * length); // A^T * A
  std::vector<double> invATAplusRhoI(length *
                                     length); // (A^T * A + rho * I)^(-1)
  std::vector<double> Atb(length);            // A^T * b
  std::vector<double> x_prev(length),
      z_prev(length); // Store previous values for convergence check

  // Compute A^T * A and A^T * b
  dmat_B_ATA(length, length, A, AtA.data());
  dmat_yATx(length, length, A, b, Atb.data());

  // Compute (A^T * A + rho * I)^(-1)
  // This step may require a specialized matrix inversion function
  // Assuming it exists and is named dmat_inv
  for (int i = 0; i < length; ++i) {
    AtA[i * length + i] += rho;
  }
  dmat_inv(length, AtA.data(), invATAplusRhoI.data());

  for (int iter = 0; iter < max_iter; ++iter) {
    // Store previous values
    std::copy(x, x + length, x_prev.begin());
    std::copy(z, z + length, z_prev.begin());

    // Update x
    // x = (A^T * A + rho * I)^(-1) * (A^T * b + rho * (z - y))
    std::vector<double> tempZ(length);
    std::transform(z, z + length, y, tempZ.begin(), std::minus<double>());
    dmat_waxpby(length, rho, tempZ.data(), 1.0, Atb.data(), tempZ.data());
    dmat_yAx(length, length, invATAplusRhoI.data(), tempZ.data(), x);

    // Update z using soft thresholding
    // z = S_lambda(rho)(x + y)
    std::transform(x, x + length, y, tempZ.begin(), std::plus<double>());
    soft_threshold(tempZ.data(), z, lambda / rho, length);

    // Update y
    // y = y + rho * (x - z)
    for (int i = 0; i < length; ++i) {
      y[i] += rho * (x[i] - z[i]);
    }

    // Check for convergence
    if (check_convergence(x, z, y, x_prev.data(), z_prev.data(), length, abstol,
                          reltol)) {
      break;
    }
  }
}

// [[Rcpp::export]]
Rcpp::List ADMM_genlasso(Rcpp::NumericMatrix X, Rcpp::NumericVector Y,
                            Rcpp::NumericMatrix D, double lambda, double rho,
                            double abstol, double reltol, int maxiter) {

  // Convert Rcpp types to raw pointers for ADMM_GenLasso function
  int length = X.ncol();
  double *x = new double[length](); // Initialize to zero
  double *z = new double[length](); // Initialize to zero
  double *y = new double[length](); // Initialize to zero
  double *A = X.begin();
  double *b = Y.begin();

  // Call the ADMM_GenLasso function
  ADMM_GenLasso(x, z, y, A, b, length, lambda, rho, maxiter, abstol, reltol);

  // Convert results back to Rcpp types for returning to R
  Rcpp::NumericVector rcpp_x(x, x + length);
  Rcpp::NumericVector rcpp_z(z, z + length);
  Rcpp::NumericVector rcpp_y(y, y + length);

  // Clean up dynamically allocated memory
  delete[] x;
  delete[] z;
  delete[] y;

  // Return results as a list
  return Rcpp::List::create(Rcpp::Named("x") = rcpp_x,
                            Rcpp::Named("z") = rcpp_z,
                            Rcpp::Named("y") = rcpp_y);
}