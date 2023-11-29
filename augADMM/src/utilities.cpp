#include "utilities.h"
#include <algorithm> // For std::transform
#include <cmath>     // For mathematical functions

// Assuming you have these functions declared in utilities.h
// soft_threshold, check_convergence, sign, dmat_norm2, dmat_norm1

void soft_threshold(double *data, double *result, double lambda, int length) {
  for (int i = 0; i < length; ++i) {
    result[i] = (data[i] > lambda)    ? data[i] - lambda
                : (data[i] < -lambda) ? data[i] + lambda
                                      : 0.0;
  }
}

bool check_convergence(const std::vector<double> &x,
                       const std::vector<double> &z,
                       const std::vector<double> &y,
                       const std::vector<double> &x_prev,
                       const std::vector<double> &z_prev, int length,
                       double abstol, double reltol) {
  // Implement the convergence check logic
  // This is an example and might need adjustments

  double norm_x = dmat_norm2(x);
  double norm_z = dmat_norm2(z);
  double norm_y = dmat_norm2(y);

  for (int i = 0; i < length; ++i) {
    if (std::fabs(x[i] - x_prev[i]) > (abstol + reltol * norm_x) ||
        std::fabs(z[i] - z_prev[i]) > (abstol + reltol * norm_z) ||
        std::fabs(y[i] - y_prev[i]) > (abstol + reltol * norm_y)) {
      return false;
    }
  }
  return true;
}

double sign(double value) { return (value > 0) - (value < 0); }

double dmat_norm2(const std::vector<double> &vec) {
  double sum = 0.0;
  for (double val : vec) {
    sum += val * val;
  }
  return std::sqrt(sum);
}

double dmat_norm1(const std::vector<double> &vec) {
  double sum = 0.0;
  for (double val : vec) {
    sum += std::fabs(val);
  }
  return sum;
}
