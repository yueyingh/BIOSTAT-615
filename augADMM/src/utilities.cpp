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

bool check_convergence(const double *x, const double *z, const double *y,
                       const double *x_prev, const double *z_prev, int length,
                       double abstol, double reltol) {
    double norm_x = 0.0, norm_z = 0.0, norm_y = 0.0;
    double norm_x_prev = 0.0, norm_z_prev = 0.0;

    // Calculate norms for x, z, y, and their previous values
    for (int i = 0; i < length; ++i) {
        norm_x += x[i] * x[i];
        norm_z += z[i] * z[i];
        norm_y += y[i] * y[i];
        norm_x_prev += x_prev[i] * x_prev[i];
        norm_z_prev += z_prev[i] * z_prev[i];
    }
    norm_x = sqrt(norm_x);
    norm_z = sqrt(norm_z);
    norm_y = sqrt(norm_y);
    norm_x_prev = sqrt(norm_x_prev);
    norm_z_prev = sqrt(norm_z_prev);

    // Compute primal and dual residuals
    double primal_res = 0.0, dual_res = 0.0;
    for (int i = 0; i < length; ++i) {
        primal_res += pow(x[i] - z[i], 2);
        dual_res += pow(x[i] - x_prev[i], 2);
    }
    primal_res = sqrt(primal_res);
    dual_res = sqrt(dual_res);

    // Compute tolerances
    double eps_pri = sqrt(length) * abstol + reltol * std::max(norm_x, norm_z);
    double eps_dual = sqrt(length) * abstol + reltol * norm_y;

    // Check if convergence criteria are met
    return primal_res <= eps_pri && dual_res <= eps_dual;
}
