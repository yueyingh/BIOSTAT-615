#ifndef UTILITIES_H
#define UTILITIES_H

void dmat_B_ATA(int rows, int cols, const double* A, double* result);
void dmat_yATx(int rows, int cols, const double* A, const double* b, double* result);
void dmat_inv(int size, const double* A, double* invA);
void dmat_waxpby(int size, double alpha, const double* x, double beta, const double* y, double* result);
void dmat_yAx(int rows, int cols, const double* A, const double* x, double* y);
void soft_threshold(double* x, double* z, double lambda, int length);
bool check_convergence(const double *x, const double *z, const double *y,
                       const double *x_prev, const double *z_prev, int length,
                       double abstol, double reltol);

#endif // UTILITIES_H
