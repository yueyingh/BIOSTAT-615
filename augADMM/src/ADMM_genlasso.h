#ifndef ADMM_GENLASSO_H
#define ADMM_GENLASSO_H

void ADMM_GenLasso(double *x, double *z, double *y, double *A, double *b, 
                   const int length, const double lambda, const double rho, 
                   const int max_iter, const double abstol, const double reltol);

#endif // ADMM_GENLASSO_H
