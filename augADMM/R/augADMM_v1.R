#' @title Augmented ADMM for Generalized LASSO
#' @description Solves generalized LASSO problems using augmented ADMM.
#' @param A Regressor matrix.
#' @param b Response vector.
#' @param D Regularization matrix.
#' @param lambda Regularization parameter.
#' @param rho Augmented Lagrangian parameter.
#' @param alpha Overrelaxation parameter.
#' @param abstol Absolute tolerance for convergence.
#' @param reltol Relative tolerance for convergence.
#' @param maxiter Maximum number of iterations.
#' @return A list containing the solution and iteration history.
#' @export
augADMM_genlasso <- function(A, b, D, lambda, rho, alpha, abstol, reltol, maxiter) {
  # Call the C++ function
  result <- augADMM_genlasso(as.matrix(A), as.vector(b), as.matrix(D), lambda, rho, alpha, abstol, reltol, maxiter)
  
  # Process and return the result
  return(result)
}


