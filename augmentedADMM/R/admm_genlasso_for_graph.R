#' @export
admm_genlasso_for_graph <- function(A, b, D = diag(length(b)), C, lambda1 = 1.0, lambda2 = 1.0, rho = 1.0,
                                    abstol = 1e-4, reltol = 1e-2, maxiter = 1000) {
  #-----------------------------------------------------------
  ## PREPROCESSING
  # 1. data validity
  if (!check_data_matrix(A)) {
    stop("* ADMM.GENLASSO : input 'A' is invalid data matrix.")
  }
  if (!check_data_vector(b)) {
    stop("* ADMM.GENLASSO : input 'b' is invalid data vector")
  }
  b <- as.vector(b)
  # 2. data size
  if (nrow(A) != length(b)) {
    stop("* ADMM.GENLASSO : two inputs 'A' and 'b' have non-matching dimension.")
  }
  # 3. D : regularization matrix
  if (!check_data_matrix(D)) {
    stop("* ADMM.GENLASSO : input 'D' is invalid regularization matrix.")
  }
  if (ncol(A) != ncol(D)) {
    stop("* ADMM.GENLASSO : input 'D' has invalid size.")
  }
  if (!check_data_matrix(C)) {
    stop("* ADMM.GENLASSO : input 'C' is invalid data matrix.")
  }
  # 4. other parameters
  if (!check_param_constant_multiple(c(abstol, reltol))) {
    stop("* ADMM.GENLASSO : tolerance level is invalid.")
  }
  if (!check_param_integer(maxiter, 2)) {
    stop("* ADMM.GENLASSO : 'maxiter' should be a positive integer.")
  }
  maxiter <- as.integer(maxiter)
  rho <- as.double(rho)
  if (!check_param_constant(rho, 0)) {
    stop("* ADMM.GENLASSO : 'rho' should be a positive real number.")
  }

  #-----------------------------------------------------------
  ## MAIN COMPUTATION
  #   1. lambda=0 case; pseudoinverse
  meps <- (.Machine$double.eps)
  negsmall <- -meps
  lambda1 <- as.double(lambda1)
  lambda2 <- as.double(lambda2)
  if (!check_param_constant(lambda1, negsmall)) {
    stop("* ADMM.GENLASSO : 'lambda1' is invalid; should be a nonnegative real number.")
  }
  if (!check_param_constant(lambda2, negsmall)) {
    stop("* ADMM.GENLASSO : 'lambda2' is invalid; should be a nonnegative real number.")
  }
  if (lambda1 < meps & lambda2 < meps) {
    message("* ADMM.GENLASSO : since both regularization parameters are effectively zero, a least-squares solution is returned.")
    xsol <- as.vector(aux_pinv(A) %*% matrix(b))
    output <- list()
    output$x <- xsol
    return(output)
  }

  #   2. main computation : Xiaozhi's work
  result <- admm_genlasso_for_graph_CPP(A, b, D, C, lambda1, lambda2, reltol, abstol, maxiter, rho)


  #-----------------------------------------------------------
  ## RESULT RETURN
  kk <- result$k
  output <- list()
  output$x <- result$x
  output$history <- data.frame(
    objval = result$objval[1:kk],
    r_norm = result$r_norm[1:kk],
    s_norm = result$s_norm[1:kk],
    eps_pri = result$eps_pri[1:kk],
    eps_dual = result$eps_dual[1:kk]
  )
  return(output)
}
