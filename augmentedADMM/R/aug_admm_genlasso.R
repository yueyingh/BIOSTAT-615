#' @export
aug_admm_genlasso <- function(A, b, D = diag(length(b)), M, lambda = 1.0, rho = 1.0, alpha = 1.0,
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
    if (!check_data_matrix(M)) {
        stop("* ADMM.GENLASSO : input 'M' is invalid data matrix.")
    }
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
    lambda <- as.double(lambda)
    if (!check_param_constant(lambda, negsmall)) {
        stop("* ADMM.GENLASSO : 'lambda' is invalid; should be a nonnegative real number.")
    }
    if (lambda < meps) {
        message("* ADMM.GENLASSO : since both regularization parameters are effectively zero, a least-squares solution is returned.")
        xsol <- as.vector(aux_pinv(A) %*% matrix(b))
        output <- list()
        output$x <- xsol
        return(output)
    }

    #   2. main computation : Xiaozhi's work
    result <- aug_admm_genlasso_CPP(A, b, D, M, lambda, reltol, abstol, maxiter, rho)


    #-----------------------------------------------------------
    ## RESULT RETURN
    kk <- result$k
    output <- list()
    output$x <- result$x
    output$x_iter <- result$x_iter
    output$history <- data.frame(
        objval = result$objval[1:kk],
        r_norm = result$r_norm[1:kk],
        s_norm = result$s_norm[1:kk],
        eps_pri = result$eps_pri[1:kk],
        eps_dual = result$eps_dual[1:kk]
    )
    return(output)
}
