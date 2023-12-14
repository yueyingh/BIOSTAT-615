#' Augmented ADMM Solver for Generalized Lasso Problems
#'
#' Solves the generalized lasso problem using the augmented Alternating Direction Method of Multipliers (ADMM). 
#' This function extends the standard ADMM approach by incorporating an additional matrix 'M', which can lead to more efficient computations in certain scenarios.
#' 
#' @param A Data matrix.
#' @param b Response vector.
#' @param D Regularization matrix, defaults to the identity matrix of the length of `b`.
#' @param M Augmentation matrix used in the augmented ADMM.
#' @param lambda Regularization parameter.
#' @param rho Augmentation parameter for ADMM.
#' @param abstol Absolute tolerance level for convergence.
#' @param reltol Relative tolerance level for convergence.
#' @param maxiter Maximum number of iterations for the algorithm.
#'
#' @return A list containing the following components:
#'   - `x`: The solution vector.
#'   - `x_iter`: Solution vector at each iteration.
#'   - `history`: A data frame containing the objective value, primal and dual norms, 
#'     and epsilon values for primal and dual convergence checks at each iteration.
#'
#' @examples
#' # Generate synthetic data
#' n <- 200; m <- 100; p <- 0.1
#' A <- matrix(rnorm(m * n), nrow = m)
#' b <- matrix(rnorm(m))
#' D <- diag(n)
#' M <- diag(n) # Example augmentation matrix
#'
#' # Solve using aug_admm_genlasso
#' result <- aug_admm_genlasso(A, b, D, M)
#'
#' @export
aug_admm_genlasso <- function(A, b, D = diag(length(b)), M, lambda = 1.0, rho = 1.0,
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
