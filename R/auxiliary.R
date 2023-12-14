# CHECKERS ----------------------------------------------------------------
#' @keywords internal
#' @noRd
check_data_matrix <- function(A) {
  cond1 <- (is.matrix(A)) # matrix
  cond2 <- (!(any(is.infinite(A)) || any(is.na(A))))
  if (cond1 && cond2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
check_data_vector <- function(b) {
  cond1 <- ((is.vector(b)) || ((is.matrix(b)) &&
    (length(b) == nrow(b)) || (length(b) == ncol(b))))
  cond2 <- (!(any(is.infinite(b)) || any(is.na(b))))
  if (cond1 && cond2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
check_param_constant <- function(num, lowerbound = 0) {
  cond1 <- (length(num) == 1)
  cond2 <- ((!is.infinite(num)) && (!is.na(num)))
  cond3 <- (num > lowerbound)

  if (cond1 && cond2 && cond3) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @keywords internal
#' @noRd
check_param_constant_multiple <- function(numvec, lowerbound = 0) {
  for (i in 1:length(numvec)) {
    if (!check_param_constant(numvec[i], lowerbound)) {
      return(FALSE)
    }
  }
  return(TRUE)
}


#' @keywords internal
#' @noRd
check_param_integer <- function(num, lowerbound = 1) {
    # Check if num is a single value and not NA or infinite
    cond1 <- (length(num) == 1) && (!is.infinite(num)) && (!is.na(num))

    # Check if num is greater than lowerbound
    cond2 <- (num >= lowerbound)

    # Check if num is an integer
    cond3 <- is.integer(num) || (num == floor(num))

    return(cond1 && cond2 && cond3)
}






# AUXILIARY COMPUTATIONS --------------------------------------------------
#   -----------------------------------------------------------------------
# 1. PseudoInverse using SVD and NumPy Scheme
# https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)
#' @keywords internal
#' @noRd
aux_pinv <- function(A) {
  svdA <- base::svd(A)
  tolerance <- (.Machine$double.eps) * max(c(nrow(A), ncol(A))) * as.double(max(svdA$d))

  idxcut <- which(svdA$d <= tolerance)
  invDvec <- (1 / svdA$d)
  invDvec[idxcut] <- 0

  output <- (svdA$v %*% diag(invDvec) %*% t(svdA$u))
  return(output)
}
#   -----------------------------------------------------------------------
