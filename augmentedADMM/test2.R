###
# Testing with data generated from ADMM package
# Standard ADMM vs augADMM
###

library(augmentedADMM)
set.seed(123)

# Generate sample data
m <- 100
n <- 200
p <- 0.1 # Percentage of non-zero elements

x0 <- matrix(Matrix::rsparsematrix(n, 1, p))
A <- matrix(rnorm(m * n), nrow = m)
for (i in 1:ncol(A)) {
  A[, i] <- A[, i] / sqrt(sum(A[, i] * A[, i]))
}
b <- A %*% x0 + sqrt(0.001) * matrix(rnorm(m))
D <- diag(n)

# Set regularization lambda value
regval <- 0.1 * Matrix::norm(t(A) %*% b, "I")

# Solve LASSO via reducing from Generalized LASSO
output <- genlasso_admm(A, b, D, lambda = regval) # Set D as identity matrix
niter <- length(output$history$s_norm)
history <- output$history
output
pdf("convergence_plots_test2.pdf")

# Report convergence plot
opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(1:niter, history$objval, "b", main = "Cost Function")
plot(1:niter, history$r_norm, "b", main = "Primal Residual")
plot(1:niter, history$s_norm, "b", main = "Dual Residual")
par(opar)

M <- diag(ncol(A))

output_with_M <- genlasso_admm_aug(A, b, D, M, lambda = regval)
niter_with_M <- length(output_with_M$history$s_norm)
history_with_M <- output_with_M$history

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(1:niter_with_M, history_with_M$objval, "b", main = "Cost Function (with M)")
plot(1:niter_with_M, history_with_M$r_norm, "b", main = "Primal Residual (with M)")
plot(1:niter_with_M, history_with_M$s_norm, "b", main = "Dual Residual (with M)")
par(opar)

dev.off()


###
# test from ADMM package
###
library(ADMM)
set.seed(123)
m = 100
n = 200
p = 0.1   # percentange of non-zero elements

x0 = matrix(Matrix::rsparsematrix(n,1,p))
A  = matrix(rnorm(m*n),nrow=m)
for (i in 1:ncol(A)){
  A[,i] = A[,i]/sqrt(sum(A[,i]*A[,i]))
}
b = A%*%x0 + sqrt(0.001)*matrix(rnorm(m))
D = diag(n);
 ## set regularization lambda value
regval = 0.1*Matrix::norm(t(A)%*%b, 'I')

## solve LASSO via reducing from Generalized LASSO
output_admm  = admm.genlasso(A,b,D,lambda=regval) # set D as identity matrix
niter   = length(output$history$s_norm)
history = output$history

output_admm

 ## report convergence plot
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3))
plot(1:niter, history$objval, "b", main="cost function")
plot(1:niter, history$r_norm, "b", main="primal residual")
plot(1:niter, history$s_norm, "b", main="dual residual")
par(opar)


