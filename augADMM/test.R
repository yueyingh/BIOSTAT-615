library(ADMM)
library(augADMM)
augADMMExample(1, 3)

# Generate sample data
set.seed(123) # for reproducibility
m <- 100 # Number of observations
n <- 200 # Number of features
p <- 0.1 # Sparsity of the true coefficient

# Generate a sparse solution vector
x_true <- matrix(rnorm(n) * (runif(n) < p))

# Create a regressor matrix A
A <- matrix(rnorm(m * n), nrow = m)
for (i in 1:ncol(A)) {
    A[, i] <- A[, i] / sqrt(sum(A[, i] * A[, i]))
}

# Generate the response vector b
b <- A %*% x_true + rnorm(m) * 0.01


D <- diag(n)

# Set regularization parameter lambda
lambda <- 0.1 * max(abs(t(A) %*% b))
# Solve the generalized LASSO problem
result <- admm.genlasso(A, b, D, lambda = lambda)

# Extract the solution
solution <- result$x

# Return a summary of the solution
summary(solution)

summary(x_true)


# Coefficient Comparison:
# You can plot the true coefficients (x_true) against the estimated coefficients (solution) to visually inspect how closely they match.


plot(x_true, solution, xlab = "True Coefficients", ylab = "Estimated Coefficients")
abline(0, 1, col = "red") # Line for perfect agreement

true_nonzero_indices <- which(x_true != 0)
estimated_nonzero_indices <- which(solution != 0)
# Check overlap
common_indices <- intersect(true_nonzero_indices, estimated_nonzero_indices)
length(common_indices)

# Convergence Analysis
plot(result$history$objval, type = "l", xlab = "Iteration", ylab = "Objective Function Value")

# Primal and Dual Residuals:
# The primal and dual residuals indicate how close the current solution is to satisfying the constraints of the optimization problem. Their trends over iterations can help assess convergence.
plot(result$history$r_norm, type = "l", col = "blue", xlab = "Iteration", ylab = "Norm of Residuals")
lines(result$history$s_norm, col = "red")
legend("topright", legend = c("Primal Residual", "Dual Residual"), col = c("blue", "red"), lty = 1)

# Feasibility Tolerances:
# Comparing primal and dual residuals to their respective feasibility tolerances (eps_pri and eps_dual) can help you understand when the algorithm has effectively converged.
plot(result$history$r_norm, type = "l", col = "blue", xlab = "Iteration", ylab = "Norm of Residuals")
lines(result$history$eps_pri, col = "green", lty = 2)
lines(result$history$s_norm, col = "red")
lines(result$history$eps_dual, col = "orange", lty = 2)
legend("topright", legend = c("Primal Residual", "Dual Residual", "Eps Pri", "Eps Dual"), col = c("blue", "red", "green", "orange"), lty = 1:2)


# 3. Other Metrics
# Runtime: You can measure the time taken by the algorithm to solve the problem. This is particularly useful when comparing the performance of different algorithms.
start_time <- Sys.time()
# Call to admm.genlasso
end_time <- Sys.time()
runtime <- end_time - start_time
print(runtime)


