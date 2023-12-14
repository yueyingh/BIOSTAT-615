###
# Script Purpose: Compare the precision of standard ADMM and augmented ADMM
# on a generated sparse matrix problem. The script outputs a plot 
# comparing the Mean Squared Error (MSE) at each iteration for both methods.
###
 
setwd("augmentedADMM")

library(ggplot2)
library(Matrix)
library(augmentedADMM)
library(reshape2)

set.seed(123)

# Define dimensions of the problem for demonstration purposes
n <- 200 # number of rows
m <- 100 # number of columns
p <- 0.1 # sparsity level

# Generate a sparse matrix x0 as the true value
x0 <- matrix(Matrix::rsparsematrix(n, 1, p))

# Generate a random matrix A and normalize its columns
A <- matrix(rnorm(m * n), nrow = m)
for (i in 1:ncol(A)) {
    A[, i] <- A[, i] / sqrt(sum(A[, i]^2))
}

# Generate response vector b with some noise
b <- A %*% x0 + sqrt(0.001) * matrix(rnorm(m))

# Define diagonal matrix D and matrix M for augmented ADMM
D <- diag(n)
M <- norm(D, "2") * diag(n)

# Set regularization lambda value
regval <- 0.1 * max(abs(t(A) %*% b))

# Run standard ADMM
output_std <- admm_genlasso(A, b, D, lambda = regval)

# Run augmented ADMM
output_aug <- aug_admm_genlasso(A, b, D, M, lambda = regval)

# Function to calculate the difference at each iteration
calculate_mse_per_iteration <- function(x_iter, true_val) {
    apply(x_iter, 2, function(x) mean((x - true_val)^2))
}

# Calculate difference for standard and augmented ADMM at each iteration
mse_std_iter <- calculate_mse_per_iteration(output_std$x_iter, x0)
mse_aug_iter <- calculate_mse_per_iteration(output_aug$x_iter, x0)

# Print the length of iterations
cat("Number of iterations for Standard ADMM:", length(mse_std_iter), "\n")
cat("Number of iterations for Augmented ADMM:", length(mse_aug_iter), "\n")

# Create a data frame for plotting
iterations <- 1:max(length(mse_std_iter), length(mse_aug_iter))
mse_data

mse_data <- data.frame(
    iteration = iterations,
    MSE_std = c(mse_std_iter, rep(NA, length(iterations) - length(mse_std_iter))),
    MSE_aug = c(mse_aug_iter, rep(NA, length(iterations) - length(mse_aug_iter)))
)
# Reshape data for ggplot
mse_data_long <- melt(mse_data, id.vars = 'iteration', variable.name = 'method', value.name = 'Difference')

# Plot the results
ggplot(mse_data_long, aes(x = iteration, y = Difference, color = method)) +
    geom_line() +
    labs(title = "Mean Squared Error between Standard ADMM and Augmented ADMM",
         x = "Iteration", y = "Mean Squared Error") +
    theme_minimal()

# Save the plot to a file
ggsave("output/std_vs_aug_admm_precision_comparison.png", device = "png", width = 10, height = 6, dpi = 300)
