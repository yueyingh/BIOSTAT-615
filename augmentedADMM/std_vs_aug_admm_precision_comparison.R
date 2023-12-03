###
# Testing precision with data generated from ADMM package
# Standard ADMM vs augADMM
###

library(ggplot2)
library(Matrix)
library(augmentedADMM)
library(reshape2)

set.seed(123)

# Define a problem size for demonstration
n <- 200
m <- 100
p <- 0.1

x0 <- matrix(Matrix::rsparsematrix(n, 1, p))
A <- matrix(rnorm(m * n), nrow = m)
for (i in 1:ncol(A)) {
    A[, i] <- A[, i] / sqrt(sum(A[, i]^2))
}
b <- A %*% x0 + sqrt(0.001) * matrix(rnorm(m))
D <- diag(n)
M <- norm(D, "2") * diag(n)

# Set regularization lambda value
regval <- 0.1 * max(abs(t(A) %*% b))

# Run standard ADMM
output_std <- admm_genlasso(A, b, D, lambda = regval)
output_std
# Run augmented ADMM
output_aug <- aug_admm_genlasso(A, b, D, M, lambda = regval)

calculate_mse_per_iteration <- function(x_iter, true_val) {
    colSums((x_iter - true_val)^2) / n
}

# Calculate MSE for standard and augmented ADMM at each iteration
mse_std_iter <- calculate_mse_per_iteration(output_std$x_iter, x0)
mse_aug_iter <- calculate_mse_per_iteration(output_aug$x_iter, x0)

# Create a data frame for plotting
iterations <- 1:max(length(mse_std), length(mse_aug))

mse_data <- data.frame(
    iteration = iterations,
    MSE_std = c(mse_std, rep(NA, length(iterations) - length(mse_std))),
    MSE_aug = c(mse_aug, rep(NA, length(iterations) - length(mse_aug)))
)
mse_data
# Reshape data for ggplot
mse_data_long <- melt(mse_data, id.vars = 'iteration', variable.name = 'method', value.name = 'MSE')

# Plot the results
ggplot(mse_data_long, aes(x = iteration, y = MSE, color = method)) +
    geom_line() +
    labs(title = "MSE Comparison between Standard ADMM and Augmented ADMM",
         x = "Iteration", y = "Mean Squared Error (MSE)") +
    theme_minimal()

# Save the plot to a file
ggsave("std_vs_aug_admm_mse_comparison.pdf", width = 10, height = 6)

