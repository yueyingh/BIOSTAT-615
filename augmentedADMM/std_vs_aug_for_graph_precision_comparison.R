###
# Script Purpose: Compare the precision of standard ADMM and augmented ADMM
# specifically for graph-structured data. The script produces a plot of the 
# Mean Squared Error (MSE) at each iteration for both methods.
###

setwd("augmentedADMM")

library(ggplot2)
library(Matrix)
library(augmentedADMM)
library(reshape2)
source("gen_data_v2.R") 

set.seed(123)

# Define parameters for generating the dataset
data_params <- list(
  num.groups = 100,              # Number of groups in the graph
  num.vars.per.group = 11,       # Number of variables per group
  n = 200,                       # Number of observations
  num.active.groups = 4,         # Number of active groups
  cor = 0.7,                     # Correlation coefficient
  err.var = 0.1                  # Error variance
)

data <- do.call(gen_data, data_params)

# Extracting matrices and vectors from the generated data
X <- data$X
Y <- data$Y
D <- data$A  # Incidence matrix for the graph
C <- data$C  # Constraint matrix
M <- data$M  # Matrix for augmented ADMM
true <- data$beta.true  # True regression coefficients

lambda1 <- 0.5
lambda2 <- 0.1

# Run standard ADMM
output_std <- admm_genlasso_for_graph(X, Y, D, C, lambda1, lambda2)

# Run augmented ADMM
output_aug <- aug_admm_genlasso_for_graph(X, Y, D, M, C, lambda1, lambda2)

# Function to calculate the difference at each iteration
calculate_mse_per_iteration <- function(x_iter, true_val) {
  apply(x_iter, 2, function(x) mean((x - true_val)^2))
}

# Calculate difference for standard and augmented ADMM at each iteration
mse_std_iter <- calculate_mse_per_iteration(output_std$x_iter, true)
mse_aug_iter <- calculate_mse_per_iteration(output_aug$x_iter, true)

# Print the length of iterations
cat("Number of iterations for Standard ADMM:", length(mse_std_iter), "\n")
cat("Number of iterations for Augmented ADMM:", length(mse_aug_iter), "\n")

# Create a data frame for plotting
iterations <- 1:max(length(mse_std_iter), length(mse_aug_iter))


mse_data <- data.frame(
  iteration = iterations,
  MSE_std = c(mse_std_iter, rep(NA, length(iterations) - length(mse_std_iter))),
  MSE_aug = c(mse_aug_iter, rep(NA, length(iterations) - length(mse_aug_iter)))
)

mse_data

# Reshape data for ggplot
mse_data_long <- melt(mse_data, id.vars = 'iteration', variable.name = 'method', value.name = 'Difference')

# Plot the results
ggplot(mse_data_long, aes(x = iteration, y = Difference, color = method)) +
  geom_line() +
  labs(title = "Mean Squared Error between Standard ADMM and Augmented ADMM for graph",
       x = "Iteration", y = "Mean Squared Error") +
  theme_minimal()

# Save the plot to a file
ggsave("output/std_vs_aug_admm_for_graph_precision_comparison.png", device = "png", width = 10, height = 6, dpi=300)

