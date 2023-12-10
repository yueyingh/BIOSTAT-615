###
# Script Purpose: Evaluate and compare the runtime of standard ADMM and augmented ADMM 
# specifically for graph-structured problems. The script generates a plot showing 
# runtime differences between the two methods across varying problem sizes.
###

setwd("augmentedADMM")

library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)
source("gen_data_v2.R")

set.seed(123)

# Define a range of problem sizes for benchmarking
problem_sizes <- seq(100, 200, by = 25)

# Data frame to store the results
benchmark_results <- data.frame(size = integer(), method = character(), time = numeric())

for (dim in problem_sizes) {
  data_params <- list(num.groups = dim, num.vars.per.group=11,
                      n = 200, num.active.groups=4, cor=0.7, err.var=0.1)
  data <- do.call(gen_data, data_params)

  # Extract matrices and vectors from generated data
  X <- data$X
  Y <- data$Y
  D <- data$A  # Incidence matrix for the graph
  C <- data$C  # Constraint matrix
  M <- data$M  # Matrix for augmented ADMM


  # Set regularization lambda value
  lambda1 <- 0.5
  lambda2 <- 0.1

  # Benchmark standard ADMM
  cat(sprintf("Benchmarking standard ADMM for graph for size: %d...\n", dim))
  admm_std_for_graph_bench <- microbenchmark(
    admm_genlasso_for_graph(X, Y, D, C, lambda1, lambda2),
    times = 10
  )

  # Benchmark augmented ADMM
  cat(sprintf("Benchmarking augmented ADMM for graph for size: %d...\n", dim))
  admm_aug_for_graph_bench <- microbenchmark(
    aug_admm_genlasso_for_graph(X, Y, D, M, C, lambda1, lambda2),
    times = 10
  )

  # Store the results
  benchmark_results <- rbind(benchmark_results,
                             data.frame(size = dim, method = "Standard ADMM for graph", time = median(admm_std_for_graph_bench$time)),
                             data.frame(size = dim, method = "Augmented ADMM for graph", time = median(admm_aug_for_graph_bench$time)))
}

benchmark_results$time <- benchmark_results$time / 1e9

# Plot the results
ggplot(benchmark_results, aes(x = size * 11, y = time, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Run Time Comparison between Standard ADMM and Augmented ADMM for graph",
       x = "Problem Size (p)", y = "Time (seconds)") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) # To format the y axis with commas

# Save the plot to a file
ggsave("output/std_vs_aug_admm_for_graph_runtime_comparison.png", device = "png", width = 10, height = 6, dpi = 300)
