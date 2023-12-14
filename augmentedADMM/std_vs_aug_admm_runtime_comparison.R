###
# Script Purpose: Compare the runtime of standard ADMM and augmented ADMM
# on varying problem sizes. The script generates a plot showing the runtime
# differences between the two methods across different problem sizes.
###
 
setwd("augmentedADMM")

library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)

set.seed(123)

# Define a range of problem sizes for the comparison
problem_sizes <- seq(500, 1000, by = 50)
lambda_values <- 10^seq(-4, 0, length.out = 20)

# Data frame to store the results
benchmark_results <- data.frame(size = integer(), method = character(), time = numeric())

# Iterate over each problem size to perform benchmarking
for (n in problem_sizes) {
    m <- 100 # fixed number of columns
    p <- 0.1 # sparsity level

    # Generate a sparse matrix and corresponding response vector
    x0 <- matrix(Matrix::rsparsematrix(n, 1, p))
    A <- matrix(rnorm(m * n), nrow = m)
    for (i in 1:ncol(A)) {
        A[, i] <- A[, i] / sqrt(sum(A[, i]^2))
    }
    b <- A %*% x0 + sqrt(0.001) * matrix(rnorm(m))

    # Define matrices D and M for the augmented ADMM
    D <- diag(n)
    M <- norm(D, "2") * diag(n)

    # Set regularization lambda value
    regval <- 0.1 * max(abs(t(A) %*% b))

    # Benchmark standard ADMM
    cat(sprintf("Benchmarking standard ADMM for size: %d...\n", n))
    admm_std_bench <- microbenchmark(
        admm_genlasso(A, b, D, lambda = regval),
        times = 10
    )

    # Benchmark augmented ADMM
    cat(sprintf("Benchmarking augmented ADMM for size: %d...\n", n))
    admm_aug_bench <- microbenchmark(
        aug_admm_genlasso(A, b, D, M, lambda = regval),
        times = 10
    )

    # Store the results
    benchmark_results <- rbind(
        benchmark_results,
        data.frame(size = n, method = "Standard ADMM", time = median(admm_std_bench$time)),
        data.frame(size = n, method = "Augmented ADMM", time = median(admm_aug_bench$time))
    )
}

benchmark_results$time <- benchmark_results$time / 1e9

# Plot the results
ggplot(benchmark_results, aes(x = size * 11, y = time, color = method)) +
    geom_line() +
    geom_point() +
    labs(
        title = "Run Time Comparison between Standard ADMM and Augmented ADMM",
        x = "Problem Size (p)", y = "Time (seconds)"
    ) +
    theme_minimal() +
    scale_y_continuous(labels = scales::comma) # To format the y axis with commas

# Save the plot to a file
ggsave("output/std_vs_aug_admm_runtime_comparison.png", device = "png", width = 10, height = 6, dpi = 300)
