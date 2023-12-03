###
# Testing with data generated from ADMM package
# Standard ADMM vs augADMM
###

library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)

set.seed(123)

# Define different problem sizes
problem_sizes <- seq(500, 1000, by = 50)
lambda_values <- 10^seq(-4, 0, length.out = 20)

# Data frame to store the results
benchmark_results <- data.frame(size = integer(), method = character(), time = numeric())

for (n in problem_sizes) {
    cat(sprintf("Starting loop iteration for problem size: %d\n", n))
    
    m <- 100
    p <- 0.1 # Percentage of non-zero elements

    cat("Generating x0...\n")
    x0 <- matrix(Matrix::rsparsematrix(n, 1, p))
    
    cat("Generating matrix A...\n")
    A <- matrix(rnorm(m * n), nrow = m)
    for (i in 1:ncol(A)) {
        A[, i] <- A[, i] / sqrt(sum(A[, i]^2))
    }

    cat("Generating vector b...\n")
    b <- A %*% x0 + sqrt(0.001) * matrix(rnorm(m))

    cat("Generating matrix D...\n")
    D <- diag(n)

    cat("Generating matrix M...\n")
    M <- norm(D, '2') * diag(n)

    cat("Setting regularization lambda value...\n")
    regval <- 0.1 * max(abs(t(A) %*% b))

    cat(sprintf("Benchmarking standard ADMM for size: %d...\n", n))
    admm_std_bench <- microbenchmark(
        admm_genlasso(A, b, D, lambda = regval),
        times = 10
    )
    cat("Standard ADMM benchmarking completed.\n")

    cat(sprintf("Benchmarking augmented ADMM for size: %d...\n", n))
    admm_aug_bench <- microbenchmark(
        aug_admm_genlasso(A, b, D, M, lambda = regval),
        times = 10
    )
    cat("Augmented ADMM benchmarking completed.\n")

    cat("Storing the results...\n")
    benchmark_results <- rbind(benchmark_results,
                               data.frame(size = n, method = "Standard ADMM", time = median(admm_std_bench$time)),
                               data.frame(size = n, method = "Augmented ADMM", time = median(admm_aug_bench$time)))
    cat(sprintf("Completed loop iteration for problem size: %d\n", n))
}


benchmark_results$time <- benchmark_results$time / 1e9

# Plot the results
ggplot(benchmark_results, aes(x = size, y = time, color = method)) +
    geom_line() +
    geom_point() +
    labs(title = "Run Time Comparison between Standard ADMM and Augmented ADMM",
         x = "Problem Size (n)", y = "Time (seconds)") +
    theme_minimal() +
    scale_y_continuous(labels = scales::comma) # To format the y axis with commas

# Save the plot to a file
ggsave("std_vs_aug_admm_runtime_comparison.pdf", width = 10, height = 6)



