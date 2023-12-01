library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)
source("gen_data_v2.R")

# Define the range of lambda values and other parameters for the experiment
lambda_values <- 10^seq(-4, 0, length.out = 20)
nu_values <- c(5, 10)

data_params <- list(
    num.groups = 5, num.vars.per.group = 10, n = 100,
    num.active.groups = 2, cor = 0.7, err.var = 0.5
)

runtime_results <- list()

# Generate the data
data <- do.call(gen_data, data_params)
X <- data$X
Y <- data$Y
D <- data$A
M <- data$M
C <- data$C

rho <- 1.0
abstol <- 1e-4
reltol <- 1e-2
maxiter <- 1000

# Loop over nu and lambda values
for (nu in nu_values) {
    for (lambda in lambda_values) {
        # Generate the D matrix based on the problem setup
        # D <- sparseMatrix(i = data$idx, j = data$jdx, x = data$val)


        # Run the standard ADMM with M matrix
        # Assuming M matrix has the same dimension as D
        # M <- sparseMatrix(i = data$idx, j = data$jdx, dims = c(nrow(D), ncol(D)))


        # Update the M matrix as required by your problem specifics
        # If the 'M' matrix needs to be different from 'D', adjust its construction accordingly

        genlasso_admm_aug_res <- microbenchmark(
            genlasso_admm_aug(X, Y, D, M, lambda, rho, alpha, abstol, reltol, maxiter),
            times = 10 # Increase if more precision is needed
        )

         # Run the augmented ADMM
        genlasso_admm_graph_res <- microbenchmark(
            genlasso_admm_for_graph(X, Y, D, M, C, lambda, lambda, rho, alpha, abstol, reltol, maxiter),
            times = 10 # Increase if more precision is needed
        )

        # Store results
        runtime_results[[paste("genlasso_admm_aug", "lambda", lambda, "nu", nu, sep = "_")]] <- genlasso_admm_aug_res
        runtime_results[[paste("genlasso_admm_for_graph", "lambda", lambda, "nu", nu, sep = "_")]] <- genlasso_admm_graph_res
    }
}

# Convert results to a data frame for plotting
results_df <- do.call(rbind, lapply(runtime_results, summary))
results_df$expr <- as.character(results_df$expr)

# Plotting the results
ggplot(results_df, aes(x = lambda, y = time, color = expr)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_continuous(trans = "log10") +
    labs(x = "Lambda", y = "Time (seconds)", color = "Method") +
    ggtitle("Runtime Comparison of ADMM Methods") +
    theme_minimal()
