library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)
source("gen_data_v2.R")

# Define the range of lambda values and other parameters for the experiment
lambda_values <- 10^seq(-4, 0, length.out = 20)
nu_values <- c(5, 10)

data_params <- list(
    num.groups = 5, num.vars.per.group = 10, n = 1000,
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

# head(names(runtime_results))
results_df <- data.frame(lambda = numeric(), nu = numeric(), mean_time = numeric(), method = character())

# Extract and store data from runtime_results
for (key in names(runtime_results)) {
    res <- runtime_results[[key]]
    parts <- unlist(strsplit(key, "_"))

    # Extract lambda and nu values, accounting for scientific notation in lambda
    lambda_index <- which(parts == "lambda") + 1
    nu_index <- which(parts == "nu") + 1
    lambda_str <- parts[lambda_index]
    nu <- as.numeric(parts[nu_index])

    # Handle scientific notation for lambda
    if (grepl("e", lambda_str)) {
        lambda <- as.numeric(sub("e", "E", lambda_str))
    } else {
        lambda <- as.numeric(lambda_str)
    }

    # Determine method based on key pattern
    method <- ifelse(grepl("for_graph", key), "ADMM for Graph", "Augmented ADMM")
    mean_time <- mean(res$time) / 1e6 # Convert to milliseconds

    # Append to the data frame
    results_df <- rbind(results_df, data.frame(lambda, nu, mean_time, method))
}

# Display the first few rows of the data frame to confirm correct extraction
head(results_df)

# Plotting the results
# Save the plot to a file
ggplot2::ggsave(
    plot = ggplot(results_df, aes(x = lambda, y = mean_time, color = method)) +
        geom_line() +
        geom_point() +
        scale_x_log10() +
        scale_y_continuous(name = "Time (milliseconds)") +
        labs(x = "Lambda", color = "Method") +
        ggtitle("Runtime Comparison of ADMM Methods") +
        theme_minimal(),
    filename = "admm_runtime_comparison.png",
    width = 10,
    height = 6,
    units = "in"
)
