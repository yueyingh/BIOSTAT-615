library(microbenchmark)
library(ggplot2)
library(Matrix)
library(augmentedADMM)
source("gen_data_v2.R")

set.seed(123)
lambda_values <- 10^seq(-4, 0, length.out = 20)
nu_values <- c(5, 10)
group_counts <- c(20, 30)
num_vars_per_group <- 11

runtime_results <- list()

for (num_groups in group_counts) {
    data_params <- list(
        num.groups = num_groups, num.vars.per.group = num_vars_per_group, n = 200,
        num.active.groups = 4, cor = 0.7, err.var = 0.1
    )

    # Generate the data
    cat("Generating data...\n")
    data <- do.call(gen_data_v2, data_params)
    X <- data$X
    Y <- data$Y
    D <- data$A

    rho <- 1.0
    abstol <- 1e-4
    reltol <- 1e-2
    maxiter <- 1000

    for (nu in nu_values) {
        for (lambda in lambda_values) {
            cat(sprintf("Running benchmarks for nu: %d, lambda: %f...\n", nu, lambda))
            # genlasso_admm_std_res <- microbenchmark(
            #     genlasso_admm(X, Y, D, lambda, rho, alpha, abstol, reltol, maxiter),
            #     times = 10 # Adjust as needed for precision
            # )

            genlasso_admm_aug_res <- microbenchmark(
                genlasso_admm_aug(X, Y, D, M, lambda, rho, alpha, abstol, reltol, maxiter),
                times = 10 # Increase if more precision is needed
            )

            # Store results for standard ADMM
            # runtime_results[[paste("standard_ADMM", "p", num_groups * num_vars_per_group, "nu", nu, "lambda", lambda, sep = "_")]] <- genlasso_admm_std_res

            # Store results for augmented ADMM for graph
            runtime_results[[paste("augmented_ADMM", "p", num_groups * num_vars_per_group, "nu", nu, "lambda", lambda, sep = "_")]] <- genlasso_admm_aug_res
        }
    }
    cat(sprintf("Finished processing group count: %d.\n", num_groups))
}

cat("Benchmarking process completed. Preparing data for plotting...\n")

# Prepare data for plotting
results_df <- data.frame(p = numeric(), nu = numeric(), lambda = numeric(), mean_time = numeric(), method = character())

for (key in names(runtime_results)) {
    cat("Key being processed:", key, "\n")
    res <- runtime_results[[key]]
    parts <- unlist(strsplit(key, "_"))

    # Print parts to see their structure
    cat("Parts extracted from key:", paste(parts, collapse = ", "), "\n")

    # Identify the correct indices for p, nu, and lambda based on method
    if (parts[1] == "standard" && parts[2] == "ADMM") {
        p_index <- 4
        nu_index <- 6
        lambda_index <- 8
    } else if (parts[1] == "augmented" && parts[2] == "ADMM" && parts[3] == "for") {
        p_index <- 6
        nu_index <- 8
        lambda_index <- 10
    } else {
        cat("Warning: Key structure is not as expected:", key, "\n")
        next
    }

    # Extract the numeric values for p, nu, and lambda
    p <- as.numeric(parts[p_index])
    nu <- as.numeric(parts[nu_index])
    lambda <- as.numeric(sub("e", "E", parts[lambda_index]))

    # Verify the conversion
    cat(sprintf("Converted p: %d, nu: %d, lambda: %f\n", p, nu, lambda))

    mean_time <- mean(res$time) / 1e6 # Convert to milliseconds

    # Determine the method based on key pattern
    method <- ifelse(grepl("for_graph", key), "ADMM for Graph", "Augmented ADMM")

    # Append to the data frame with the correct method label
    results_df <- rbind(results_df, data.frame(p, nu, lambda, mean_time, method))
}


cat("Data preparation complete. Starting to plot...\n")

# Plotting the results
pdf("admm_runtime_comparison_1.pdf", width = 10, height = 6)

# Define p values for plotting (calculated based on the number of groups)
p_values <- c(220, 330)

# Loop through each p and nu value combination and plot
for (p_value in p_values) {
    for (nu_value in nu_values) {
        cat(sprintf("Plotting for p = %d and nu = %d...\n", p_value, nu_value))
        # Filter the data for the current p and nu combination
        plot_data <- subset(results_df, p == p_value & nu == nu_value)

        # Create the plot
        p <- ggplot(plot_data, aes(x = lambda, y = mean_time, color = method)) +
            geom_line() +
            geom_point() +
            scale_x_log10() +
            scale_y_continuous(name = "Time (milliseconds)") +
            labs(
                title = paste("Runtime Comparison for p =", p_value, "and nu =", nu_value),
                x = "Lambda", color = "Method"
            ) +
            theme_minimal()

        # Print the plot to the PDF
        print(p)
    }
}

# Close the PDF device
dev.off()

cat("Plotting complete. PDF saved.\n")
