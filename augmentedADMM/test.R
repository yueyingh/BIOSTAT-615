# Load the function
source("gen_data.R")

set.seed(123)
# Generate data
data_params <- list(
    num.groups = 5, num.vars.per.group = 10, n = 100,
    num.active.groups = 2, cor = 0.7, err.var = 0.5
)

data <- do.call(gen_data, data_params)

# Access the generated data
X <- data$X
Y <- data$Y
D <- diag(ncol(X))

lambda <- 0.1 # Example, adjust as needed
rho <- 1.0 # Example, adjust as needed
alpha <- 1.5 # Over-relaxation parameter, adjust as needed
abstol <- 1e-4 # Absolute tolerance for convergence
reltol <- 1e-2 # Relative tolerance for convergence
maxiter <- as.integer(1000) # Maximum number of iterations


xinit <- rep(0, ncol(X)) # Initial solution vector

# Load your package
library(augmentedADMM)

# result <- .Call('_augmentedADMM_admm_genlasso', X, Y, D, lambda, reltol, abstol, maxiter, rho)

# Call the function from your package
result <- genlasso_admm(X, Y, D, lambda, rho, alpha, abstol, reltol, maxiter)
result <- genlasso_admm_with_M()(X, Y, D, lambda, rho, alpha, abstol, reltol, maxiter)

# View the results
print(result)


