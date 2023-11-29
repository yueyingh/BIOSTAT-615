# Load the function
source("gen_data.R")

set.seed(123) 
# Generate data
data_params <- list(num.groups = 5, num.vars.per.group = 10, n = 100,
                    num.active.groups = 2, cor = 0.7, err.var = 0.5)

data <- do.call(gen_data, data_params)
# Access the generated data
X <- data$X
Y <- data$Y
D <- diag(ncol(X))

lambda <- 0.1  # Example, adjust as needed
rho <- 1.0     # Example, adjust as needed
abstol <- 1e-4 # Absolute tolerance for convergence
reltol <- 1e-2 # Relative tolerance for convergence
maxiter <- 1000 # Maximum number of iterations

# Load your package
library(augADMM)

# Call the function from your package
result <- ADMM_genlasso(X, Y, D, lambda, rho, abstol, reltol, maxiter)

# View the results
print(result)
