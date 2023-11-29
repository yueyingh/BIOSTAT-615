library(augADMM)

set.seed(123) # For reproducibility
n <- 100 # Number of observations
p <- 1000 # Number of features

# Generate a sparse true coefficient vector
beta_true <- rep(0, p)
non_zero_indices <- sample(1:p, size = 100) # 10% sparsity
beta_true[non_zero_indices] <- rnorm(length(non_zero_indices))

# Generate the predictor matrix X
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# Generate the response vector y
noise <- rnorm(n, mean = 0, sd = 0.5)
y <- X %*% beta_true + noise

# Generate a graph structure (for simplicity, a random graph)
graph <- matrix(sample(0:1, p * p, replace = TRUE, prob = c(0.95, 0.05)), nrow = p)
graph <- ifelse(graph == t(graph), 1, 0) # Making it symmetric
diag(graph) <- 0 # No self-loops

# Convert the graph to an incidence matrix
edges <- which(graph == 1, arr.ind = TRUE)
num_edges <- nrow(edges)
A <- matrix(0, nrow = num_edges, ncol = p)

for (i in 1:num_edges) {
    A[i, edges[i, "row"]] <- -1
    A[i, edges[i, "col"]] <- 1
}

# Create the regularization matrix D (for simplicity, use the identity matrix)
D <- diag(p)

# Parameters for augADMM_genlasso
lambda <- 0.1 # Example value for lambda
rho <- 1.0 # Example value for rho
alpha <- 1.5 # Example value for alpha (should be in [1,2])
abstol <- 1e-4 # Example value for absolute tolerance
reltol <- 1e-2 # Example value for relative tolerance
maxiter <- 1000 # Maximum number of iterations

# Call the augADMM_genlasso function
result <- augADMM_genlasso(X, y, D, lambda, rho, alpha, abstol, reltol, maxiter)

# Check the results
print(result)
