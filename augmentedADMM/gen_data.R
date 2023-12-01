gen_data <- function(num.groups, num.vars.per.group, n, num.active.groups, cor, err.var) {
    # Initialize parameters
    p <- num.groups * num.vars.per.group # Total number of variables
    beta.true <- rep(0, p) # True beta coefficients
    X <- matrix(rnorm(n * p), n, p) # Generate random data matrix
    num.active.vars <- num.active.groups * num.vars.per.group # Number of active variables

    # Assign coefficients and induce correlation within groups
    for (group_idx in 1:num.groups) {
        # Assign non-zero coefficients for active groups
        if (group_idx <= num.active.groups) {
            start_idx <- 1 + (group_idx - 1) * num.vars.per.group
            end_idx <- group_idx * num.vars.per.group
            beta.true[start_idx:end_idx] <- floor((group_idx + 1) / 2) * (-1)^(group_idx + 1)
        }

        # Induce correlation within each group of variables
        for (var_idx in 2:num.vars.per.group) {
            base_var <- X[, 1 + (group_idx - 1) * num.vars.per.group]
            additional_var <- X[, var_idx + (group_idx - 1) * num.vars.per.group]
            X[, var_idx + (group_idx - 1) * num.vars.per.group] <- sqrt(cor) * base_var + sqrt(1 - cor) * additional_var
        }
    }

    # Normalize columns of X
    X <- sweep(X, 2, apply(X, 2, mean), "-")
    X <- sweep(X, 2, sqrt(n - 1) * apply(X, 2, sd), "/")

    # Generate response variable Y
    Y <- X %*% beta.true + sqrt(err.var) * rnorm(n)

    # Create erroneous edges for the graph
    num.err.edges <- (num.vars.per.group - 1) * num.active.vars
    m <- num.groups * num.vars.per.group * (num.vars.per.group - 1) / 2 + num.err.edges
    jdx <- rep(0, 2 * m)

    # Create edges between random nodes
    for (i in 1:num.active.vars) {
        a <- rep(i, num.vars.per.group - 1)
        b <- sample((1 + num.active.vars):p, num.vars.per.group - 1)
        jdx_indices <- 1 + 2 * (i - 1) * (num.vars.per.group - 1)
        jdx[jdx_indices:(jdx_indices + 2 * (num.vars.per.group - 1) - 1)] <- c(rbind(a, b))
    }

    # Generate edges within groups
    c <- NULL
    for (i in 1:(num.vars.per.group - 1)) {
        for (j in (i + 1):num.vars.per.group) {
            c <- c(c, c(i, j))
        }
    }
    for (group_idx in 1:num.groups) {
        jdx_indices <- 1 + 2 * num.err.edges + (group_idx - 1) * num.vars.per.group * (num.vars.per.group - 1)
        jdx[jdx_indices:(jdx_indices + 2 * num.vars.per.group * (num.vars.per.group - 1) - 1)] <- c + (group_idx - 1) * num.vars.per.group
    }

    # Create the incidence matrix C
    inc_mat <- matrix(0, nrow = m, ncol = p)
    for (i in 1:m) {
        inc_mat[i, jdx[2 * i - 1]] <- -1 # Outgoing edge
        inc_mat[i, jdx[2 * i]] <- 1 # Incoming edge
    }

    # Calculate node degrees and create M matrix
    node_degrees <- as.vector(table(jdx))
    M <- diag(2 * node_degrees + 1)

    # Combine identity and incidence matrices to create A
    A <- rbind(diag(p), inc_mat)

    return(list(Y = Y, X = X, A = A, M = M, val = val, idx = idx, jdx = jdx, p = p, m = m, beta.true = beta.true))
}
