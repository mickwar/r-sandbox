### Density function and obtaining a random draw from the Wishart distribution

# get current list of objects in the console
TEMP_LIST_O = ls()

### dwishart()
# Calculate the density function of a Wishart(V, df) random variable
#
# Params: X   - p x p positive definite matrix (the support)
#         V   - p x p positive definite matrix (scale matrix)
#         df  - numeric that satisfies df > p - 1
#         log - logical, should the log density be evaluated?
dwishart = function(X, V, df, log = FALSE){
    p = nrow(V)
    if (df <= p - 1)
        stop("df must be greater than p - 1.")

    # trace function
    tr = function(x)
        sum(diag(x))

    if (log){
        # log multi gamma
        lmgamma = function(a, p)
            p*(p-1)/4*log(pi) + sum(lgamma(a+(1-1:p)/2))
        return ((df-p-1)/2*determinant(X)$modulus[1] - (1/2)*tr(solve(V) %*% X) - (df*p/2)*log(2) -
            (df/2)*determinant(V)$modulus[1] - lmgamma(df/2, p))
    } else {
        # multivariate gamma function
        mgamma = function(a, p)
            pi^(p*(p-1)/4) * prod(gamma(a+(1-1:p)/2))
        return (det(X)^((df-p-1)/2) * exp(-0.5*tr(solve(V) %*% X)) /
            (2^(df*p/2) * det(V)^(df/2) * mgamma(df/2, p)))
        }
    }


### random draws from a Wishart distribution (several implementations)
# Bartlett decomposition (for random draws)
# get one random draw
rwishart = function(V, df){
    p = nrow(V)
    # function to randomize lower triangle of a matrix with independent N(0,1)s
    lower.rand = function(x){
        n = nrow(x)
        z = sequence(n)
        index = cbind(
            row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
            col = rep(z[-length(z)], times = rev(tail(z, -1))-1))
        x[index] = rnorm(nrow(index))
        return (x)
        }
    L = t(chol(V))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A = lower.rand(A)
    X = L %*% A
    return (X %*% t(X))
#   return (L %*% A %*% t(A) %*% t(L))
    }

# Obtaining a random draw using the definition (these are not as fast as rwishart)
rwish2 = function(V, df){
    require(MASS)
    p = nrow(V)
    X = mvrnorm(df, rep(0, p), V)
    t(X) %*% X
    }

rwish3 = function(V, df, U = NULL){
    if (is.null(U))
        U = chol(V)
    p = NROW(U)
    X = matrix(rnorm(p*df, 0, 1), df, p) %*% U
    t(X) %*% X
    }


###
# compare TEMP_LIST with ls() which now has additional functions
TEMP_LIST_N = ls()
TEMP_LIST_N = TEMP_LIST_N[-which(TEMP_LIST_N == "TEMP_LIST_O")] # remove TEMP_LIST_O
TEMP_LIST_N = TEMP_LIST_N[!(TEMP_LIST_N %in% TEMP_LIST_O)] # get only added functions

# print the sourced functions
if (length(TEMP_LIST_N) > 0){
    cat("Sourced functions:\n    ")
    cat(TEMP_LIST_N, sep="\n    ")
    }

rm(TEMP_LIST_O, TEMP_LIST_N) # remove temporary variables
