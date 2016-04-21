##########

# density of Wishart(V, df) at X
dwishart = function(X, V, df, log = FALSE){
    # X is p x p, positive definite (the support)
    # V is p x p, positive definite (scale matrix)
    # df > p - 1 is degrees of freedom

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

V = matrix(c(1, 0.9, 0.9, 1), 2, 2)
X = diag(2)
p = nrow(V)
df = 10

dwishart(X, V, df)
dwishart(X, V, 1, TRUE)

### random draws from a Wishart distribution (several implementations)
# Bartlett decomposition (for random draws)
# get one random draw
rwishart = function(V, df){
    p = nrow(V)
    L = t(chol(V))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[lower.tri(A)] = rnorm(p*(p-1)/2)
    X = L %*% A
    return (X %*% t(X))
#   return (L %*% A %*% t(A) %*% t(L))
    }

riwishart = function(V, df)
    solve(rwishart(V, df))

riwishart2 = function(V, df){
    p = nrow(V)
    L = chol(V)
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[upper.tri(A)] = rnorm(p*(p-1)/2)
    A = solve(A)
    L = solve(L)
    X = L %*% A
    return (X %*% t(X))
    }

riwishart3 = function(V, df){
    p = nrow(V)
    L = backsolve(chol(V), diag(p))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[upper.tri(A)] = rnorm(p*(p-1)/2)
    A = backsolve(A, diag(p))
    X = L %*% A
    return (X %*% t(X))
    }


# rwish2 = function(V, df){
#     require(MASS)
#     p = nrow(V)
#     X = mvrnorm(df, rep(0, p), V)
#     t(X) %*% X
#     }
# 
# rwish3 = function(V, df, U = NULL){
#     if (is.null(U))
#         U = chol(V)
#     p = NROW(U)
#     X = matrix(rnorm(p*df, 0, 1), df, p) %*% U
#     t(X) %*% X
#     }

rwishart(V, df)
rwish2(V, df)
rwish3(V, df)



p = 20
V = matrix(0.5, p, p)
diag(V) = 1
df = p

set.seed(1)
X1 = riwishart(V, df)

set.seed(1)
X2 = riwishart2(V, df)

set.seed(1)
X3 = riwishart3(V, df)





n = 10000
out = rep(list(matrix(0, p, p)), n)
out2 = rep(list(matrix(0, p, p)), n)
out3 = rep(list(matrix(0, p, p)), n)
system.time(
    for (i in 1:n){
        out[[i]] = riwishart(V, df)
        }
    )
system.time(
    for (i in 1:n){
        out2[[i]] = riwishart2(V, df)
        }
    )
system.time({
    U <- chol(V)
    for (i in 1:n){
        out3[[i]] = riwishart3(V, df)
        }
    })
mean.out = matrix(0, p, p)
mean.out2 = matrix(0, p, p)
mean.out3 = matrix(0, p, p)
for (i in 1:n){
    mean.out = mean.out + out[[i]]
    mean.out2 = mean.out2 + out2[[i]]
    mean.out3 = mean.out3 + out3[[i]]
    }

#VAR1 = matrix(0, p, p)
#VAR2 = matrix(0, p, p)
#for (i in 1:p){
#    for (j in 1:p){
#        VAR1[i,j] = var(sapply(out, function(x) x[i,j]))
#        VAR2[i,j] = var(sapply(out2, function(x) x[i,j]))
#        }
#    }


mean.out = mean.out / n
mean.out2 = mean.out2 / n
mean.out3 = mean.out3 / n

# theoretical mean: df*V
df * V
mean.out
mean.out2
mean.out3
##########
