##########
# random draws from a Wishart distribution
# trace function
tr = function(x)
sum(diag(x))

# multivariate gamma function
# p = 1 is univariate gamma
mgamma = function(a, p)
    pi^(p*(p-1)/4) * prod(gamma(a+(1-1:p)/2))
# log multi gamma
lmgamma = function(a, p)
    p*(p-1)/4*log(pi) + sum(lgamma(a+(1-1:p)/2))

# density of Wishart(V, n) at X
dwishart = function(X, V, df, log = FALSE){
    # X is p x p, positive definite (the support)
    # V is p x p, positive definite
    # df > p - 1 is degrees of freedom
    p = nrow(V)
    if (log){
        return ((df-p-1)/2*determinant(X)$modulus[1] - (1/2)*tr(solve(V) %*% X) - (df*p/2)*log(2) -
            (df/2)*determinant(V)$modulus[1] - lmgamma(df/2, p))
    } else {
        return (det(X)^((df-p-1)/2) * exp(-0.5*tr(solve(V) %*% X)) /
            (2^(df*p/2) * det(V)^(df/2) * mgamma(df/2, p)))
        }
    }

V = matrix(c(1, 0.9, 0.9, 1), 2, 2)
X = diag(2)
p = nrow(V)
df = 10

dwishart(X, V, df)
dwishart(X, V, df, TRUE)

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
    return (L %*% A %*% t(A) %*% t(L))
    }

rwish2 = function(V, df){
    require(MASS)
    p = nrow(V)
    X = mvrnorm(df, rep(0, p), V)
    t(X) %*% X
    }

rwishart(V, df)
rwish2(V, df)


p = 100
V = matrix(0.5, p, p)
diag(V) = 1
df = p


n = 1000
out = rep(list(matrix(0, p, p)), n)
out2 = rep(list(matrix(0, p, p)), n)
system.time(
    for (i in 1:n){
        out[[i]] = rwishart(V, df)
        }
    )
system.time(
    for (i in 1:n){
        out2[[i]] = rwish2(V, df)
        }
    )
mean.out = matrix(0, 2, 2)
mean.out2 = matrix(0, 2, 2)
for (i in 1:n){
    mean.out[1,1] = mean.out[1,1] + out[[i]][1,1]
    mean.out[1,2] = mean.out[1,2] + out[[i]][1,2]
    mean.out[2,1] = mean.out[2,1] + out[[i]][2,1]
    mean.out[2,2] = mean.out[2,2] + out[[i]][2,2]
    mean.out2[1,1] = mean.out2[1,1] + out2[[i]][1,1]
    mean.out2[1,2] = mean.out2[1,2] + out2[[i]][1,2]
    mean.out2[2,1] = mean.out2[2,1] + out2[[i]][2,1]
    mean.out2[2,2] = mean.out2[2,2] + out2[[i]][2,2]
    }

mean.out = mean.out / n
mean.out2 = mean.out2 / n

# theoretical mean: df*V
df * V
mean.out
mean.out2
##########
