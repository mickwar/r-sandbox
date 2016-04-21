### Multivariate analysis of the Iris data
library(MASS)
head(iris)

y = iris[,c(1,2)]
n = nrow(y)

rwishart = function(V, df){
    p = nrow(V)
    L = t(chol(V))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A[lower.tri(A)] = rnorm(p*(p-1)/2)
    X = L %*% A
    return (X %*% t(X))
    }

# Model:
# Likelihood:   y_i = (Sepal.Length_i, Sepal.Width_i), i=1:150
#               y_i | mu, V ~ MVN(mu, V), i=1:150
# Prior:        mu, V | m, b, S, r ~ NIW = N(mu | m, V/b) x IW(V | S, r)

m = apply(y, 2, mean)
b = 1
S = diag(2)
r = 10

S0 = matrix(0, 2, 2)
for (i in 1:n)
    S0 = S0 + matrix(as.numeric(y[i,] - m), 2, 1) %*% matrix(as.numeric(y[i,] - m), 1, 2)

nburn = 1000
nmcmc = 5000

params.mu = matrix(0, nburn+nmcmc, 2)
params.V = rep(list(matrix(0, 2, 2)), nburn + nmcmc)

for (i in 1:(nburn + nmcmc)){
    params.V[[i]] = solve(rwishart(solve(S + S0), r+n))
    params.mu[i,] = mvrnorm(1, m, (params.V[[i]])/(b+n))
    }
params.mu = tail(params.mu, nmcmc)
params.V = tail(params.V, nmcmc)

rbind(
    t(apply(params.mu, 2, quantile, c(0.025, 0.5, 0.975))),
    quantile(sapply(params.V, function(x) x[1,1]), c(0.025, 0.5, 0.975)),
    quantile(sapply(params.V, function(x) x[2,2]), c(0.025, 0.5, 0.975)),
    quantile(sapply(params.V, function(x) x[1,2]), c(0.025, 0.5, 0.975))
    )

