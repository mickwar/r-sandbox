### logit function
logit = function(x)
    log(x/(1-x))
### covariance(s) for simulator input (eq. 1)
make.cov = function(y, rho, lam){
    n = nrow(y)
    out = diag(n)
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            out[i,j] <- out[j,i] <- prod(rho^(4*(y[i,] - y[j,])^2))
            }
        }
    return (1/lam*out)
    }


make2 = function(y, rho, lam){
    n = nrow(y)
    out = matrix(1, n, n)
    for (i in 1:length(rho))
        out = out * rho[i]^(4 * as.matrix(dist(y[,i]))^2)
#   out = out + diag(n)
    return (1/lam*out)
    }

# with initialized distances (should work for x's and t's, but not thetas)
y.dist = rep(list(matrix(0, m, m)), p)
for (i in 1:p)
    y.dist[[i]] = 4*as.matrix(dist(y[,i]))^2

make3 = function(y, rho, lam){
    n = nrow(y)
    out = matrix(1, n, n)
    for (i in 1:length(rho))
        out = out * rho[i]^y.dist[[i]]
    return (1/lam*out)
    }


