### simulate data

set.seed(1)
n = 40
g = 10
p = n/g

# Y = X*beta + Z*u + epsilon
X = cbind(rep(1, n), sort(rep(runif(g, -1, 1), each=p)))#,
#   sort(rep(runif(g), each=p)))
true.beta = c(0.4, -0.75)#, 0.0)

Z = kronecker(diag(g), rep(1, p))
g.sig = sqrt(0.25)
true.u = rnorm(g, 0, g.sig)
# u ~ N(0, 0.25*I)


# epsilon ~ N(0, 0.10*I)
e.sig = sqrt(0.10)
V = g.sig^2*Z %*% t(Z) + e.sig^2*diag(n)
# Y1 = rnorm(n, X %*% true.beta + Z %*% true.u, e.sig)

library(MASS)
Y = mvrnorm(1, X %*% true.beta, V)


# plot(X[,2], Y1, pch=20, col=rep(1:g, each=p))
# par(mfrow=c(1,2))
# Y = Y1
# plot(1:p, Y[1:p], type='l', col=1, ylim=range(Y))
# for (i in 2:g){
#     lines(1:p, Y[(p*i-(p-1)):(p*i)], col=i)
#     }
# Y = Y2
# plot(1:p, Y[1:p], type='l', col=1, ylim=range(Y))
# for (i in 2:g){
#     lines(1:p, Y[(p*i-(p-1)):(p*i)], col=i)
#     }

# bhat = solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% Y
# uhat = 0.25 * t(Z) %*% solve(V) %*% (Y - X %*% bhat)

autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

calc.post = function(params){
    beta = params[1:2]
    g.sig2 = params[3]
    r.sig2 = params[4]

    # V = ZGZ' + R
    V = g.sig2 * Z %*% t(Z) + r.sig2 * diag(n)
    A = matrix(g.sig2, p, p) + 

    # likelihood (multivariate normal)
    out = -0.5*determinant(V)$modulus[1] - 0.5*
        t(Y - X %*% beta) %*% solve(V) %*% (Y - X %*% beta)

    # priors
    out = out + sum(dnorm(beta, 0, 10, log = TRUE))
    out = out + dgamma(g.sig2, 0.25, 1, log = TRUE)
    out = out + dgamma(r.sig2, 0.15, 1, log = TRUE)
    
    return (out)
    }


nburn = 20000
nmcmc = 25000
nparams = 4
params = matrix(0, nburn+nmcmc, nparams)
# starting values
params[1,] = 1
sigs = rep(1, nparams)

lower = c(-Inf, -Inf, 0, 0)
upper = rep(Inf, 4)

# initialize log posterior value
post = calc.post(params[1,])
cand.param = params[1,]
keep.post = double(nburn+nmcmc)
keep.post[1] = post

accept = matrix(0, nburn+nmcmc, nparams)
window = 200

# mcmc loop
for (i in 2:(nburn+nmcmc)){
    cat("\rIteration",i,"/",nburn+nmcmc)
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            if (log(runif(1)) < cand.post - post){
                post = cand.post
                params[i,j] = cand
                accept[i,j] = 1
            } else {
                cand.param[j] = params[i,j]
                }
        } else {
            cand.param[j] = params[i,j]
            }
        }
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    keep.post[i] = post
    if (i == (nburn+nmcmc))
        cat("\n")
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]
keep.post = keep.post[(nburn+1):(nburn+nmcmc)]

plot(params[,1], type='l')
apply(accept, 2, mean)

plot(params[,c(1,3)], type='l')
plot(keep.post, type='l')

apply(params, 2, mean)
apply(params, 2, quantile, c(0.025, 0.50, 0.975))
