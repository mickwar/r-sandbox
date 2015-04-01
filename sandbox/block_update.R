###########
# normal mixture model
library(MASS)
source("~/files/R/mcmc/bayes_functions.R")
rho = 0.90
mu1 = 2
mu2 = 10
sig1 = 2.0
sig2 = 1


gen = function(n){
    r = runif(n)
    l = sum(r < rho)
    out = double(n)
    out[r < rho] = rnorm(l, mu1, sig1)
    out[r >= rho] = rnorm(n-l, mu2, sig2)
    return (out)
    }

set.seed(4)
y = gen(100)

calc.post = function(params){
    rho = params[1]
    mu1 = params[2]
    mu2 = params[3]
    sig1 = params[4]
    sig2 = params[5]
    out = 0
    # likelihood
    out = sum(log(rho*dnorm(y, mu1, sig1) +
        (1-rho)*dnorm(y, mu2, sig2)))
    # priors
    out = out + dbeta(rho, 7, 3, log = TRUE)
    out = out + dnorm(mu1, 0, 5, log = TRUE)
    out = out + dnorm(mu2, 0, 5, log = TRUE)
    out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
    out = out + dgamma(sig2, 1.5, 0.5, log = TRUE)
    return (out)
    }

nparams = 5
sigs = diag(rep(1, nparams))
nburn = 15000
nmcmc = 5000
upper = c(1, Inf, Inf, Inf, Inf)
lower = c(0, -Inf, -Inf, 0, 0)
window = 500

params = matrix(0, nburn + nmcmc, nparams)
params[1,] = c(0.5, 0, 0, 1, 1)
cand.param = params[1,]
accept = double(nburn + nmcmc)

post = calc.post(params[1,])
cand.post = post
keep.sig = matrix(0, nburn/window, 5)
keep.sig[1,] = diag(sigs)

for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i,], sigs)
    if (all(cand >= lower & cand <= upper)){
        cand.post = calc.post(cand)
        # check whether to accept draw or not
        if (log(runif(1)) < cand.post - post){
            post = cand.post
            params[i,] = cand
            accept[i] = 1
            }
        }
    # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn){
        sigs = autotune(mean(accept[(i-window+1):i]), k = max(window/50, 1.1)) *
            var(params[(i-window+1):i,]) / 4
        keep.sig[i/window,] = diag(sigs)
        }
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }

plot(c(0.01, keep.sig[,4]), type='l')

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc)]

for (i in 1:nparams){
    plot(params[,i], type='l')
    if (i != nparams)
        readline()
    }

apply(params, 2, mean)
mean(accept)
# it's okay to estimate rho as 1-rho, just make the
# changes with the other parameters and it all
# works out
pairs(params, pch = 20)
# plot(params[,c(1,4)], type='l')

m = nrow(params)
pred.y = double(m)
mix = rbinom(m, 1, params[,1])
for (i in 1:m){
    if (mix[i] == 1){
        pred.y[i] = rnorm(1, params[i,2], params[i,4])
    } else {
        pred.y[i] = rnorm(1, params[i,3], params[i,5])
        }
    }

plot(density(pred.y), col='green', lwd=3)
points(density(y), col='black', type='l', lwd=3)
legend("topleft", col=c("green", "black"), legend=c("Predictive", "Data"), lty=1, lwd=3)

##########

apply(params, 2, sd) / sigs
