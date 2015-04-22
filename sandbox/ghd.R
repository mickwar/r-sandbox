### generalized hyperbolic distribution
source("~/files/R/mcmc/bayes_functions.R")

# modified bessel function of the second kind (called the third kind by R)
# besselK()

# density function
dgh = function(x, lambda, alpha, beta, delta, mu, log = TRUE){
    if (alpha > 0 && abs(alpha) > abs(beta) && delta > 0){
        gamma = sqrt(alpha^2 - beta^2)
        out = lambda * (log(gamma) - log(delta)) - 0.5*log(2*pi) - log(besselK(delta*gamma, lambda)) +
            beta*(x - mu) + log(besselK(alpha*sqrt(delta^2+(x-mu)^2), lambda-0.5)) - 
            (0.5-lambda)*(0.5*log(delta^2+(x-mu)^2) - log(alpha))
    } else {
        out = rep(-Inf, length(x))
        }
    if (log){
        return (out)
    } else {
        return (exp(out))
        }
    }

xx = seq(-15, 45, length = 1000)
lambda = 7
alpha = 1.2
beta = 0.8
delta = 0.5
mu = -9
plot(xx, dgh(xx, lambda, alpha, beta, delta, mu, log = FALSE), type='l')

### random draws using metropolis hastings
calc.post = function(params)
    dgh(params, lambda, alpha, beta, delta, mu)

nparams = 1

nburn = 50000
nmcmc = 2000000
upper = c(Inf)
lower = c(-Inf)
window = 100

sigs = old.sigs
#old.sigs = sigs

params = matrix(0, nburn + nmcmc, nparams)
params[1,] = 5
cand.param = params[1,]
accept = matrix(0, nburn + nmcmc, nparams)

post = calc.post(params[1,])
cand.post = post

set.seed(1)
for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            # check whether to accept draw or not
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
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(mean(accept[(i-window+1):i,]), k = max(window/50, 1.1))
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }


params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

mean(accept)

#plot(params, type='l')
plot(density(params))
lines(xx, dgh(xx, lambda, alpha, beta, delta, mu, log = FALSE), type='l', col='red')


y = sample(params, 1000)

### ### posterior parameter draws
calc.post = function(params)
    sum(dgh(y, params[1], params[2], params[3], params[4], params[5], log = TRUE))

nparams = 5
#sigs = rep(1, nparams)
sigs = old.sigs
#old.sigs = sigs

nburn = 15000
nmcmc = 100000
upper = c(Inf, Inf, Inf, Inf, Inf)
lower = c(-Inf, 0, 0, 0, -Inf)
window = 500

params = matrix(0, nburn + nmcmc, nparams)
#params[1,] = c(1, 1.5, 1, 1, 1)
params[1,] = old.params
#old.params= params[i,]
cand.param = params[1,]
accept = matrix(0, nburn + nmcmc, nparams)

post = calc.post(params[1,])
cand.post = post

for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            # check whether to accept draw or not
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
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

apply(accept, 2, mean)

# truth:
c(lambda, alpha, beta, delta, mu)

# estimated means
apply(params, 2, mean)

for (k in 1:nparams){
    plot(params[,k], type='l')
    readline()
    }

plot(params[,c(1,2)], type='l')
