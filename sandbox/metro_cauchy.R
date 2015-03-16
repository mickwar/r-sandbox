### The metropolis algorithm doesn't perform too well with
### the cauchy distribution. Why would this be the case?
### Because of the undefined moments? Because of the heavy tails?

source("~/files/R/mcmc/bayes_functions.R")
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
window = 500

f.star = function(x)
    dcauchy(x)

xx = seq(-10, 10, length=100)
plot(xx, f.star(xx), type='l')

nmcmc = 100000
nburn = 500000
params = double(nmcmc+nburn)
accept = double(nmcmc+nburn)
sig = 1

curr.post = log(f.star(params[1]))
cand.post = curr.post

count = 1
keep.sig = NULL

for (i in 2:(nmcmc+nburn)){
    params[i] = params[i-1]
    cand = rnorm(1, params[i-1], sig)
    cand.post = log(f.star(cand))
    if (log(runif(1)) < cand.post - curr.post){
        curr.post = cand.post
        accept[i] = 1
        params[i] = cand
        }
    if (i <= nburn && floor(i/window) == i/window){
        sig = sig * autotune(mean(accept[(i-window+1):i]),
            target = 0.25, k = window/50)
        keep.sig[count] = sig
        count = count + 1
        }
    if (i == nburn)
        sig = calc.mode(density(keep.sig, n = 10000))
    }

params = params[(nburn+1):(nburn+nmcmc)]
accept = accept[(nburn+1):(nburn+nmcmc)]

mean(accept)
sig

plot(keep.sig, type = 'l')

plot(params, type='l', ylab="x", cex.lab = 1.5,
    main = "Trace Plot", cex.main = 1.5)

plot(density(params, n = 10000), xlim=c(-10, 10))
curve(dcauchy(x), from = -10, to = 10, add = TRUE,
    col = 'red')

qs = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5,
    0.9, 0.95, 0.99, 0.995, 0.999)
qparam = quantile(params, qs)
qtheory = qcauchy(qs)

cbind(qparam, qtheory)
