y = scan("~/files/R/651/data/faculty.dat")

n = length(y)

autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
window = 200

calc.post = function(params){
    require(truncnorm)

    # split
    thetas = params[1:n]
    sig2 = params[n+1]
    mu = params[n+2]
    tau2 = params[n+3]
    
    # likelihood
    out = sum(log(dtruncnorm(y, a = 1, b = 7, thetas, sqrt(sig2))))

    # priors
    out = out + sum(dnorm(thetas, mu, sqrt(tau2), log = TRUE))
    out = out + dgamma(sig2, a_sig, b_sig, log = TRUE)
    out = out + dnorm(mu, m, sqrt(s2), log = TRUE)
    out = out + dgamma(tau2, a_tau, b_tau, log = TRUE)

    return (out)
    }

m = 6
s2 = 0.5^2
a_sig = 1
b_sig = 1
a_tau = 1
b_tau = 1

nparams = n+3
nburn = 5000
nmcmc = 10000

upper = rep(Inf, nparams)
lower = rep(-Inf, nparams)
lower[(n+1):(n+3)] = 0

params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
sig = rep(1, nparams)

params[1, 1:n] = m
params[1, n+1] = a_sig * b_sig
params[1, n+2] = m
params[1, n+3] = a_tau * b_tau

curr.post = calc.post(params[1,])
cand.post = curr.post

for (i in 2:(nburn + nmcmc)){
    params[i,] = params[i-1,]
    cand.vec = params[i,]
    for (j in 1:nparams){
        cand.vec[j] = rnorm(1, params[i-1,j], sig[j])
        if (cand.vec[j] > lower[j] && cand.vec[j] < upper[j]){
            cand.post = calc.post(cand.vec)
            if (log(runif(1)) < cand.post - curr.post){
                curr.post = cand.post
                params[i,j] = cand.vec[j]
                accept[i,j] = 1
            } else {
                cand.vec[j] = params[i-1,j]
                }
        } else {
            cand.vec[j] = params[i-1,j]
            }
        }
    if (floor(i/window) == i/window && i <= nburn)
        sig = sig * autotune(apply(accept[(i-window+1):i,], 2, mean),
            k = max(1.1, window/50))
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

range(apply(accept, 2, mean))

for (i in 1:n)
    plot(params[,i], type='l')

plot(params[,n+1], type='l')
plot(params[,n+2], type='l')
plot(params[,n+3], type='l')

preds = double(nmcmc)
for (i in 1:nmcmc){
    temp.theta = rnorm(1, params[i,n+2], sqrt(params[i,n+3]))
    preds[i] = rtruncnorm(1, a=1, b=7, temp.theta,
        sqrt(params[i,n+1]))
    }

plot(density(y), col='blue')
points(density(preds), lwd=3, type='l')

rr = cor(params)
filled.contour(1:nparams, 1:nparams, rr)
