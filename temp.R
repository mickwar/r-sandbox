nmcmc = 10000
params = double(nmcmc)
accept = double(nmcmc)

mu  = c(-2, 2)
sig = c(1, 1)

rprop = function(x){
    k = sample(length(mu)+1, 1)
    if (k == 1)
        return (rnorm(1, x, 1))
    return (rnorm(1, mu[k-1], sig[k-1]))
    }

dprop = function(x, y)
    sum(1/(length(mu)+1)*dnorm(x, mu, sig)) + 1/(length(mu)+1)*dnorm(x,y,1)
dprop.vec = function(x){
    out = double(length(x))
    for (i in 1:length(x))
        out[i] = sum(1/length(mu)*dnorm(x[i], mu, sig))
    return (out)
    }

for (i in 2:nmcmc){
    params[i] = params[i-1]
#   cand = rnorm(1, mu, sig)
#   A = (dnorm(cand, 0, 1, log = TRUE) - dnorm(params[i], 0, 1, log = TRUE)) +
#       (dnorm(params[i], mu, sig, log = TRUE) - dnorm(cand, mu, sig, log = TRUE))
    cand = rprop(params[i])
    A = (dnorm(cand, 0, 1, log = TRUE) - dnorm(params[i], 0, 1, log = TRUE)) +
        (log(dprop(params[i], cand)) - log(dprop(cand, params[i])))
    if (log(runif(1)) < A){
        params[i] = cand
        accept[i] = 1
        }
    }

par(mfrow = c(2, 1), mar = c(3.1, 2.1, 2.1, 1.1))
plot(tail(params, 10000), type='l')
plot(density(params))
curve(dnorm(x, 0, 1), col = 'red', add = TRUE)
curve(dprop.vec(x), col = 'blue', add = TRUE)
mean(accept)

quants = c(0.001, 0.005, 0.01, 0.02, 0.03, 0.45, 0.5, 0.55, 0.97, 0.98, 0.99, 0.995, 0.999)
cbind(quantile(params, quants), qnorm(quants))


dmix = function(x, y, p, n, sigma){
    out = double(length(x))
    for (i in 1:length(x))
        out[i] = sum(p * dnorm(x[i], y + 2*n*sigma, sigma))
    return (out)
    }
rmix = function(y, p, n, sigma){
    k = sample(n, 1, prob = p)
    rnorm(1, 

f(0, 2, p, n, sigma)
f(2, 0, p, n, sigma)

y = 0
n = -5:5
p = (1 / (1 + abs(n)^1)) / sum(1 / (1 + abs(n)^1))
p = 1/length(n)
sigma = 0.2

xx = seq(y - 20*sigma, y + 20*sigma, length = 500)
plot(xx, dmix(xx, y, p, n, sigma), type='l')
curve(dt(x, 1), add = TRUE, col = 'red')

A = 5
B = rgamma(1, A, 1)
dgamma(B, A, 1)
dgamma(A, B, 1)

curve(dgamma(x, A, 1), xlim = c(0, 15))
curve(dgamma(x, B, 1), xlim = c(0, 15), add = TRUE, col = 'red')
lines(rep(A, 2), c(0, dgamma(A, B, 1)))
lines(rep(B, 2), c(0, dgamma(B, A, 1)), col = 'red')
