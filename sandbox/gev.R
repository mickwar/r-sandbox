dgev = function(x, mu, sigma, ksi){
    supp = c(-Inf, Inf)
    if (ksi < 0)
        supp[2] = mu - sigma / ksi
    if (ksi > 0)
        supp[1] = mu - sigma / ksi

    tx = if (ksi == 0){
        exp(-(x - mu)/sigma)
    } else {
        (1 + ksi * (x - mu)/sigma)^(-1/ksi)
        }

    y = 1/sigma * tx^(ksi + 1) * exp(-tx)
    ifelse(x >= supp[1] & x <= supp[2], y, 0)
    }

pgev = function(x, mu, sigma, ksi){
    supp = c(-Inf, Inf)
    if (ksi < 0)
        supp[2] = mu - sigma / ksi
    if (ksi > 0)
        supp[1] = mu - sigma / ksi

    y = exp(-if (ksi == 0){
        exp(-(x - mu)/sigma)
    } else {
        (1 + ksi * (x - mu)/sigma)^(-1/ksi)
        })

    y[x <= supp[1]] = 0
    y[x >= supp[2]] = 1
    return (y)
    }

rgev = function(n, mu, sigma, ksi){
#   if (ksi == 0)
#       return (mu - sigma*log(-log(runif(n)))) # Gumbel
#   return (mu + sigma/ksi*((-log(runif(n)))^(-ksi) - 1))   # Frechet or Weibull
    if (length(ksi) == 1)
        ksi = rep(ksi, n)
    return (ifelse(ksi == 0, mu - sigma*log(-log(runif(n))),
        mu + sigma/ksi*((-log(runif(n)))^(-ksi) - 1)))
    }

#set.seed(5)
#n = 1000
#size = 100
#y = apply(matrix(rnorm(n*size), n, size), 1, max)
#
#n = 200
#y = rnorm(n)
#
#
#plot(density(y))
dat = read.table("~/files/data/coles_sea_level.txt", header = TRUE)
plot(dat, pch = 20)

y = dat$sea.level

calc.post = function(param){
    # likelihood
    out = sum(log(dgev(y, param[1], param[2], param[3])))
    # priors
    out = out + dnorm(param[1], prior.mu.a, prior.mu.b, log = TRUE)
    out = out + dgamma(param[2], prior.sigma.a, prior.sigma.b, log = TRUE)
    out = out + dnorm(param[3], prior.ksi.a, prior.ksi.b, log = TRUE)
    return (out)
    }

prior.mu.a = 0
prior.mu.b = 3
prior.sigma.a = 1/4
prior.sigma.b = 1/4
prior.ksi.a = 0
prior.ksi.b = 3


source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

nburn = 50000
nmcmc = 100000

nparam = 3
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

set.seed(1)
params[1,] = c(rnorm(1, prior.mu.a, prior.mu.b),
    rgamma(1, prior.sigma.a, prior.sigma.b),
    rnorm(1, prior.ksi.a, prior.ksi.b))

lower = c(-Inf, 0, -Inf)
upper = c(Inf, Inf, Inf)
window = 100

post = calc.post(params[1,])


for (i in 2:(nburn + nmcmc)){
    if (floor(i/window) == i/window)
        cat("\r", i, "/", nburn+nmcmc)
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i-1,], cand.sig)
    if (all(cand > lower) && all(cand < upper)){
        cand.post = calc.post(cand)
        if (log(runif(1)) <= cand.post - post){
            post = cand.post
            params[i,] = cand
            accept[i] = 1
            }
        }
    if ((floor(i/window) == i/window) && (i <= nburn))
        cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
            (cand.sig + window * var(params[(i-window+1):i,]) / i)
    }

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)

apply(params, 2, mean)
cov(params)

mean(accept)

hpds = apply(params, 2, hpd.uni)

par(mfrow = c(3,1))
plot(density(params[,1]), main = expression(mu), xlab = "", cex.main = 2)
abline(v=c(3.87, mean(params[,1])), col = c(2, 1))
abline(v=hpds[,1], lty = 2)
abline(v=c(3.82, 3.93), lty = 2, col = 2)
plot(density(params[,2]), main = expression(sigma), xlab = "", cex.main = 2)
abline(v=c(0.198, mean(params[,2])), col = c(2, 1))
abline(v=hpds[,2], lty = 2)
abline(v=c(.158,.238), lty = 2, col = 2)
plot(density(params[,3]), main = expression(xi), xlab = "", cex.main = 2)
abline(v=c(-0.05, mean(params[,3])), col = c(2, 1))
abline(v=hpds[,3], lty = 2)
abline(v=c(-.242, .142), lty = 2, col = 2)

par(mfrow = c(1,1))
plot(density(y))
xx = seq(min(density(y)$x), max(density(y)$x), length = 1000)
lines(xx, dgev(xx, mean(params[,1]), mean(params[,2]), mean(params[,3])), col = 'red')
lines(xx, dgev(xx, median(params[,1]), median(params[,2]), median(params[,3])), col = 'blue')

plot(density(y), ylim = c(0,2.4))
for (i in seq(1, nmcmc, by = 50))
    lines(xx, dgev(xx, params[i,1], params[i,2], params[i,3]), col = rgb(1,0,0,0.1))
lines(density(y), lwd=2)

