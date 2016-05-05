### Multivariate extrems for Fremantle and Port Pirie annual maximum sea-levels
### 8.2.3 of Coles (p. 148-150)
source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

dgev = function(x, mu, sigma, ksi){
    n = length(x)
    lower.supp = rep(-Inf, n)
    upper.supp = rep(Inf, n)
    if (ksi < 0)
        upper.supp = mu - sigma / ksi
    if (ksi > 0)
        lower.supp = mu - sigma / ksi

    tx = if (ksi == 0){
        exp(-(x - mu)/sigma)
    } else {
        (1 + ksi * (x - mu)/sigma)^(-1/ksi)
        }

    y = 1/sigma * tx^(ksi + 1) * exp(-tx)
    ifelse(x >= lower.supp & x <= upper.supp, y, 0)
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

qgev = function(p, mu, sigma, ksi){
    if (ksi == 0){
        return (mu - sigma*log(-log(1-p)))
    } else {
        return (mu - sigma/ksi*(1-(-log(1-p))^(-ksi)))
        }
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


x = read.table("~/files/data/coles/fremantle.txt", header = TRUE)
y = read.table("~/files/data/coles/port_pirie.txt", header = TRUE)

# Get data with matching years
ind1 = which(x[,1] %in% y[,1])
ind2 = which(y[,1] %in% x[,1])
times = x[,1] - min(x[,1]) + 1
x = x[,2]
y = y[,2]
plot(x[ind1], y[ind2], pch = 20)


param = c(1.51, 0, 0.117, -0.149, 3.87, 0.197, -0.043, 0.922)

calc.post = function(param){
    beta0 = param[1]
    beta1 = param[2]
    mu1 = beta0 + beta1*times
    sig1 = param[3]
    ksi1 = param[4]
    mu2 = param[5]
    sig2 = param[6]
    ksi2 = param[7]
    alpha = param[8]
    if (ksi1 == 0){
#       z1 = exp((x - mu1) / sig1)
    } else {
        if (any(1 + ksi1*(x-mu1)/sig1 <= 0))
            return (-Inf)
#       z1 = (1+ksi1*(x - mu1)/sig1)^(1/ksi1)
        }
    if (ksi2 == 0){
#       z2 = exp((y - mu2)/ sig2)
    } else {
        if (any(1 + ksi2*(y-mu2)/sig2 <= 0))
            return (-Inf)
#       z2 = (1+ksi2*(y - mu2)/sig2)^(1/ksi2)
        }
    z1 = x
    z2 = y

    V = function(x, y, alpha)
        (x^(-1/alpha)+y^(-1/alpha))^alpha
    Vx = function(x, y, alpha)
        -x^(-1-1/alpha)*(x^(-1/alpha) + y^(-1/alpha))^(alpha-1)
    Vy = function(x, y, alpha)
        -y^(-1-1/alpha)*(x^(-1/alpha) + y^(-1/alpha))^(alpha-1)
    Vxy = function(x, y, alpha)
        (alpha - 1) / alpha*(x*y)^(-1-1/alpha)*(x^(-1/alpha)+y^(-1/alpha))^(alpha-2)

    # likelihood (joint part)
    out = (Vx(z1[ind1], z2[ind2], alpha) * Vy(z1[ind1], z2[ind2], alpha) -
        Vxy(z1[ind1], z2[ind2], alpha)) * exp(-V(z1[ind1], z2[ind2], alpha))
#   log(Vx(z1[ind1], z2[ind2], alpha) * Vy(z1[ind1], z2[ind2], alpha) -
#       Vxy(z1[ind1], z2[ind2], alpha)) 
#   sum(-V(z1[ind1], z2[ind2], alpha))
    if (any(is.na(out)))
        out = -Inf
    if (any(out == 0))
        out = -Inf
    if (out[1] != -Inf)
        out = sum(log(out))

    # likelihood (marginal parts)
    out = out + sum(log(dgev(x, mu1, sig1, ksi1)))
    out = out + sum(log(dgev(y, mu2, sig2, ksi2)))

    # priors
    out = out +  dnorm(beta0, beta0.a, beta0.b, log = TRUE)
    out = out +  dnorm(beta1, beta1.a, beta1.b, log = TRUE)
    out = out + dgamma(sig1,  sig1.a,  sig1.b, log = TRUE)
    out = out +  dnorm(ksi1,  ksi1.a,  ksi1.b, log = TRUE)
    out = out +  dnorm(mu2,   mu2.a,   mu2.b, log = TRUE)
    out = out + dgamma(sig2,  sig2.a,  sig2.b, log = TRUE)
    out = out +  dnorm(ksi2,  ksi2.a,  ksi2.b, log = TRUE)
    out = out +  dbeta(alpha, alpha.a, alpha.b, log = TRUE)
    return (out)
    }

# Hyperpriors
beta0.a = 0; beta0.b = 3; beta1.a = 0; beta1.b = 3; sig1.a = 1/4; sig1.b = 1/4; ksi1.a = 0; ksi1.b = 3
mu2.a = 0; mu2.b = 3; sig2.a = 1/4; sig2.b = 1/4; ksi2.a = 0; ksi2.b = 3
alpha.a = 1; alpha.b = 1

nburn = 20000
nmcmc = 10000

nparam = 8
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

set.seed(1)
params[1,] = c(0, 0, 1, 0, 0, 1, 0, 0.5)

lower = c(-Inf, -Inf, 0, -Inf, -Inf, 0, -Inf, 0)
upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, 1)
window = 200

ind1 = 1:60
ind2 = 1:60

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
mean(accept)

par(mfrow = c(3,3), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:nparam){
    plot(params[,i], type='l')
    abline(h = mean(params[,i]), col = 'red', lwd = 2)
    abline(h = quantile(params[,i], c(0.025, 0.975)), lty=2, col = 'red')
    }
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

par(mfrow = c(3,3), mar = c(4.1, 2.1, 2.1, 1.1))
for (i in 1:nparam)
    acf(params[,i])
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))

apply(params, 2, mean)
sqrt(diag(cov(params)))

apply(params, 2, quantile, c(0.025, 0.5, 0.975))


#####
V = function(x, y, alpha)
    (x^(-1/alpha)+y^(-1/alpha))^alpha
Vx = function(x, y, alpha)
    -x^(-1-1/alpha)*(x^(-1/alpha) + y^(-1/alpha))^(alpha-1)
Vy = function(x, y, alpha)
    -y^(-1-1/alpha)*(x^(-1/alpha) + y^(-1/alpha))^(alpha-1)
Vxy = function(x, y, alpha)
    (alpha - 1) / alpha*(x*y)^(-1-1/alpha)*(x^(-1/alpha)+y^(-1/alpha))^(alpha-2)


like = function(param){
    x = param[1]
    y = param[2]
    log((Vx(x, y, alpha) * Vy(x, y, alpha) - Vxy(x, y, alpha)) * exp(-V(x, y, alpha)))
    }


nburn = 50000
nmcmc = 50000

nparam = 2
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

set.seed(1)
params[1,] = c(1, 1)

lower = c(0, 0)
upper = c(Inf, Inf)
window = 200

alpha = 0.1
post = like(params[1,])

for (i in 2:(nburn + nmcmc)){
    if (floor(i/window) == i/window)
        cat("\r", i, "/", nburn+nmcmc)
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i-1,], cand.sig)
    if (all(cand > lower) && all(cand < upper)){
        cand.post = like(cand)
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

mean(accept)

plot(params, pch = 20, type='l')
segments(params[-nmcmc,1], params[-nmcmc,2], params[-1,1], params[-1,2],
    col = rgb(seq(0, 1, length = nmcmc - 1), 0, 0))

plot(params[,1], type='l')
plot(params[,2], type='l')

params = params[sample(nmcmc, 60),]
x = params[,1]
y = params[,2]

