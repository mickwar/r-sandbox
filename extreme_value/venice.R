### Modeling the r Largest Order Statistics
### 3.5.3 of Coles (p. 69-72)
dat = read.table("~/files/data/pirazzoli_venice.txt", header = TRUE)
dat = tail(dat, 51)

year = dat$year

y = dat[,-1]
n = NROW(y)
ri = apply(y, 1, function(x) sum(!is.na(x)))

plot(0, type='n', xlim = range(year), ylim = range(y, na.rm = TRUE))
for (i in 1:n)
    points(rep(year[i], 10), y[i,], pch = 20)



calc.post = function(param){
    mu = param[1]
    sigma = param[2]
    ksi = param[3]
    # Check support
    if (!all(1+ksi*(y-mu)/sigma > 0, na.rm = TRUE))
        return (-Inf)

    # Likelihood
    if (ksi != 0){
        out = sum(apply(y, 1, function(x){
            x = x[!is.na(x)]
            r = length(x)
            -(1 + ksi * (x[r] - mu) / sigma) ^ (-1/ksi) - r*log(sigma) - (1/ksi + 1) * 
                sum(log(1 + ksi * (x - mu)/sigma))
            }))
    } else {
        out = sum(apply(y, 1, function(x){
            x = x[!is.na(x)]
            r = length(x)
            -exp(-(x[r] - mu) / sigma) - r*log(sigma) - sum((x - mu)/sigma)
            }))
        }

    # Priors
    out = out + dnorm(mu, prior.mu.a, prior.mu.b, log = TRUE)
    out = out + dgamma(sigma, prior.sigma.a, prior.sigma.b, log = TRUE)
    out = out + dnorm(ksi, prior.ksi.a, prior.ksi.b, log = TRUE)

    return (out)
    }

prior.mu.a = 100
prior.mu.b = 20
prior.sigma.a = 20
prior.sigma.b = 1
prior.ksi.a = 0
prior.ksi.b = 3


source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

nburn = 5000
nmcmc = 10000

nparam = 3
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

params[1,] = c(prior.mu.a, prior.sigma.a/prior.sigma.b, prior.ksi.a)

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

mean(accept)

### Moment estimates
apply(params, 2, mean)
var(params)

# Marginal SD
sqrt(diag(var(params)))

par(mfrow = c(3,1))
for (i in 1:3)
    plot(params[,i], type='l')

pairs(params, pch = 20)


muhat = mean(params[,1])
sighat = mean(params[,2])
ksihat = mean(params[,3])
z = sort(apply(y, 1, max, na.rm = TRUE))

par(mfrow=c(1,1))
# Probability plot
plot((1:n)/(n+1), exp(-(1+ksihat*(z - muhat)/sighat)^(-1/ksihat)), pch = 20)
abline(0, 1)

pplot = apply(params, 1, function(x) exp(-(1+x[3]*(z - x[1])/x[2])^(-1/x[3])))
qqlines = apply(pplot, 1, quantile, c(0.025, 0.975))
#matplot((1:n)/(n+1), pplot, type='l')

# Probability plot with bounds
plot((1:n)/(n+1), exp(-(1+ksihat*(z - muhat)/sighat)^(-1/ksihat)), pch = 20)
segments((1:n)/(n+1), y0 = qqlines[1,], y1 = qqlines[2,])
#lines((1:n)/(n+1), qqlines[1,], col = 'red')
#lines((1:n)/(n+1), qqlines[2,], col = 'red')
abline(0, 1)

# Quantile plot
plot(muhat - sighat/ksihat*(1-(-log((1:n)/(n+1)))^(-ksihat)), z, pch = 20)
# pplot = apply(params, 1, function(x) x[1] - x[2]/x[3]*(1-(-log((1:n)/(n+1)))^(-x[3])))
# qqlines = apply(pplot, 1, quantile, c(0.025, 0.975))
# lines(qqlines[1,], z, col = 'red')
# lines(qqlines[2,], z, col = 'red')
abline(0, 1)

# # Quantile plot (axes swapped)
# plot(z, muhat - sighat/ksihat*(1-(-log((1:n)/(n+1)))^(-ksihat)), pch = 20)
# pplot = apply(params, 1, function(x) x[1] - x[2]/x[3]*(1-(-log((1:n)/(n+1)))^(-x[3])))
# qqlines = apply(pplot, 1, quantile, c(0.025, 0.975))
# lines(z, qqlines[1,], col = 'red')
# lines(z, qqlines[2,], col = 'red')
# abline(0, 1)

muhat - sighat/ksihat*(1-(-log((1:n)/(n+1)))^(-ksihat))
returns = apply(params, 1, function(x)
    x[1] - x[2]/x[3]*(1-(-log((1:n)/(n+1)))^(-x[3])))

p = (1:n)/(n+1)
yp = -log(1-p)
muhat - sighat/ksihat*(1-yp^(-ksihat))
plot(log(yp), muhat - sighat/ksihat*(1-yp^(-ksihat)), pch = 20)

plot(returns[,9], z)
