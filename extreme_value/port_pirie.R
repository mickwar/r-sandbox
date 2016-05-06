### Annual Maximum Sea-levels at Port Pirie
### 3.4.1 of Coles (p. 59-64)
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




dat = read.table("~/files/data/coles/port_pirie.txt", header = TRUE)
plot(dat, pch = 20)

y = dat$sea.level
n = length(y)
plot(density(y))

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
window = 200

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

# Mean and covariances
apply(params, 2, mean)
cov(params)

# Standard errors
sqrt(diag(cov(params)))

# Get hpds (might take a while)
#hpds = apply(params, 2, hpd.uni)

# Equal-tailed 95%
hpds = apply(params, 2, quantile, c(0.025, 0.975))


# Marginal posteriors (black: my estimates / red: coles estimates)
#par(mfrow = c(3,1))
#par("mfg" = c(1,1,3,1))
#plot.post(params[,1], density(params[,1]), hpds[,1],
#    main = expression(mu), xlab = "", cex.main = 2)
#par("mfg" = c(2,1,3,1))
#plot.post(params[,2], density(params[,2]), hpds[,2],
#    main = expression(sigma), xlab = "", cex.main = 2)
#par("mfg" = c(3,1,3,1))
#plot.post(params[,3], density(params[,3]), hpds[,3],
#    main = expression(xi), xlab = "", cex.main = 2)
#par("new" = TRUE)


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


# Density estimates based on the mean (red) and median (blue) of the posteriors
#par(mfrow = c(1,1))
#plot(density(y))
#xx = seq(min(density(y)$x), max(density(y)$x), length = 1000)
#lines(xx, dgev(xx, mean(params[,1]), mean(params[,2]), mean(params[,3])), col = 'red')
#lines(xx, dgev(xx, median(params[,1]), median(params[,2]), median(params[,3])), col = 'blue')

#plot(density(y), ylim = c(0,2.4))
#for (i in seq(1, nmcmc, by = 50))
#    lines(xx, dgev(xx, params[i,1], params[i,2], params[i,3]), col = rgb(1,0,0,0.1))
#lines(density(y), lwd=2)

### Predictive distribution
pred.y = rgev(nmcmc, params[,1], params[,2], params[,3])
#plot(density(y), lwd = 2)
#lines(density(pred.y), col = 'red', lwd = 2)

### Diagnostic Plots
# Computations
pplot = apply(params, 1, function(x) exp(-(1+x[3]*(sort(y) - x[1])/x[2])^(-1/x[3])))
qplot = apply(params, 1, function(x) x[1] - x[2]/x[3]*(1-(-log((1:n)/(n+1)))^(-x[3])))
lq1 = apply(pplot, 1, quantile, c(0.025, 0.975))
lm1 = apply(pplot, 1, mean)
lq2 = apply(qplot, 1, quantile, c(0.025, 0.975))
lm2 = apply(qplot, 1, mean)

pp = rev(seq(0.001, 0.99, length = 100))
ZZ = apply(params, 1, function(x) x[1] - x[2]/x[3] * (1 - (-log(1-pp))^(-x[3])))
Zm = apply(ZZ, 1, mean)
Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))


par(mfrow=c(2,2))
# Probability plot  (i/(n+1), squig(G)(z_(i)))
plot((1:n)/(n+1), lm1, pch = 20, main = "Probability plot",
    xlab = "Empirical", ylab = "Model")
segments((1:n)/(n+1), y0 = lq1[1,], y1 = lq1[2,], col = 'gray50')
abline(0, 1)

# Quantile plot     (z_(i), hat(G)^(-1)(z_(i)))
plot(sort(y), lm2, pch = 20, main = "Quantile plot",
    xlab = "Empirical", ylab = "Model")
segments(sort(y), y0 = lq2[1,], y1 = lq2[2,], col = 'gray50')
abline(0, 1)

# Return level plot
plot(-log(-log(1-pp)), Zm, ylim = range(Zq), type='l', axes = FALSE,
    xlim = c(log(0.1), log(1000)), xlab = "Return Period", ylab = "Return Level",
    main = "Return Level Plot")
points(-log(-log((1:n)/(n+1))), sort(y), pch = 20)
axis(1, at = seq(log(0.1), log(1000), length = 5),
    labels = round(exp(seq(log(0.1), log(1000), length = 5)), 1))
axis(2)
box()
lines(-log(-log(1-pp)), Zq[1,], col = 'gray50')
lines(-log(-log(1-pp)), Zq[2,], col = 'gray50')

#plot(1/pp, Zm, ylim = range(Zq), type='l', axes = TRUE,
#    xlab = "Return Period", ylab = "Return Level",
#lines(-log(-log(1-pp)), Zq[1,], col = 'gray50')
#lines(-log(-log(1-pp)), Zq[2,], col = 'gray50')
#    main = "Return Level Plot")
#points(1/((n:1)/(n+1)), sort(y), pch = 20)
#lines(1/pp, Zq[1,], col = 'gray50')
#lines(1/pp, Zq[2,], col = 'gray50')

# Density plot
hist(y, col = 'gray', freq=  FALSE, breaks = 10, main = "Density Plot")
points(y, rep(0, length(y)), pch = 20)
lines(density(pred.y), col = 'red', lwd = 2)



###### FIX
# plot(0, type='n', axes = FALSE, xlab = "", ylab = "")
# # Return level plot
# 
# pp = seq(0.0005, 0.99, length = 1000)
# zp
# plot(-(log(1-pp)), qgev(1-pp, 0, 1, 0))
# 
# plot(-log(-log(1-pp)), qgev(pp, 0, 1, 0), type='l', ylim=c(-2, 15), axes = FALSE)
# axis(1, at = seq(-2, 8, by = 1), labels = round(exp(seq(-2, 8, by = 1)), 1))
# lines(-log(-log(1-pp)), qgev(pp, 0, 1, -0.2))
# lines(-log(-log(1-pp)), qgev(pp, 0, 1, 0.2))
# 
# 
# muhat = mean(params[,1])
# sighat = mean(params[,2])
# ksihat = mean(params[,3])
# yy = seq(min(y), max(y), length = 100)
# plot(yy, pgev(yy, muhat, sighat, ksihat))
# plot(return.period, return.period, log='y')
# 
# xx = seq(-4, 4, length = 1000)
# 
# 
# pp = seq(0.001, 0.99, length = 100)
# logyp = log(-log(1-pp))
# zp = muhat - sighat/ksihat * (1 - (-log(1-pp))^(-ksihat))
# 
# muhat - sighat/ksihat * (1 - (-log(1-0.01))^(-ksihat))
# 
# plot(-logyp, zp, pch = 20)
# 
# plot(exp(logyp), zp, pch = 20)
# plot(1/pp, zp, pch = 20)


# pp = seq(0.001, 0.99, length = 1000)
# plot(-log(-log(1-pp)), qgev(pp, 3.87, 0.198, -0.05), ylim=c(3.5, 5), type='l', axes = FALSE,
#     xlim = c(log(0.1), log(1000)))
# axis(1, at = seq(log(0.1), log(1000), length = 5),
#     labels = round(exp(seq(log(0.1), log(1000), length = 5)), 1))
# axis(2)
# box()
# lines(-log(-log(1-pp)), qgev(pp, 3.82, 0.158, -0.21), col = 'blue')
# lines(-log(-log(1-pp)), qgev(pp, 3.93, 0.238, 0.17), col = 'blue')
# 
# 
# logyp = seq(-2, 7, length = 100)
# pp = 
# 1-exp(-exp(logyp))
# zp = 0 - 1*log(-log(1-pp))
# 
# zp = qgev(1-exp(-logyp), 0, 1, 0)
# 
# plot(pp, zp)
# 
# qgev(pp, muhat, sighat, ksihat)
# plot(pp, qgev(1-pp, 0, 1, 0))
# 
# plot(log(1-pp), zp)
# plot(-log(1-pp), muhat-sighat/ksihat*(1-(-log(1-pp))^(-ksihat)), type='l')
# plot(muhat-sighat/ksihat*(1-(-log(1-pp))^(-ksihat)), log(1/pp), type='l')
######

