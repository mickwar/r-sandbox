### Threshold exceedance example
### Daily rainfall at a location in south-west England from 1914-1962
### measured in mm

source("~/files/R/mcmc/bayes_functions.R")
library(MASS)

dgpd = function(x, mu, sigma, ksi){
    y = 1/sigma * (1 + ksi * (x - mu) / sigma)^(-1/ksi - 1)
    z = 1*(x >= mu)
    if (ksi < 0)
        z = z * (x <= mu - sigma/ksi)
    return (y*z)
    }
rgpd = function(n, mu, sigma, ksi){
    if (length(ksi) == 1)
        ksi = rep(ksi, n)
    return (ifelse(ksi == 0, mu - sigma * log(runif(n)),
        mu + sigma * (runif(n)^(-ksi) - 1) / ksi))
    }

dat = read.table("~/files/data/coles/daily_rain.txt", header = FALSE)[,1]
n = length(dat)
plot(dat, pch = 20)

threshold = 30 # Chosen in Coles

# Get the exceeded data
y = dat[dat > threshold] - threshold
k = length(y)
k/n

plot(density(y), pch = 20)

calc.post = function(params){
    sigma = params[1]
    ksi = params[2]
    zeta = params[3]
    # (lower) boundary check
    if (any(1 + ksi*y/sigma < 0))
        return (-Inf)
    # Likelihood
    if (ksi != 0){
        out = k*log(zeta) + (n-k)*log(1-zeta) - k*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*y/sigma))
#       out = -k*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*y/sigma))
    } else {
        out = k*log(zeta) + (n-k)*log(1-zeta) - k*log(sigma) - sum(y)/sigma
#       out = -k*log(sigma) - sum(y)/sigma
        }
    # Priors
    out = out + dgamma(sigma, a.sigma, b.sigma, log = TRUE)
    out = out + dnorm(ksi, m.ksi, s.ksi, log = TRUE)
    out = out + dbeta(zeta, a.zeta, b.zeta, log = TRUE)
    return (out)
    }

# Prior values
a.sigma = 2
b.sigma = 0.5
m.ksi = 0
s.ksi = 1
a.zeta = 1
b.zeta = 1

### MCMC
nburn = 25000
nmcmc = 100000

nparam = 3
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

set.seed(1)
params[1,] = c(
    rgamma(1, a.sigma, b.sigma),
    abs(rnorm(1, m.ksi, s.ksi)),
    rbeta(1, a.zeta, b.zeta))

#lower = c(0, -Inf)
#upper = c(Inf, Inf)
lower = c(0, -Inf, 0)
upper = c(Inf, Inf, 1)
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

apply(params, 2, mean)
cov(params)

# Standard errors
sqrt(diag(cov(params)))

# Get hpds (might take a while)
#hpds = apply(params, 2, hpd.uni)

# Equal-tailed 95%
hpds = apply(params, 2, quantile, c(0.025, 0.975))

par(mfrow = c(nparam,1))
plot(density(params[,1]), main = expression(sigma), xlab = "", cex.main = 2)
abline(v=mean(params[,1]), col = 1)
abline(v=hpds[,1], lty = 2)
#abline(v=7.44, col = 2)
#abline(v=c(3.82, 3.93), lty = 2, col = 2)

plot(density(params[,2]), main = expression(xi), xlab = "", cex.main = 2)
abline(v=mean(params[,2]), col = 1)
abline(v=hpds[,2], lty = 2)
#abline(v=0.184, col = 2)
#abline(v=c(-0.014,0.383), lty = 2, col = 2)

plot(density(params[,3]), main = expression(zeta), xlab = "", cex.main = 2)
abline(v=mean(params[,3]), col = 1)
abline(v=hpds[,3], lty = 2)
par(mfrow=c(1,1))

### Estimate zeta_u = Pr(X > u) with bootstrap
# zeta = double(nmcmc)
# for (i in 1:nmcmc)
#     zeta[i] = mean(sample(dat, replace = TRUE) > threshold)
# 
# mean(zeta*n)
# var(zeta*n)
# 

### 100-year return level
m = 365*100
xm = threshold + params[,1]/params[,2]*((m*params[,3])^params[,2]-1)
plot(density(xm))

mean(xm)
var(xm)
(int = quantile(xm, c(0.025, 0.975)))

mean(dat > int[1])
sum(dat > int[1])

### Diagnostics plots
# Computations

z = sort(y)
pplot = apply(params, 1, function(x) 1-(1+x[2]*z/x[1])^(-1/x[2]))
qplot = apply(params, 1, function(x) threshold+x[1]/x[2]*(((k:1)/(k+1))^(-x[2])-1))
lq1 = apply(pplot, 1, quantile, c(0.025, 0.975))
lm1 = apply(pplot, 1, mean)
lq2 = apply(qplot, 1, quantile, c(0.025, 0.975))
lm2 = apply(qplot, 1, mean)

#m = 365*seq(0.2, 1000, length = 100)
m = 365*exp(seq(-2.0, 7, length = 100))
#jj = seq(-1, 3.75+log10(365), by = 0.1)
#m = c(1/mean(params[,3]), 10^jj)
ZZ = apply(params, 1, function(x) threshold + x[1]/x[2] * ((m*x[3])^x[2] - 1))
Zm = apply(ZZ, 1, mean)
Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))


par(mfrow=c(2,2))
# Probability plot  (i/(n+1), squig(G)(z_(i)))
plot((1:k)/(k+1), lm1, pch = 20, main = "Probability plot",
    xlab = "Empirical", ylab = "Model")
segments((1:k)/(k+1), y0 = lq1[1,], y1 = lq1[2,], col = 'gray50')
abline(0, 1)

# Quantile plot     (z_(i), hat(G)^(-1)(z_(i)))
plot(threshold+z, lm2, pch = 20, main = "Quantile plot",
    xlab = "Empirical", ylab = "Model")
segments(threshold+z, y0 = lq2[1,], y1 = lq2[2,], col = 'gray50')
abline(0, 1)

# Return level plot
plot(log(m/365), Zm, ylim = range(Zq), type='l', axes = FALSE,
    xlim = c(log(0.1), log(1000)), xlab = "Return Period", ylab = "Return Level",
    main = "Return Level Plot")
sdat = sort(dat)
points(log((1/(1 - (1:n)/(n+1)) / 365)[sdat > threshold]),
    sdat[sdat > threshold], pch = 20)
axis(1, at = seq(log(0.1), log(1000), length = 5),
    labels = round(exp(seq(log(0.1), log(1000), length = 5)), 1))
axis(2)
box()
lines(log(m/365), Zq[1,], col = 'gray50')
lines(log(m/365), Zq[2,], col = 'gray50')

# Density plot
pred.y = rgpd(nmcmc, rep(threshold, nmcmc), params[,1], params[,2])
hist(threshold + y, col = 'gray', freq=  FALSE, breaks = 10, main = "Density Plot")
points(threshold + y, rep(0, length(y)), pch = 20)
lines(density(pred.y, n = 10000), col = 'red', lwd = 2)

