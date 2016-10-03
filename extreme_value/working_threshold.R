### Treating the threshold as a parameter
### Behrens, et. al. (2003) Bayesian Analysis of Extreme Events with Threshold Estimation
library(MASS)
library(mwBASE)
library(truncnorm)
rgpd = function(n, mu, sigma, ksi){
    if (length(ksi) == 1)
        ksi = rep(ksi, n)
    return (ifelse(ksi == 0, mu - sigma * log(runif(n)),
        mu + sigma * (runif(n)^(-ksi) - 1) / ksi))
    }

n = 10000
set.seed(1)
x = rnorm(n, 0, 1)

#set.seed(1)
#u = 2.0
#x = c(rtruncnorm(round(n*0.90), a = -Inf, b = u, sd = 2), u + rgpd(n - round(n*0.90), 0, 3.2, -0.5))
plot(density(x))


nburn = 20000
nmcmc = 20000
window = 500

nparam = 5
params = matrix(0, nburn + nmcmc, nparam)
accept = double(nburn + nmcmc)
cand.sig = diag(0.1, nparam)

params[1,] = c(1, 0, quantile(x, 0.90), 0, 1)

lower = c(0, -Inf, -Inf, -Inf, 0)
upper = c(Inf, Inf, quantile(x, (n-20)/n), Inf, Inf)

calc.post = function(x, param){
    sigma = param[1]
    ksi = param[2]
    u = param[3]
    mu = param[4]
    tau = param[5]

    y = x[x > u] - u
    k = length(y)

    # (lower) boundary check
    if (any(1 + ksi*y/sigma < 0))
        return (-Inf)

    ### Likelihood
    # For exceedances
    if (ksi != 0){
        out = -k*log(sigma) - (1 + 1/ksi)*sum(log(1+ksi*y/sigma))
    } else {
        out = -k*log(sigma) - sum(y)/sigma
        }

#   # Original
#   out = out + k*log(pnorm(u, mu, tau, lower = FALSE))
#   # For center
#   out = out + sum(dnorm(x[x <= u], mu, tau, log = TRUE))

    # My version with continuity at u
    out = out + k*log(sigma) + k*dnorm(u, mu, tau, log = TRUE) 
    # For center
    out = out + sum(dnorm(x[x <= u], mu, tau, log = TRUE))
    out = out - n*log(pnorm(u, mu, tau) + sigma*dnorm(u, mu, tau))

    ### Priors
    out = out - log(sigma)
    out = out + dnorm(ksi, 0, 10, log = TRUE)
    out = out + dnorm(u, quantile(x, 0.90), 1, log = TRUE)
    out = out + dnorm(mu, 0, 10, log = TRUE)
    out = out - log(tau)
    return (out)
    }

post = calc.post(x, params[1,])

for (i in 2:(nburn + nmcmc)){
    if (floor(i/window) == i/window)
        cat("\r   ", i, "/", nburn+nmcmc)
    params[i,] = params[i-1,]
    cand = mvrnorm(1, params[i-1,], cand.sig)
    if (all(cand > lower) && all(cand < upper)){
        cand.post = calc.post(x, cand)
        if (log(runif(1)) <= cand.post - post){
            post = cand.post
            params[i,] = cand
            accept[i] = 1
            }
        }
    if ((floor(i/window) == i/window) && (i <= nburn))
        cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
            (cand.sig + window * var(params[(i-window+1):i,]) / i)
    if (i == (nburn + nmcmc))
        cat("\n")
    }

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)

par(mfrow = c(3,2))
for (i in 1:nparam)
    plot(params[,i], type='l')
par(mfrow = c(1,1))

# pairs(params, pch = 16, cex = 0.5)

mean(accept)
apply(params, 2, mean)

k = sapply(params[,3], function(z) sum(x > z))
hist(k, col = 'gray', border = 'white', freq = FALSE)

zeta = rbeta(nmcmc, k + 1/2, n - k + 1/2)
params = cbind(params, k)

pred.y = double(nmcmc)
v = (runif(nmcmc) < pnorm(params[,3], params[,4], params[,5]))
pred.y[v] = rtruncnorm(sum(v), a = -Inf, b = params[v,3], mean = params[v, 4], sd = params[v, 5])
pred.y[!v] = params[!v, 3] + rgpd(sum(!v), 0, params[!v, 1], params[!v, 2])

hist(x, col = 'gray', border = 'white', freq = FALSE, breaks = 30)
lines(density(pred.y, n = 10000), lwd = 3)


hpds = apply(params, 2, hpd_mult)
par(mfrow = c(3,2))
for (i in 1:ncol(params)){
    if (is.list(hpds)){
        plot_hpd(density(params[,i]), hpds[[i]], "dodgerblue")
    } else {
        plot_hpd(density(params[,i]), hpds[,i], "dodgerblue")
        }
    }
par(mfrow = c(1,1))

m = exp(seq(log(10), log(1000), length = 100))
rl = sapply(m, function(z) params[,3] + params[,1] / params[,2] * ( (z * zeta)^params[,2] - 1))

rq = apply(rl, 2, quantile, c(0.025, 0.975))
rm = apply(rl, 2, mean)

plot(log(m), rm, type='l', lwd = 3, ylim = range(rq), col ='blue')
polygon(c(log(m), rev(log(m))), c(rq[1,], rev(rq[2,])), col=col_fade("dodgerblue", 0.5), border = NA)
points(tail(log(1/((n:1)/(n+1))), mean(k)), tail(sort(x), mean(k)), pch = 16)
abline(h = params[,3], col = col_fade("gray99", 0.01))

