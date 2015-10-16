###########
# normal mixture model
library(date)
source("~/files/R/mcmc/bayes_functions.R")
rho = 0.7
mu1 = 2
mu2 = 6
sig1 = 0.8
sig2 = 1

gen = function(n){
    r = runif(n)
    l = sum(r < rho)
    out = double(n)
    out[r < rho] = rnorm(l, mu1, sig1)
    out[r >= rho] = rnorm(n-l, mu2, sig2)
    return (out)
    }

calc.post = function(params){
    rho = params[1]
    mu1 = params[2]
    mu2 = params[3]
    sig1 = params[4]
    sig2 = params[5]
    out = 0
    # likelihood
    out = sum(log(rho*dnorm(y, mu1, sig1) +
        (1-rho)*dnorm(y, mu2, sig2)))
    # priors
#   out = out + dbeta(rho, 7, 3, log = TRUE)
#   out = out + dnorm(mu1, 3, 10, log = TRUE)
#   out = out + dnorm(mu2, 3, 10, log = TRUE)
#   out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
#   out = out + dgamma(sig2, 1.5, 0.5, log = TRUE)
    out = out + dbeta(rho, 10, 15, log = TRUE)
    out = out + dnorm(mu1, 24, 1, log = TRUE)
    out = out + dnorm(mu2, 31, 1, log = TRUE)
    out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
    out = out + dgamma(sig2, 1.75, 0.5, log = TRUE)
    return (out)
    }

y = gen(500)

# mpg data
dat = read.table("~/files/data/MPG_saturn.txt", header = TRUE)
y = dat$miles / dat$gallons
#y = y[-c(13, 20)] # removing the most extreme values

### Natural splines (number of days between fill ups to predict MPG)
t = as.numeric(as.date(gsub("-","", dat[,1])))
x = diff(t)
y = y[-1]
t = t[-1]

plot(as.date(t), y, type='b', pch = 20)

ord = order(x)
mod = lm(y[ord] ~ ns(x[ord], knots= c(10, 30, 60))) #Some ad hoc grouping, I figure <10 was road trip
deviance(mod)
summary(mod)

plot(x, y, pch=20, main = "MPG on days since previous fill up"); lines(x[ord], predict(mod), lwd = 3)

plot(rstudent(mod), pch = 20, main = "Standardized Residuals"); abline(h = 0, lwd = 2)

plot(predict(mod), y[ord], pch = 20, main = "Fitted vs Observed"); abline(0, 1)

### MCMC
y = dat$miles / dat$gallons
plot(density(y))
nparams = 5
sigs = rep(1, nparams)
nburn = 15000
nmcmc = 25000
upper = c(1, Inf, Inf, Inf, Inf)
lower = c(0, -Inf, -Inf, 0, 0)
window = 500

params = matrix(0, nburn + nmcmc, nparams)
params[1,] = c(0.5, 0, 0, 1, 1)
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

for (i in 1:nparams){
    plot(params[,i], type='l')
    if (i != nparams)
        readline()
    }

apply(params, 2, mean)
apply(accept, 2, mean)
# it's okay to estimate rho as 1-rho, just make the
# changes with the other parameters and it all
# works out

m = nrow(params)
pred.y = double(m)
mix = rbinom(m, 1, params[,1])
for (i in 1:m){
    if (mix[i] == 1){
        pred.y[i] = rnorm(1, params[i,2], sqrt(params[i,4]))
    } else {
        pred.y[i] = rnorm(1, params[i,3], sqrt(params[i,5]))
        }
    }

plot(density(pred.y), col='green', lwd=3)
points(density(y), col='black', type='l', lwd=3)
legend("topleft", col=c("green", "black"), legend=c("Predictive", "Data"), lty=1, lwd=3)

##########
