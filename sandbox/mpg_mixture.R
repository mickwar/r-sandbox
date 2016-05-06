###########
# normal mixture model
library(date)
source("~/files/R/mcmc/bayes_functions.R")

gen = function(n){
    r = runif(n)
    l = sum(r < rho)
    out = double(n)
    out[r < rho] = rnorm(l, mu1, sig1)
    out[r >= rho] = rnorm(n-l, mu2, sig2)
    return (out)
    }


y = gen(500)

# mpg data
dat = read.table("~/files/data/MPG_saturn.txt", header = TRUE)
y = dat$miles / dat$gallons
#y = y[-c(13, 20)] # removing the most extreme values

### Simple linear regression
t = as.numeric(as.date(gsub("-","", dat[,1])))
x = diff(t)
y = y[-1]
t = t[-1]

par(mfrow = c(2,3), mar = c(4.1, 2.1,2.1, 1.1))
ord = order(x)
mod = lm(y[ord] ~ x[ord]) #Some ad hoc grouping, I figure <10 was road trip
summary(mod)
# predict on data
plot(x, y, pch=20, main = "MPG", xlab = "Days before last fill up", ylab = "MPG"); lines(x[ord], predict(mod), lwd = 3)

# residuals
plot(rstudent(mod), pch = 20, main = "Standardized Residuals"); abline(h = 0, lwd = 2)

# fitted vs observed
plot(y[ord], predict(mod), pch = 20, main = "Fitted vs Observed"); abline(0, 1)

### Natural splines (number of days between fill ups to predict MPG)
library(splines)
ord = order(x)
mod = lm(y[ord] ~ ns(x[ord], knots= c(10, 30, 60))) #Some ad hoc grouping, I figure <10 was road trip
deviance(mod)
summary(mod)

# predict on data
plot(x, y, pch=20, main = "MPG", xlab = "Days before last fill up", ylab = "MPG"); lines(x[ord], predict(mod), lwd = 3)

# residuals
plot(rstudent(mod), pch = 20, main = "Standardized Residuals"); abline(h = 0, lwd = 2)

# fitted vs observed
plot(y[ord], predict(mod), pch = 20, main = "Fitted vs Observed"); abline(0, 1)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

### MCMC
calc.post = function(params){
    rho1 = params[1]
    rho2 = params[2]
    mu1 = params[3]
    mu2 = params[4]
    mu3 = params[5]
    sig1 = params[6]
    sig2 = params[7]
    sig3 = params[8]
    out = 0
    # likelihood
    out = sum(log(
        rho1*dnorm(y, mu1, sig1) +
        rho2*dnorm(y, mu2, sig2) +
        (1-rho1-rho2)*dnorm(y, mu3, sig3)))
    # priors
    out = out + dbeta(rho1, 2, 3, log = TRUE)
    out = out + dbeta(rho2, 3, 2, log = TRUE)
    out = out + dnorm(mu1, 24, 1, log = TRUE)
    out = out + dnorm(mu2, 31, 1, log = TRUE)
    out = out + dnorm(mu3, 28, 5, log = TRUE)
    out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
    out = out + dgamma(sig2, 1.75, 0.5, log = TRUE)
    out = out + dgamma(sig3, 2.25, 0.5, log = TRUE)
    return (out)
    }
#rho1 = 0.7
#mu1 = 2
#mu2 = 6
#sig1 = 0.8
#sig2 = 1

y = dat$miles / dat$gallons
plot(density(y))
nparams = 8
sigs = rep(1, nparams)
nburn = 15000
nmcmc = 25000
upper = c(1, 1, Inf, Inf, Inf, Inf, Inf, Inf)
lower = c(0, 0, -Inf, -Inf, -Inf, 0, 0, 0)
window = 500

params = matrix(0, nburn + nmcmc, nparams)
params[1,] = c(0.3, 0.6, 0, 0, 0, 1, 1, 1)
cand.param = params[1,]
accept = matrix(0, nburn + nmcmc, nparams)

post = calc.post(params[1,])
cand.post = post

for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        if (j != 1 && j != 2){
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
        } else {
            if (cand >= 0 && cand <= (1-params[i, 3-j])){
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
        }
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }

params = tail(params, nmcmc)
accept = tail(accept, nmcmc)

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
mix = runif(nmcmc)
for (i in 1:nmcmc){
    if (mix[i] <= params[i,1])
        pred.y[i] = rnorm(1, params[i,3], sqrt(params[i,6]))
    if (mix[i] > params[i,1] && mix[i] <= sum(params[i,1:2]))
        pred.y[i] = rnorm(1, params[i,4], sqrt(params[i,7]))
    if (mix[i] > sum(params[i,1:2]))
        pred.y[i] = rnorm(1, params[i,5], sqrt(params[i,8]))
    }

plot(density(pred.y), col='green', lwd=3)
points(density(y), col='black', type='l', lwd=3)
legend("topleft", col=c("green", "black"), legend=c("Predictive", "Data"), lty=1, lwd=3)

##########
