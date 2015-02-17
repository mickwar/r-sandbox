### functions
source("~/files/R/mcmc/bayes_functions.R")

calc.post = function(param.vec, like.only = FALSE){
    # likelihood
    px = pnorm(x %*% param.vec)
    if (any(px == 0 | px == 1)) # when likelihood is 0
        return (-Inf)
    out = sum(y*log(px) + (1-y)*log(1-px))
    # priors (beta ~ N(0, 1))
    if (!like.only)
        out = out + sum(dnorm(param.vec, beta.mean, beta.var, log = TRUE))
    return (out)
    }

### read in data from arthurs github
library(RCurl)
dat = read.table(text = getURL("https://raw.githubusercontent.com/luiarthur/Fall2014/master/Stat637/5/sore.txt"),
    header = TRUE)
dat = dat[,-1]
dat = dat[order(dat[,3]),]
dat = dat[order(dat[,2]),]
dat = dat[order(dat[,1]),]

y = dat$Y
x = cbind(1, dat$D, dat$T)

### priors
beta.mean = c(0, 0, 0)
beta.cov = diag(3)
beta.var = diag(beta.cov)


### MCMC part
### initialize mcmc settings
nburn = 25000
nmcmc = 100000
window = 1000
nparams = 3
params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
lower = rep(-Inf, 3)
upper = rep(Inf, 3)

sigs = rep(1, nparams)
sigs = c(1.0878, 0.0293, 1.5996)

params[1,] = rep(0, 3)
post = calc.post(params[1,])
cand.param = params[1,]
post.mat = matrix(-Inf, nburn + nmcmc, nparams)
post.mat[1, 1] = post

for (i in 2:(nburn+nmcmc)){
#   cat("\rIteration",i,"/",nburn+nmcmc)
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
        post.mat[i, j] = post
        }
    # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
#   if (i == (nburn+nmcmc))
#       cat("\n")
    }

### remove burnin
params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

### posterior means
apply(params, 2, mean)

### acceptance rates
apply(accept, 2, mean)

### deviance (the full log likelihood is 0)
devs = -2*apply(params, 1, calc.post, like.only = TRUE)

### hpd
hpds = t(apply(params, 2, hpd.uni))
dev.hpd = hpd.uni(devs)

### trace and posterior density plots
pdf("./figs/mcmc_post.pdf", height = 12, width = 9)
main.vec = c("Intercept", "Duration", "Type")
par(mfrow=c(4,1), mar = c(2.1, 4.1, 4.1, 2.1))
plot.post(params[,1], density(params[,1]), hpds[1,], main = main.vec[1], cex = 2.5, xlab = "")
par(mfg=c(2,1,4,1))
plot.post(params[,2], density(params[,2]), hpds[2,], main = main.vec[2], cex = 2.5, xlab = "")
par(mfg=c(3,1,4,1))
plot.post(params[,3], density(params[,3]), hpds[3,], main = main.vec[3], cex = 2.5, xlab = "")
par(mfg=c(4,1,4,1))
plot.post(devs, density(devs), dev.hpd, main = "Deviance", cex = 2.5, xlab = "")
dev.off()

hpds
dev.hpd

# 44 minute surgery with T=1 tracheal tube
pred.x = c(1, 44, 1)
pred.p = pnorm(params %*% pred.x)

hpd.uni(pred.p)
mean(pred.p)
