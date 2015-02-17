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

y = dat$Y
x = cbind(1, dat$D, dat$T)

### priors
beta.mean = c(0, 0, 0)
beta.cov = diag(3)
beta.var = diag(beta.cov)
sig.beta = solve(t(x) %*% x + solve(beta.cov))


### Albert and Chib Gibbs Part
library(truncnorm)
library(MASS)

nburn = 1000
ngibb = 100000

params.beta = matrix(0, nburn + ngibb, 3)
params.z = matrix(0, nburn + ngibb, length(y))

# requires that the data be ordered by the response
z.lower = c(rep(-Inf, sum(y == 0)), rep(0, sum(y == 1)))
z.upper = c(rep(0, sum(y == 0)), rep(Inf, sum(y == 1)))

for (i in 2:(nburn + ngibb)){
#   cat("\rIteration",i,"/",nburn+ngibb)

    # draw beta
    params.beta[i,] = mvrnorm(1, sig.beta %*% t(x) %*% params.z[i-1,], sig.beta)

    # draw latent variable (z)
    params.z[i,] = rtruncnorm(1, a = z.lower, b = z.upper, mean = x %*% params.beta[i,])

#   if (i == (nburn+ngibb))
#       cat("\n")
    }

### remove burnin
params.beta = params.beta[(nburn + 1):(nburn + ngibb),]
params.z = params.z[(nburn + 1):(nburn + ngibb),]

### posterior means
apply(params.beta, 2, mean)

### deviance (the full log likelihood is 0) (is this correct?)
devs = -2*apply(params.beta, 1, calc.post, like.only = TRUE)

### hpd
hpds = t(apply(params.beta, 2, hpd.uni))
dev.hpd = hpd.uni(devs)

### trace and posterior density plots
pdf("./figs/gibb_post.pdf", height = 12, width = 9)
main.vec = c("Intercept", "Duration", "Type")
par(mfrow=c(4,1), mar = c(2.1, 4.1, 4.1, 2.1))
plot.post(params.beta[,1], density(params.beta[,1]), hpds[1,],
    main = main.vec[1], cex = 2.5, xlab = "")
par(mfg=c(2,1,4,1))
plot.post(params.beta[,2], density(params.beta[,2]), hpds[2,],
    main = main.vec[2], cex = 2.5, xlab = "")
par(mfg=c(3,1,4,1))
plot.post(params.beta[,3], density(params.beta[,3]), hpds[3,],
    main = main.vec[3], cex = 2.5, xlab = "")
par(mfg=c(4,1,4,1))
plot.post(devs, density(devs), dev.hpd, main = "Deviance", cex = 2.5, xlab = "")
dev.off()

hpds
dev.hpd

# 44 minute surgery with T=1 tracheal tube
pred.x = c(1, 44, 1)
pred.p = pnorm(params.beta %*% pred.x)

hpd.uni(pred.p)
mean(pred.p)
