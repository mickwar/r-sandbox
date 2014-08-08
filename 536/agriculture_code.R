# non-linearity
# gaussian process (via mcmc)

library(MASS)
source("~/files/R/pb_linux.R")

dat = read.table("./cwsi.csv", sep=",", header=T)
dat = dat[,-1]
y = as.matrix(dat$SWC)
x = as.matrix(dat$CWSI)

matern = function(x, y=NULL, nu, phi){
    n = 0
    m = 0    
    if (!is.null(x))
        n = nrow(x)
    if (!is.null(y))
        m = nrow(y)
    D = dist(rbind(x, y), upper=TRUE, diag=TRUE)
    as.matrix(1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*D)^nu *
        besselK(2*phi*sqrt(nu)*D, nu))+diag(n+m)
    }
new.seq = function(x, back){
    out = seq(x-back+1, x)
    return (out[out > 0])
    }
trailing = function(x, digits=3)
    formatC(x, digits=digits, format="f")
last.na = function(x) # not really, more like last non na
    x[max(which(!is.na(x)))]
cov.function = function(x, phi, y=NULL, i=NA){
    SIG2 = last.na(draws[["sig2"]])
    TAU2 = last.na(draws[["tau2"]])
    if (!is.na(i)){
        SIG2 = draws[["sig2"]][i]
        TAU2 = draws[["tau2"]][i]
        }
    if (is.null(y))
        return (SIG2 * matern(x, y, nu, phi) + diag(TAU2, nrow(x)))
    return (SIG2 * matern(x, y, nu, phi) + diag(TAU2, nrow(x)+nrow(y)))
    }

# choose phi
matern.cor.to.phi = function(d, nu, cor.target = 0.05, guess = 1,
    eps = 1e-6){
    phi = guess
    out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
        besselK(2*phi*sqrt(nu)*d, nu)
    adjust = 1
    while (out > cor.target + eps){
        phi = phi + adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out < cor.target){
            phi = phi - adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    adjust = 1
    if (phi <= 1)
        phi = phi + 1.1
    while (out < cor.target - eps){
        phi = phi - adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out > cor.target){
            phi = phi + adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    return (phi)
    }

nu = 2 # fix nu
compute.posterior = function(G){
    MU = last.na(draws[["mu"]])
    SIG2 = last.na(draws[["sig2"]])
    TAU2 = last.na(draws[["tau2"]])
    out = -1/2*determinant(G)$modulus[1] - 1/2*t(y-MU)%*%solve(G)%*%(y-MU) # log likelihood
    out = out + -1/(2*a.mu) * (MU - b.mu)^2 # log of normal
    out = out + (a.sig2-1)*log(SIG2) - b.sig2*SIG2 # log gamma (rate parametrization)
    out = out + (a.tau2-1)*log(TAU2) - b.tau2*TAU2
    return (as.vector(out))
    }

label = c("mu", "sig2", "tau2", "phi")
nburn = 2000
nmcmc = 10000
draws = NULL
accepts = NULL
sigmas = NULL
dampen = NULL
for (var.name in label[-which(label=="phi")]){
    draws[[var.name]] = NA*double(nburn+nmcmc)
    accepts[[var.name]] = double(nburn+nmcmc)
    sigmas[[var.name]] = 1
    dampen[[var.name]] = 1
    }
draws[["mu"]][1] = mean(y)
draws[["sig2"]][1] = var(y)
draws[["tau2"]][1] = 1
# support
upper.bound = list("mu"=Inf, "sig2"=Inf, "tau2"=Inf)
lower.bound = list("mu"=-Inf, "sig2"=0, "tau2"=0)
# prior values
a.mu = 100
b.mu = mean(y)
a.sig2 = 1.1
b.sig2 = 0.0001
a.tau2 = 1.1
b.tau2 = 0.0001
phi.vals = c(matern.cor.to.phi(d=0.25, nu=nu),
    matern.cor.to.phi(d=0.25, nu=nu),
    matern.cor.to.phi(d=0.5, nu=nu),
    matern.cor.to.phi(d=0.75, nu=nu),
    matern.cor.to.phi(d=1, nu=nu),
    matern.cor.to.phi(d=1.25, nu=nu),
    matern.cor.to.phi(d=1.50, nu=nu),
    matern.cor.to.phi(d=1.75, nu=nu),
    matern.cor.to.phi(d=2, nu=nu),
    matern.cor.to.phi(d=2.25, nu=nu),
    matern.cor.to.phi(d=2.50, nu=nu),
    matern.cor.to.phi(d=2.75, nu=nu),
    matern.cor.to.phi(d=3, nu=nu))
phi.vals = seq(0.01, 3, length=25)
draws[["phi"]] = phi.vals[1]
    


adjust = 0.5
back = 100
acc.bounds = c(0.2, 0.3) # the acceptance rates the sigmas will try
    # to reach for

cov.mat = cov.function(x, phi.vals[1])
post = compute.posterior(cov.mat)
for (iter in 2:(nburn+nmcmc)){
    pb.linux(iter-1, nburn+nmcmc-1)
    for (var.name in label){
        if (var.name == "phi"){
            llike = double(length(phi.vals))
            for(j in 1:length(phi.vals)){
                phi.cov = cov.function(x, phi.vals[j])
                mu.cov = last.na(draws[["mu"]])
                llike[j] = -1/2*determinant(phi.cov)$modulus[1] -
                    1/2*t(y-mu.cov)%*%solve(phi.cov)%*%(y-mu.cov)
                }
            # logarithm identity for log of sums
            log.sum = max(llike)+log(sum(exp(llike-max(llike)))) 
            probs = exp(llike-log.sum)
            draws[[var.name]][iter] = sample(phi.vals, 1, prob=probs)
        } else {
            hold = draws[[var.name]][iter-1]
            draws[[var.name]][iter] = hold
            cand = rnorm(1, hold, sigmas[[var.name]])
            if (cand > lower.bound[[var.name]] && cand < upper.bound[[var.name]]){
                draws[[var.name]][iter] = cand
                cand.covmat = cov.function(x, draws[["phi"]][iter-1])
                cand.post = compute.posterior(cand.covmat)
                if (log(runif(1)) < cand.post - post){
                    accepts[[var.name]][iter] = 1
                    cov.mat = cand.covmat
                    post = cand.post
                    if (mean(accepts[[var.name]][new.seq(iter, back)]) > acc.bounds[2])
                        sigmas[[var.name]] = sigmas[[var.name]]*(1+adjust/dampen[[var.name]])
                } else {
                    draws[[var.name]][iter] = hold
                    if (mean(accepts[[var.name]][new.seq(iter, back)]) < acc.bounds[1])
                        sigmas[[var.name]] = sigmas[[var.name]]*(1-adjust/dampen[[var.name]])
                    }
            } else {
                if (mean(accepts[[var.name]][new.seq(iter, back)]) < acc.bounds[1])
                    sigmas[[var.name]] = sigmas[[var.name]]*(1-adjust/dampen[[var.name]])
                }
            if (mean(accepts[[var.name]][new.seq(iter, back)]) >= acc.bounds[1] &&
                mean(accepts[[var.name]][new.seq(iter, back)]) <= acc.bounds[2])
                dampen[[var.name]] = dampen[[var.name]] + 1
#           if (post == cand.post){
#               cat("\n", var.name, rep(" ", 7-nchar(var.name)), trailing(mean(accepts[[var.name]][1:iter])),
#                   " ", trailing(post), " ", trailing(sigmas[[var.name]]), sep="")
#               }
            if (iter > nburn)
                dampen[[var.name]] = Inf
            }
        }

#   par(mfrow=c(length(label),1), mar=double(4)+0.2, oma=c(3,5,0,0))
#   for (var.name in label)
#       plot(draws[[var.name]][1:iter], type='l', xaxt='n')
#   mtext(iter, 1, 2)
#   mtext(label, 2, 3, at=seq(0.95, 0.05, length=length(label)),
#       outer=TRUE, cex=0.8)
    }

mean(draws[["mu"]][(nburn+1):(nburn+nmcmc)])
mean(draws[["sig2"]][(nburn+1):(nburn+nmcmc)])
mean(draws[["tau2"]][(nburn+1):(nburn+nmcmc)])
mean(draws[["phi"]][(nburn+1):(nburn+nmcmc)])


# predictions and intervals
n = nrow(x)
m = 100
predx = as.matrix(seq(0, 1, length=m))
post.pred = matrix(0, nmcmc, m)
for (i in 1:nmcmc){
    pb.linux(i, nmcmc)
    R.big = cov.function(predx, draws[["phi"]][i], x, i+nburn)
    R.11 = R.big[1:m, 1:m]
    R.12 = R.big[1:m, (m+1):(n+m)]
    R.22 = R.big[(m+1):(n+m), (m+1):(n+m)]
    mu = draws[["mu"]][i+nburn]
    inv = solve(R.22)
    post.pred[i,] = mvrnorm(1, mu+R.12%*%inv%*%(y-mu),
        R.11 - R.12%*%inv%*%t(R.12))
    }
upper = apply(post.pred, 2, quantile, 0.975)
mean.curve = apply(post.pred, 2, mean)
lower = apply(post.pred, 2, quantile, 0.025)

par(mfrow=c(1,1), mar=c(4,4,3,3))#, oma=c(5,3,3,1))
plot(dat, xlab="CWSI", ylab="SWC", ylim=c(min(lower),
    max(upper)), main="Fitted Model")
polygon(c(predx, predx[m:1]), c(upper, lower[m:1]), col='lightblue',
    border=NA)
points(predx, mean.curve, type='l', lwd=3)
points(dat, col='red', pch=1, lwd=2, cex=1.5)
points(predx, post.pred[1,], type='l', lwd=1, lty=2)
points(predx, post.pred[2,], type='l', lwd=1, lty=2)
