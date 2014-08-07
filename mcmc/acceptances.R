library(foreach)
library(doMC)
registerDoMC(2)

dat = read.csv("../data/fgm.jazz.csv")
dat = dat[-c(67,68),]

# log-posterior stuff
like = function(lam)
    -sum(lam)+sum(y*log(lam))
prior = function(x)
    dnorm(x, 0, 10, log=TRUE)
lam = function(x, beta)
    exp(beta[1] + beta[2]*x[,1] + beta[3]*x[,2])

y = dat[,1]
x = cbind(dat[,2], 1*dat[,4])

# already burned in
mcmc = function(sigs = c(0.05, 0.0005, 0.07)){
    nmcmc = 10000
    params = matrix(0, nmcmc, 3)
    params[1,] = c(3.6132, -0.0000840, 0.040105)
    accept = matrix(0, nmcmc, 3)
    post = like(lam(x, params[1,]))+sum(prior(params[1,]))
    for (i in 2:nmcmc){
        params[i,] = params[i-1,]
        for (j in 1:3){
            cand = rnorm(1, params[i,j], sigs[j])
            cand.vec = params[i,]
            cand.vec[j] = cand
            cand.post = like(lam(x, cand.vec))+sum(prior(cand.vec))
            if (log(runif(1)) < cand.post - post){
                accept[i,j] = 1
                params[i,j] = cand
                post = cand.post
                }
            }
        }
    return (apply(accept, 2, mean))
    }

nreps = 100
out = matrix(0, nreps, 3)
sig1 = seq(0.001, 1, length=nreps)
sig2 = seq(0.00001, 0.01, length=nreps)
sig3 = seq(0.001, 1, length=nreps)
for (i in 1:nreps)
    out[i,] = mcmc(c(sig1[i], sig2[i], sig3[i]))

plot(sig1, out[,1], pch=20, ylim=c(0,1)); abline(h=0.25, col='red')
points(sig1, min(sig1)/(sig1), type='l')
plot(sig2, out[,2], pch=20, ylim=c(0,1)); abline(h=0.25, col='red')
plot(sig3, out[,3], pch=20, ylim=c(0,1)); abline(h=0.25, col='red')

plot(out[,1], type='l', ylim=c(0,1)); abline(h=0.25, col='red')
points(out[,2], type='l', col='green')
points(out[,3], type='l', col='blue')
xx=1:100
points(xx, 2/xx, type='l', col='red')

# get the distribution of acceptance rates
nreps = 100
overall.accept = foreach(1:nreps, .combine=rbind) %dopar% mcmc()

apply(overall.accept, 2, quantile, c(0.025, 0.975))
