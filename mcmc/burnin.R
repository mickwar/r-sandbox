# playing around with different k values for
# autotune(), window size, and burn-in length
# to see if i can find good values

# thoughts: smaller windows should be associated with small k
# and larger windows with higher k.

# The larger window gives a more accurate estimate of the current
# acceptance rate, higher k causes greater adjustments to smaller
# changes in acceptance rates from window to window, so there is
# more fluctuation.

# smaller windows give less accurate acceptance rate estimates, so
# the sigmas should only be allowed to change by little

# it takes longer for high acceptance rates to be pulled down than
# it does for lower ones to catch up

# k = window/50 seems to work generally well


dat = read.csv("~/files/R/data/fgm.jazz.csv")
dat = dat[-c(67,68),]

datY = dat[,1]
datX = cbind(1, dat[,2], 1*dat[,4])
nparams = ncol(datX)

# log-posterior
calc.post = function(params){
    XB = datX %*% params
    -sum(exp(XB)) + sum(datY * XB) + 
        sum(dnorm(params, 0, 10, log=TRUE))
    }

# autotune candidate sigmas
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

nburn = 10000
nmcmc = 10000
window = 500
k = window/50

params = matrix(0, nburn+nmcmc, 3)
sigs = c(0.0000001, 1, 1)
accept = matrix(0, nburn+nmcmc, 3)
window.accept = matrix(0, nburn/window, 3)

post = calc.post(params[1,])


for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        cand.vec = params[i,]
        cand.vec[j] = cand
        cand.post = calc.post(cand.vec)
        if (log(runif(1)) < cand.post - post){
            accept[i,j] = 1
            params[i,j] = cand
            post = cand.post
            }
        }
    vv = i/window
    if (floor(vv) == vv && i <= nburn){
        window.accept[vv,] = apply(accept[(i-window+1):i,],2,mean)
        sigs = sigs * autotune(window.accept[vv,], k=k)
        }
    }

par(mfrow=c(3,1), mar=c(2.1,2.1,1.1,1.1))
for (i in 1:3){
    # x axis is the window chunk number
    # y axis is acceptance rate at that window
    plot(window.accept[,i], type='l', ylim=c(0, 1))
    abline(h=c(0.15, 0.40), col='blue', lty=2)
    abline(h=0.25, col='red', lty=2)
    }
apply(accept[(nburn+1):(nburn+nmcmc),], 2, mean)
sigs
