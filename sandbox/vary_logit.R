##########
# varying base logistic
logit = function(x, b = exp(1))
    log(x/(1-x), base = b)
logistic = function(x, b = exp(1))
    1 / (1 + b^(-x))
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
calc.mode = function(dens, method="loess"){
    dx = dens$x
    dy = dens$y
    if (method == "loess"){
        l = loess(dy ~ dx)
        return (l$x[which.max(l$y)])
        }
    if (method == "density")
        return (dx[which.max(dy)])
    }

calc.post = function(params){
    beta = params[1:2]
    b = 1/params[3]
    # likelihood
    lx = logistic(x %*% beta, b)
    keep = which(lx != 1 & lx != 0)
    out = sum(y[keep]*log(lx[keep]) + (1-y[keep])*log(1-lx[keep]))
    # priors (beta ~ N(0, 1), b ~ Unif(1.5, 10))
    out = out + sum(dnorm(beta, 0, 1, log = TRUE))
    out = out + dunif(1/b, b.lower, b.upper, log = TRUE)
    return (out)
    }

# (y*log(logistic(x %*% beta, b)) + (1-y)*log(1-logistic(x %*% beta, b)))
# which(logistic(x %*% beta, b) == 1)

set.seed(1)
n = 1000
true.beta = c(-1.00, 1.50)
true.b = 2
b.lower = 0
b.upper = 1
x = cbind(1, sort(rnorm(n, 0, 1)))
y = rbinom(n, 1, prob = logistic(x %*% true.beta, true.b))

par(mfrow=c(1,1))
plot(x[,2], y, pch = 20)
points(x[,2], logistic(x[,2], exp(1)), type='l')
points(x[,2], logistic(x[,2], 2), type='l')
points(x[,2], logistic(x[,2], 5), type='l')

# initialize mcmc settings
nburn = 5000
nmcmc = 10000
window = 500
nparams = 3
params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
lower = c(-Inf, -Inf, b.lower)
upper = c(Inf, Inf, b.upper)

params[1,] = c(true.beta, (b.lower + b.upper)/2)
sigs = c(0.1, 0.1, 0.01)

post = calc.post(params[1,])
cand.param = params[1,]
post.mat = matrix(-Inf, nburn + nmcmc, nparams)
post.mat[1, 1] = post

for (i in 2:(nburn+nmcmc)){
    cat("\rIteration",i,"/",nburn+nmcmc)
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
    if (i == (nburn+nmcmc))
        cat("\n")
    }

# plot(x[,2], y, pch = 20)
# points(x[,2], logistic(x %*% cand.param[1:2], cand.param[3]), type='l')
# points(x[,2], logistic(x %*% params[i,1:2], params[i,3]), type='l')

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

par(mfrow=c(3,1))
plot(params[,1], type='l')
plot(params[,2], type='l')
plot(params[,3], type='l')
 
calc.mode(density(1/params[,3]))

apply(params, 2, mean)
apply(accept, 2, mean)

max(post.mat)
# with b in (1.5, 10)
# n = 200, true.beta = c(-1.0, 1.5), true.b = 5 (fixed at 5): -82.36014
# n = 200, true.beta = c(-1.0, 1.5), true.b = 5 (fixed at exp(1)): -84.26674
# n = 200, true.beta = c(-1.0, 1.5), true.b = 5 (fixed at 2): -87.12315
# n = 200, true.beta = c(-1.0, 1.5), true.b = 5 (varying between 1.5 and 10): -81.70131
# the likelihood may be higher, but the posterior distribution is essentiall uniform (the prior)

# even when letting b vary from 1.5 to 100, the distribution was
# about uniform (maybe slightly skewed), but no increase in likelihood

par(mfrow=c(1,1))
plot(params[,c(1,2)])
plot(params[,c(2,3)])
plot(params[,c(1,3)])
