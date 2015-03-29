dat = read.table("~/files/data/room_creak.txt")[,1]

plot(dat, type='l')

plot(density(dat))

# gamma moment matching
curve(dgamma(x, shape = mean(dat)^2 / var(dat), scale = var(dat) / mean(dat)),
    from = 0, to = 25, add = TRUE, col = 'red')

# normal moment matching
curve(dnorm(x, mean(dat), sd(dat)), from = 0, to = 25,
    add = TRUE, col = 'blue')

plot(density(dat))
curve(dgamma(x, 1, 1), from=0, to=25, add=TRUE, col='green')

### mcmc
source("~/files/R/mcmc/bayes_functions.R")

calc.post = function(params){
    a = params[1]
    s = params[2]
    sum(dgamma(dat, shape=a, scale=s, log=TRUE)) + 
        dgamma(a, shape=5, scale=1, log=TRUE) + 
        dgamma(s, shape=1, scale=1, log=TRUE)
    }

nburn = 50000
nmcmc = 100000
window = 500
nparams = 2
params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
lower = rep(0, nparams)
upper = rep(Inf, nparams)
params[1,] = c(5, 1)

sigs = rep(1, nparams)

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

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]
post.mat = post.mat[(nburn+1):(nburn+nmcmc),]

apply(accept, 2, mean)
plot(params[,1], type='l')
plot(params[,2], type='l')
plot(params, type='l')

est.mode = function(params, post){
    # params, post are nparams x nmcmc
    # the first nparams - 1 columns in the first row
    # of post should be -Inf
    nparams = nrow(params)
    maxx = which.max(as.numeric(post))
    end = 1 + ((maxx - 1) %% nparams) # the index of the ending parameter
    order = c((nparams - end + 1):nparams, 1:(nparams - end))[1:nparams]
    vec = as.numeric(params)[(maxx-nparams+1):maxx]
    return(list("mode"=vec[order], "height"=post[maxx]))
    }
mode = est.mode(t(params), t(post.mat))
means = apply(params, 2, mean)


post.pred = rgamma(nmcmc, shape = params[,1], scale = params[,2])
plot(density(dat))
lines(density(post.pred), col='green')

mean(post.pred)
var(post.pred)

plot(density(params[,1] * params[,2]))
plot(density(params[,1] * params[,2]^2))
