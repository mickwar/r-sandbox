source("~/files/R/mcmc/bayes_functions.R")

load("~/files/data/casella.Rdata")

str(test)

str(train)

train[1,]

new.train = train[,-which(names(train) %in% c("title.review", "review", "name"))]
new.train[is.na(new.train)] = 0

n = NROW(new.train)
Y = new.train$stars.scores
temp = matrix(0, n, 5)
for (i in 1:n)
    temp[i, Y[i]] = 1
Y = temp
rm(temp)

X = cbind(1, new.train$help.yes, new.train$help.total, new.train$year)


calc.post = function(param.vec, like.only = FALSE){
    # likelihood
    px = pnorm(X %*% param.vec)
    if (any(px == 0 | px == 1)) # when likelihood is 0
        return (-Inf)
    out = sum(y*log(px) + (1-y)*log(1-px))
    # priors (beta ~ N(0, 1))
    if (!like.only)
        out = out + sum(dnorm(param.vec, beta.mean, beta.var, log = TRUE))
    return (out)
    }

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
