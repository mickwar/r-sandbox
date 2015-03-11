source("~/files/R/mcmc/bayes_functions.R")

load("~/files/data/casella.Rdata")

### training data
new.train = train[,-which(names(train) %in% c("title.review", "review", "name"))]
new.train[is.na(new.train)] = 0

n = NROW(new.train)
Y = new.train$stars.scores
J = max(Y)
# temp = matrix(0, n, 5)
# for (i in 1:n)
#     temp[i, Y[i]] = 1
# Y = temp
# rm(temp)


plot(X[, 2], Y, pch=20, col=Y)



#X = cbind(1, new.train$help.yes, new.train$help.total, new.train$year)
X = cbind(1, train$help.yes, train$help.total, train$year)
X[is.na(X)] = 0
Xmean = apply(X, 2, mean)
Xsd = apply(X, 2, sd)
Xmean[4] = 2007
Xsd[4] = 1
X = (X - matrix(Xmean, nrow(X), length(Xmean), byrow = TRUE)) /
    matrix(Xsd, nrow(X), length(Xsd), byrow = TRUE)
X[,1] = 1

X = cbind(X, X[,2] * X[,3])
# X = cbind(1, train$help.yes, train$help.total, train$help.yes*train$help.total,
#     train$year, train$year^2)
p = ncol(X)


### test data
new.test = test[,-which(names(test) %in% c("title.review", "review", "name"))]
new.test[is.na(new.test)] = 0

n.test = NROW(new.test)
Y.test = new.test$stars.scores

#X.test = cbind(1, new.test$help.yes, new.test$help.total, new.test$year)
X.test = cbind(1, test$help.yes, test$help.total, test$year)
X.test[is.na(X.test)] = 0
X.test = (X.test - matrix(Xmean, nrow(X.test), length(Xmean), byrow = TRUE)) /
    matrix(Xsd, nrow(X.test), length(Xsd), byrow = TRUE)
X.test[,1] = 1
X.test = cbind(X.test, X.test[,2] * X.test[,3])

### compute log posterior
calc.post = function(params){
    # calculate p's sequentially
    temp.p = matrix(0, n, J)
    for (j in 1:(J-1)){
        ### cumulative logit
        # pindex = seq(p*(j-1)+1, j*p, 1)
        # temp.xp = exp(X %*% params[pindex])
        
        ### proportional odds
        temp.xp = exp(X %*% params[c(j, J:nparams)])
        if (j == 1){
            temp.p[,j] = temp.xp / (1 + temp.xp)
        } else {
            if (j > 2){
                temp.sum = apply(temp.p[,1:(j-1)], 1, sum)
            } else {
                temp.sum = temp.p[,1]
                }
            temp.p[,j] = (1 - temp.sum) * temp.xp - temp.sum
            temp.p[,j] = temp.p[,j] / (1 + temp.xp)
            }
       }
    temp.p[,J] = apply(temp.p[,1:(J-1)], 1, function(x) 1 - sum(x))
    px = diag(temp.p[1:n, Y])
    if (any(px <= 0 | px >= 1))
        return (-Inf)

    # likelihood
    out = sum(log(px))

    # priors
    out = out + sum(dnorm(params, 0, beta.sd, log = TRUE))
    return (out)
    }
beta.sd = 10


nburn = 5000
nmcmc = 10000
window = 500
# nparams = (J-1) * p # cumulative logit
nparams = (J-1) + (p-1) # porportional
params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
lower = rep(-Inf, nparams)
upper = rep(Inf, nparams)
params[1, 1:4] = seq(0, 1, length=4)

sigs = rep(0.001, nparams)
# sigs = c(1.2340859, 2.4084625, 1.9771273, 1.5098575, 0.6543709, 1.2939921, 1.1159381,
#     0.7045175, 0.5206961, 1.5336611, 1.3352886, 0.6416922, 0.6122230, 2.2898420,
#     1.9232403, 0.8846596)


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

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]
post.mat = post.mat[(nburn+1):(nburn+nmcmc),]
 
apply(accept, 2, mean)

# plot(c(t(post.mat)), type='l')
plot(post.mat[,nparams], type='l')
# 
plot(params[,1], type='l', col=1, ylim=range(params))
for (i in 1:nparams)
    lines(params[,i], type='l', col=i)
# 
# plot(params[,c(1,5)], type='l')
# pairs(params[1:1000,], pch=20)

### estimate the mode
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

### posterior predictions 
post.pred = function(params, data){
    Y = data$Y
    X = data$X
    n = nrow(X)
    mse = double(nrow(params))
    out = matrix(0, nrow(params), length(Y.test))
    for (i in 1:nrow(params)){
        temp.p = matrix(0, n, J)
        for (j in 1:(J-1)){
            ### proportional odds
            temp.xp = exp(X %*% params[i, c(j, J:nparams)])
            if (j == 1){
                temp.p[,j] = temp.xp / (1 + temp.xp)
            } else {
                if (j > 2){
                    temp.sum = apply(temp.p[,1:(j-1)], 1, sum)
                } else {
                    temp.sum = temp.p[,1]
                    }
                temp.p[,j] = (1 - temp.sum) * temp.xp - temp.sum
                temp.p[,j] = temp.p[,j] / (1 + temp.xp)
                }
           }
        temp.p[,J] = apply(temp.p[,1:(J-1)], 1, function(x) 1 - sum(x))
        # the posterior predictive
        out[i,] = apply(temp.p, 1, function(x) sample(1:J, 1, prob = x))
        mse[i] = mean((out[i,] - Y)^2)
#       mse[i] = mean(abs(out[i,] - Y))
        }
    return (list("preds"=out, "mse"=mse))
    }
pred = post.pred(params, list("Y"=Y.test, "X"=X.test))
pred = post.pred(params, list("Y"=Y, "X"=X))
plot(density(pred$mse))
calc.mode(density(pred$mse))

hist(pred$preds[,12])


library(VGAM)
mod = vglm(Y ~ X[,2:4], family = cumulative(parallel = TRUE))
summary(mod)

plot(params[,c(5,6)], type='l')

apply(params, 2, mean)
apply(params, 2, quantile, c(0.025, 0.975))
