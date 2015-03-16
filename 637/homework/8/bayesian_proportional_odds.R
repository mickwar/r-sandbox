library(splines)
source("~/files/R/mcmc/bayes_functions.R")

load("~/files/data/casella.Rdata")

### training data
n = NROW(train)
Y = train$stars.scores
J = max(Y)

X = cbind(1, train$help.yes/train$help.total)
X[is.na(X)] = 0
X = cbind(1, 100*X[,2], nchar(as.character(train$review)))
#    ns(nchar(as.character(train$review)), knots = c(1000)))
p = ncol(X)

plot(nchar(as.character(train$review)), Y, pch=20)

### test data
n.test = NROW(test)
Y.test = test$stars.scores

X.test = cbind(1, test$help.yes/test$help.total)
X.test[is.na(X.test)] = 0
X.test = cbind(1, 100*X.test[,2], nchar(as.character(test$review)))
#    ns(nchar(as.character(test$review)), knots = c(1000)))

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


nburn = 10000
nmcmc = 20000
window = 500
# nparams = (J-1) * p # cumulative logit
# nparams = (J-1) + (p-1) # porportional
nparams = (J-1) + (p-1) # proportional (with spline)
params = matrix(0, nburn + nmcmc, nparams)
accept = matrix(0, nburn + nmcmc, nparams)
lower = rep(-Inf, nparams)
upper = rep(Inf, nparams)
params[1, 1:4] = seq(0, 1, length=4)

sigs = rep(0.1, nparams)
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

plot(post.mat[,nparams], type='l')

hpds = apply(params, 2, hpd.uni)
plot(params[,1], type='l', col=1, ylim=range(params))
plot(params[,6], type='l', col=1)#, ylim=range(params))
for (i in 1:nparams)
    lines(params[,i], type='l', col=i)

pairs(params[1:1000,], pch=20)


pdf("./figs/post.pdf", height = 9, width = 9)
main.vec = c("Int1", "Int2", "Int3", "Int4", "Helpful", "Length of Review")
par(mfrow=c(3,2), mar = c(2.1, 4.1, 4.1, 2.1))
plot.post(params[,1], density(params[,1]), hpds[,1], main = main.vec[1], cex = 2.5, xlab = "")
par(mfg=c(1,2,3,2))
plot.post(params[,2], density(params[,2]), hpds[,2], main = main.vec[2], cex = 2.5, xlab = "")
par(mfg=c(2,1,3,2))
plot.post(params[,3], density(params[,3]), hpds[,3], main = main.vec[3], cex = 2.5, xlab = "")

par(mfg=c(2,2,3,2))
plot.post(params[,4], density(params[,4]), hpds[,4], main = main.vec[4], cex = 2.5, xlab = "")
par(mfg=c(3,1,3,2))
plot.post(params[,5], density(params[,5]), hpds[,5], main = main.vec[5], cex = 2.5, xlab = "")
par(mfg=c(3,2,3,2))
plot.post(params[,6], density(params[,6]), hpds[,6], main = main.vec[6], cex = 2.5, xlab = "")
dev.off()




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
    out = matrix(0, nrow(params), length(Y))
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
#pred = post.pred(params, list("Y"=Y, "X"=X))
mse.hpd = hpd.uni(pred$mse)
calc.mode(density(pred$mse))
mean(pred$mse)

pdf("figs/mse.pdf")
plot.post(pred$mse, density(pred$mse), mse.hpd, main = "Posterior MSE", xlab = "MSE")
dev.off()



# hist(pred$preds[,20])


# library(VGAM)
# mod = vglm(Y ~ X[,2:4], family = cumulative(parallel = TRUE))
# summary(mod)
# 
# plot(params[,c(5,6)], type='l')
# 
# apply(params, 2, mean)
# apply(params, 2, quantile, c(0.025, 0.975))
