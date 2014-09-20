source("~/files/R/mcmc/bayes_functions.R")
library(splines)
movie.dat = read.table("../data/movie_data.txt", header=T)
ben = read.table("../data/movie_ben.txt", header=T)
ben[,1] = as.character(ben[,1])
movies = ben[,1]


movie.dat[,1] = as.character(movie.dat[,1])
movie.dat$MPAA = as.character(movie.dat$MPAA)

n.var = ncol(movie.dat) - 1
dat = data.frame(matrix(0, nrow(ben), n.var + ncol(ben)))
dat[,1:2] = ben[,2:1]
not_found = NULL
for (i in 1:nrow(ben)){
    j = which(ben[i,1] == movie.dat[,1])
    if (length(j) == 0){
        not_found[length(not_found)+1] = ben[i,1]
    } else {
        dat[i,3:(2+n.var)] = movie.dat[j,2:(n.var+1)]
        }
    }

names(dat) = c("ben", names(movie.dat))

dat = cbind(dat, "kids_I"=is.na(dat$kids_S)*1)
for (i in 3:5)
    dat[,i] = ifelse(is.na(dat[,i]), 0, dat[,i])

# functions
# if drop = TRUE, drop the last column
mpaa.matrix = function(x){
    level = c("PG13", "PG", "G")
    # NR is the dropped level
    n = length(x)
    k = length(level)
    out = matrix(0, n, k)
    for (i in 1:k)
        out[which(x == level[i]), i] = 1
    return (out)
    }
mw.smooth = function(x, y, m, d){
    # m = number of points to do the smoothing,
    # should be greater than length(x) to look nice
    if (missing(m))
        m = max(100, 8*length(x))
    if (missing(d))
        d = diff(range(x)) / length(x)

    # gaussian kernel, delta is the bandwidth parameter
    kern = function(x, y)
        exp(-1/(2*d^2)*(x-y)^2)

    # go outside the range of the data by some
    w = diff(range(x))*0.25

    outx = seq(min(x) - w, max(x) + w, length=m)
    outy = double(m)
    # loop to compute the weighted value at each new x
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i]))/
            sum(kern(x, outx[i]))
    return (list("x"=outx, "y"=outy))
    }
x.spline = function(x, var, at){
    if (var == 1){
        out = cbind(ns(x[,var], knots = at), x[,2:ncol(x)])
    } else {
        if (var == ncol(x) && var != 1){
            out = cbind(x[,1:(ncol(x)-1)], ns(x[,var], knots = at))
        } else {
            out = cbind(x[,1:(var-1)], ns(x[,var], knots = at),
                x[,(var+1):ncol(x)])
            }
        }
    return (out)
    }

# make the X matrix
Y = dat[,1]
raw.X = as.matrix(cbind(dat[,3:5], mpaa.matrix(dat$MPAA),
    dat[,7:15]))
X = cbind(raw.X[,1:6], 100*(raw.X[,7]-raw.X[,9])/(raw.X[,7]+
    raw.X[,8]-raw.X[,9]-raw.X[,10]), 100*raw.X[,9]/(raw.X[,9]+ 
    raw.X[,10]), raw.X[,11:15])

# top critic rotten tomatoes indicator
X = cbind(X, "rot_I"=is.na(X[,8])*1)
X[,8] = ifelse(is.na(X[,8]), 0, X[,8])

var.names = c(names(dat)[3:5], as.character(unique(dat[,6])[1:3]),
    "rot_crit", "rot_top", names(dat)[11:15], "rot_I")
colnames(X) = var.names

n = nrow(X)

CC = X[,which(colnames(X) == "kids_P")]
CC = X[,3]
SS = sort(unique(CC))
means = double(length(SS))
for (i in 1:length(means))
    means[i] = mean(Y[which(CC == SS[i])])
smooth = mw.smooth(SS, means)
plot(CC+rnorm(n,0,0.05),dat$ben+rnorm(n,0,0.010),pch=20,cex=0.5)
lines(smooth$x, smooth$y, col='red')

X2 = x.spline(X, 2, c(2, 6))
X2 = x.spline(X2, 5, 3)

X = X2

# removing the correlated variables
X.un = data.frame(X[,c(1,2,3,4,5,6,7,10,12,13,14)])
X.un[,9] = log(X.un[,9])

full.mod = glm(Y ~ kids_S + kids_V + ns(kids_P, 3) + PG13 + PG + G +
    rot_crit + live_action + kids_I + rot_I,
    data=X.un, family=binomial)
null.mod = glm(Y ~ 1, data=X.un, family=binomial)

mod = step(null.mod, scope=list(lower=null.mod, upper=full.mod),
    k=log(nrow(dat)), data=dat, direction="both", family=binomial)
summary(mod)
vif(mod)
vif(full.mod)

risk.score = predict(mod, newdata=X.un)
p = exp(risk.score)/(1+exp(risk.score))

cbind(substr(movies, 1, 15), p, Y, abs(p-Y))

those = as.vector(which(p > 0.83))
cbind(substr(movies[those], 1, 15), p[those], Y[those], abs(p[those]-Y[those]))

cbind(movies[those], X.un[those,])




# the model: y_i ~ bern(p_i)
#            logit(p_i) = x_i * beta
#            beta ~ n(0, sig^2 * g * (x'x)^(-1))             
#            fix g = n
#            sig^2 ~ gamma(a, b)
# log multivariate normal density (kernel) (sig is the cholesky 
# decomposition of the covariance matrix, i.e. sig = chol(sigma))
dmvnorm = function(x, mu, sig)
    -sum(log(diag(sig))) - 1/2*t(x-mu) %*% chol2inv(sig) %*% (x - mu)
# function for log posterior
calc.post = function(params){
    p = 1 / (1 + exp(-X %*% params[ind.beta]))
    # likelihood
    out = sum(dbinom(Y, 1, p, log=TRUE))
    # priors
    # multivariate normal for betas
    out = out + dmvnorm(params[ind.beta], mu, sig_chol)
    # gamma for sigma^2
    out = out + (sig.a-1)*log(params[ind.sig]) - params[ind.sig]/sig.b
    return (out)
    }

ind.beta = 1:ncol(X)
ind.sig = max(ind.beta) + 1
nparams = ind.sig

n = nrow(X)
# train.obs = sort(sample(n, floor(n*1), replace=FALSE))
# X = full.X[train.obs,]
# #test.X = full.X[(1:n)[-train.obs],]
# test.X = full.X
# train.Y = Y[train.obs]
# #test.Y = Y[(1:n)[-train.obs]]
# test.Y = Y
# train.movies = movies[train.obs]
# #test.movies = movies[(1:n)[-train.obs]]
# test.movies = movies

# fixed parameters
g = nrow(X)
mu = rep(0, ncol(X))
sig.a = 2
sig.b = 0.1

lower = rep(-Inf, nparams)
upper = rep(+Inf, nparams)
lower[nparams] = 0

window = 200
nburn = 5000
nmcmc = 10000

params = matrix(0, nburn+nmcmc, nparams)
params[1, ind.sig] = 80
accept = matrix(0, nburn+nmcmc, nparams)
sig_chol = chol(params[1, ind.sig] * g * solve(t(X) %*% X))
sigs = rep(1, nparams)
#sigs = c(0.4478, 0.1743, 0.3559, 1.0953, 1.4509, 1.4804, 0.0094,
#    0.0112, 0.0091, 0.8352, 0.0955, 2.1764, 4.3964, 0.1486)
cand.param = params[1,]
curr.post = calc.post(params[1,])

# metropolis loop
for (iter in 2:(nburn+nmcmc)){
    # copy previous vector of parameters to current iteration
    params[iter,] = params[iter-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[iter, j], sigs[j])
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            if (j == nparams)
                sig_chol = sig_chol * sqrt(cand / params[iter, j])
            cand.post = calc.post(cand.param)
            # test for acceptance
            if (log(runif(1)) < cand.post - curr.post){
                accept[iter, j] = 1
                params[iter, j] = cand
                curr.post = cand.post
            } else {
                if (j == nparams)
                    sig_chol = sig_chol * sqrt(params[iter, j] / cand)
                cand.param[j] = params[iter, j]
                }
            }
        # make adjustments to candidate sigmas
        if (floor(iter/window) == iter/window  && iter <= nburn)
            sigs[j] = sigs[j] * autotune(mean(accept[(iter-window+1):
                iter, j]), k = max(window/50, 1.1))
        }
    }

# burn-in
bparams = params[(nburn+1):(nburn+nmcmc),]
baccept = accept[(nburn+1):(nburn+nmcmc),]

# thin
thin = seq(1, nmcmc, by = 100)
bparams = bparams[thin,]
baccept = baccept[thin,]

# acceptance rates
apply(baccept, 2, mean)

# trace plots
for (j in 1:nparams){
    plot(bparams[,j], type='l', main=var.names[j])
    readline()
    }

# posterior densities
for (j in 1:nparams){
    dens = density(bparams[,j])
    hpd.int = hpd.mult(bparams[,j], dens)
    dens.col = "gray80"
    if (j == nparams){
        plot(dens, main=var.names[j], ylab="Density",
            xlab=expression(sigma^2))
    } else {
        plot(dens, main=var.names[j], ylab="Density",
            xlab=expression(beta))
        }
    polygon(dens, col=dens.col)
    shade = col.mult(dens.col, "cyan")
    for (k in 1:(length(hpd.int)/2))
        color.den(dens, hpd.int[2*k-1], hpd.int[2*k], shade)
    lines(dens, lwd=1.5)
    lines(range(dens$x), c(0, 0))
    lines(rep(bound(0, dens), 2), c(0, bound(0,
        dens, FALSE)), col='black', lwd=2, lty=2)
    readline()
    }

# posterior predictive
# risk = X %*% t(bparams[,-14])
# p = 1/(1+exp(-risk))
# hist(p[2,])
# apply(p, 1, function(x) calc.mode(density(x), "density"))

risk = test.X %*% t(bparams[,-14])
p = 1/(1+exp(-risk))
point.p = matrix(rbinom(length(p), 1, p), nrow(p), ncol(p))

pos.rate = double(ncol(point.p))
neg.rate = double(ncol(point.p))
error = double(ncol(point.p))
for (i in 1:ncol(point.p)){
    t.pos = sum(test.Y == 1 & point.p[,i] == 1)
    t.neg = sum(test.Y == 0 & point.p[,i] == 0)
    total.neg = length(test.Y)-sum(test.Y)
    total.pos = sum(test.Y)

    pos.rate[i] = 1-t.pos/total.pos
    neg.rate[i] = 1-t.neg/total.neg
    error[i] = mean(point.p[,i] != test.Y)
    }

hist(pos.rate, col='gray', breaks=50)
hist(neg.rate, col='gray', breaks=50)
hist(error, col='gray', breaks=50)

c(mean(pos.rate), calc.mode(density(pos.rate)))
c(mean(neg.rate), calc.mode(density(neg.rate)))
c(mean(error), calc.mode(density(error)))

# pred = read.table("../data/movie_pred.txt", header=T)
# pred.movies = pred[,2]
# pred = pred[,-2]
# pred = cbind(pred, "kids_I"=is.na(pred$kids_S)*1)
# for (i in 2:4)
#     pred[,i] = ifelse(is.na(pred[,i]), 0, pred[,i])
# pred = cbind(pred, "rot_I"=is.na(pred$rot_top)*1)
# for (i in 6:8)
#     pred[,i] = ifelse(is.na(pred[,i]), 0, pred[,i])
# 
# pred.x = as.matrix(cbind(pred[,2:4],mpaa.matrix(pred[,5]),pred[,6:12]))
# colnames(pred.x)[c(4,5,6)] = c("PG13", "PG", "G")
# pred.y = pred[,1]
# 
# pred.p = 1/(1+exp(- pred.x %*% t(bparams[,-14])))
# calc.mode(density(pred.p), "density")
