# example code using the bayes functions file

# Modeling data with a gamma and point mass mixture, where the
# the gamma is reparametrized to estimate the mean and variance:
# X ~ Gamma(a, b) is equivalent to Y ~ Gamma(mu^2/sig^2, sig^2/mu)
# E(X) = ab, Var(X) = ab^2, E(Y) = mu, Var(Y) = sig^2
# Define n-vector mu = exp(X %*% beta)
# beta ~ Normal(m, v)
# sig^2 ~ Gamma(a.sig, b.sig)

# incomplete

source("./bayes_functions.R")

# generate some data
set.seed(1)
n = 1000
true.beta = c(3, 0.1, -0.3, 0)
true.sig2 = 4
true.params = c(true.beta, true.sig2)
dat.x = matrix(runif(n*(length(true.beta)-1), -1, 14), n, length(true.beta)-1)
X = cbind(1, dat.x)
Y = rgamma(n, (exp(X %*% true.beta)^2)/true.sig2, scale = true.sig2/exp(X %*% true.beta))

plot(density(Y))

# make really small values 0
Y = ifelse(Y < 0.01, 0, Y)

range(Y)
plot(X[,2], Y, pch=20)
points(X[Y == 0, 2], Y[Y == 0], pch=20, col='blue')
plot(X[,3], Y, pch=20)
points(X[Y == 0, 3], Y[Y == 0], pch=20, col='blue')
plot(X[,4], Y, pch=20)
points(X[Y == 0, 4], Y[Y == 0], pch=20, col='blue')
pairs(cbind(Y, X))

bad.mod = lm(Y ~ dat.x)
plot(fitted(bad.mod), rstudent(bad.mod), pch=20)

# calculatle log posterior
calc.post = function(params){
    # log-likelihood
    mu = exp(X %*% params[ind.beta])
    a = (mu^2) / params[ind.sig]
    b = params[ind.sig] / mu
    out = sum(-lgamma(a) - a*log(b) + (a-1)*log(Y) - Y/b)
    # priors
    # gamma on sig2
    out = out + (a.sig-1)*log(params[ind.sig])-params[ind.sig]/b.sig
    # normal on beta
    out = out - sum(1/2*log(v.beta) - 1/(2*v.beta)*(params[ind.beta]-m.beta)^2)
    return (out)
    }

ind.beta = 1:4
ind.sig = 5
nparams = ind.sig

# support bounds for each parameter
lower = double(nparams)-Inf
upper = double(nparams)+Inf
lower[ind.sig] = 0

# hyperparameter specifications
# beta
m.beta = c(0, 0.5, -1.0, 0)
v.beta = 10
# sig
a.sig = 1.5
b.sig = 5
curve(dgamma(x, a.sig, scale = b.sig), from = 0, to = 10)

nburn = 50000
nmcmc = 50000
params = matrix(0, nburn+nmcmc, nparams)
# starting values
params[1, ind.sig] = 1
sigs = rep(1, nparams)

# if burn in has already been done (these should
# be the most up-to-date values)
#sigs = read.table("./mcmc_cand_sigmas.txt")[,1]
#params[1,] = read.table("./mcmc_init_params.txt")[,1]

# initialize log posterior value
post = calc.post(params[1,])
cand.param = params[1,]

# confidence intervals on acceptance rates?
accept = matrix(0, nburn+nmcmc, nparams)
window = 500

# mcmc loop
for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:nparams){
        cand = rnorm(1, params[i,j], sigs[j])
        cand.param[j] = cand
        if (cand >= lower[j] && cand <= upper[j]){
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
        }
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }

# remove burned-in
params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

for (i in 1:nparams){
    plot(params[,i], type='l'); abline(h = true.params[i], col='red', lwd=2)
    if (i < nparams)
        readline()
    }

# check for optimal acceptance rates
apply(accept, 2, mean)
# calculate means of posterior draws
apply(params, 2, mean)

# loop the posterior density of each parameter
names.VAR = c("b0", "b1", "b2", "b3", "sig2")
for (i in 1:nparams){
    dens = density(params[,i], width=sd(params[,i]))
    # compute the hpd set
    hpd.int = hpd.mult(params[,i], dens)
    plot(dens, main=names.VAR[i], ylab="Density")
    # overall fill color
    dens.col = "gray80"
    polygon(dens, col=dens.col)
    # mutiply cyan with gray80
    shade = col.mult(dens.col, "cyan")
    # shade in each portion of the hpd (if multiple intervals)
    for (k in 1:(length(hpd.int)/2))
        color.den(dens, hpd.int[2*k-1], hpd.int[2*k], shade)
    # re-draw the density
    lines(dens, lwd=1.5)
    lines(range(dens$x), c(0,0))
    # make a line at x=0
    at.x = 0
    lines(rep(bound(at.x, dens), 2), c(0, bound(at.x,
        dens, FALSE)), col='black', lwd=2, lty=2)
    # make a line at true value
    at.x = true.params[i]
    lines(rep(bound(at.x, dens), 2), c(0, bound(at.x,
        dens, FALSE)), col='red', lwd=2)
    # pause and wait for user input (hit enter)
    if (i < nparams)
        readline()
    }
