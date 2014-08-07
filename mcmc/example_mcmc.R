# example code using the bayes functions file
# using bayesian regression
source("./bayes_functions.R")

# generate some data
set.seed(1)
n = 12
datx = runif(n, 0, 10)
b0 = 2
b1 = 1.2
sig2 = 1
daty = b0 + b1*datx + rnorm(n, 0, sqrt(sig2))
plot(datx, daty, pch=20)

# into vectors
y = as.matrix(daty)
x = cbind(1, datx)

# calculatle log posterior
calc.post = function(params){
    # log-likelihood
    out = -n/2 * log(params[ind.eps]) - 1/(2*params[ind.eps]) *
        t(y - x %*% params[ind.beta]) %*% (y - x %*% params[ind.beta])
    # priors
    # gamma on sigma^2
    out = out + (eps.a-1)*log(params[ind.eps])-params[ind.eps]/eps.b
    # normal on b0 and b1
    out = out-1/2*log(sig_0)-1/(2*sig_0)*(params[ind.beta][1]-mu_0)^2
    out = out-1/2*log(sig_1)-1/(2*sig_1)*(params[ind.beta][2]-mu_1)^2
    return (out)
    }

ind.beta = 1:2
ind.eps = 3
nparams = 3

# support bounds for each parameter
lower = double(nparams)-Inf
upper = double(nparams)+Inf
lower[ind.eps] = 0

# hyperparameter specifications
eps.a = 2
eps.b = 10
mu_0 = 0
mu_1 = 0
sig_0 = 100
sig_1 = 100

nburn = 10000
nmcmc = 20000
params = matrix(0, nburn+nmcmc, nparams)
# starting values
params[1, ind.eps] = 2.5
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
window = 100

# mcmc loop
dir = "ex_mcmc"
prefix = "mcmc_"
mcmc_time(iter = 0, dir = dir, prefix = prefix)
for (i in 2:(nburn+nmcmc)){
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
        }
        # adjust the candidate sigma
    if (floor(i/window) == i/window && i <= nburn)
        sigs = sigs*autotune(apply(accept[(i-window+1):i,], 2,
            mean), k = max(window/50, 1.1))
    mcmc_time(do = TRUE, iter = i, every = 100, params, accept,
        sigs, nburn, nmcmc, dir = dir, prefix = prefix)
    }





params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

# check acceptance rates
apply(accept, 2, mean)

# loop the posterior density of each parameter
names.VAR = c("b0", "b1", "sigma")
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
    # pause and wait for user input (hit enter)
    if (i < nparams)
        readline()
    }
 
