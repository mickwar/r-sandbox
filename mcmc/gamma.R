# example code using the bayes functions file
# using bayesian regression
source("./bayes_functions.R")

# generate some data
set.seed(1)
n = 5000
true.beta = c(2, 0.5, -1)
true.a = 0.6
true.b = 2.5
dat.x = matrix(runif(n*2, -5, 10), n, 2)
X = cbind(1, dat.x)
Y = X %*% true.beta + rgamma(n, true.a, scale=true.b)

# calculatle log posterior
ones = t(rep(1, n))
calc.post = function(params){
    # log-likelihood
    out = -n * lgamma(params[ind.a]) - params[ind.a] * n * log(params[ind.b]) +
        (params[ind.a] - 1) * ones %*% log(Y + params[ind.a] * params[ind.b] -
        X %*% params[ind.beta]) - 1/params[ind.b] * ones %*% (Y + params[ind.a] *
        params[ind.b] - X %*% params[ind.beta])
    # priors
    # gamma on a and b
    out = out + (a.a-1)*log(params[ind.a])-params[ind.a]/b.a
    out = out + (a.b-1)*log(params[ind.b])-params[ind.b]/b.b
    # normal on beta
    out = out - sum(1/2*log(v.beta) - 1/(2*v.beta)*(params[ind.beta]-m.beta)^2)
    return (out)
    }

ind.a = 1
ind.b = 2
ind.beta = 3:5
nparams = 5

# support bounds for each parameter
lower = double(nparams)-Inf
upper = double(nparams)+Inf
lower[c(ind.a, ind.b)] = 0

# hyperparameter specifications
# a
a.a = 2
b.a = 1
# b
a.b = 2
b.b = 1
# beta
m.beta = c(0, 0.5, -1.0)
v.beta = 10

nburn = 10000
nmcmc = 20000
params = matrix(0, nburn+nmcmc, nparams)
# starting values
params[1, ind.a] = abs(min(Y))
params[1, ind.b] = 1.1
params[1, 3] = 0
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
window = 200

# mcmc loop
for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in c(1,2,4,5)){
        cand = rnorm(1, params[i,j], sigs[j])
        cand.param[j] = cand
        if (cand >= lower[j] && cand <= upper[j] &&
            all(Y + prod(cand.param[1:2]) - X %*% cand.param[3:5] > 0)){
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
        cat(i, "\n")
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

plot(params[,1], type='l')
plot(params[,2], type='l')
plot(params[,3], type='l')
plot(params[,4], type='l')
plot(params[,5], type='l')

apply(accept, 2, mean)
apply(params, 2, mean)

# loop the posterior density of each parameter
names.VAR = c("a", "b", "int", "b1", "b2")
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

