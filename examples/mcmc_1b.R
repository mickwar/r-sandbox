### "Likelihood" Method that calcualtes priors and likelihood seperately

### generate some data
set.seed(1)
n = 10
datx = runif(n, 0, 10)
b0 = 2
b1 = 1.2
sig2 = 1
daty = b0 + b1*datx + rnorm(n, 0, sqrt(sig2))

# into vectors
y = as.matrix(daty)
x = cbind(1, datx)

### calculate log likelihood
calc.like = function(params){
    beta = params[ind.beta]
    sig2 = params[ind.eps]

    # log-likelihood
    out = -n/2 * log(sig2) - 1/(2*sig2) *
        t(y - x %*% beta) %*% (y - x %*% beta)

    return (out)
    }
### calculate log priors
calc.prior = function(param, j){
    if (j == 1)
        return (-1/2*log(sig_0)-1/(2*sig_0)*(param-mu_0)^2)
    if (j == 2)
        return (-1/2*log(sig_1)-1/(2*sig_1)*(param-mu_1)^2)
    if (j == 3)
        return((eps.a-1)*log(param)-param/eps.b)
    }

# hyperparameter specifications
eps.a = 2
eps.b = 10
mu_0 = 0
mu_1 = 0
sig_0 = 100
sig_1 = 100

# indices for parameters
ind.beta = 1:2
ind.eps = 3
nparams = 3

# support bounds for each parameter
lower = double(nparams)-Inf
upper = double(nparams)+Inf
lower[ind.eps] = 0


### MCMC settings
nmcmc = 10000
params = matrix(0, nmcmc, nparams)
accept = matrix(0, nmcmc, nparams)
# starting values
params[1, ind.eps] = 2.5
sigs = rep(1, nparams)

# initialize log likelihood and log prior values
like = calc.like(params[1,])
priors = double(nparams)
for (i in 1:nparams)
    priors[i] = calc.prior(params[1,i], i)

# vector used to testing proposals
cand.param = params[1,]


### MCMC loop
for (i in 2:nmcmc){
    # set current parameters to previous draws
    params[i,] = params[i-1,]

    # iterate through each parameter
    for (j in 1:nparams){
        # draw a candidate value from proposal distribtion (all Normal in this case)
        cand = rnorm(1, params[i,j], sigs[j])

        # check if candidate is in proper bounds
        if (cand >= lower[j] && cand <= upper[j]){
            cand.param[j] = cand
            cand.like = calc.like(cand.param)
            cand.prior = calc.prior(cand.param[j], j)
            # check whether to accept draw or not
            if (log(runif(1)) < cand.like + cand.prior - like - priors[j]){
                like = cand.like        # store likelihood
                priors[j] = cand.prior  # store prior
                params[i,j] = cand      # update parameter draws
                accept[i,j] = 1         # indicator for accepting
            } else {
                # candidate not accepted, replace cand.param[j] with previous draw
                cand.param[j] = params[i,j]
                }
        } else {
            # candidate out of range, replace cand.param[j] with previous draw
            cand.param[j] = params[i,j]
            }
        }

    }
