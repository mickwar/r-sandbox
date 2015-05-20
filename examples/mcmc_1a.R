### "Full" Method that calculates full posterior each update

### generate some data
set.seed(1)
n = 200
p = 12
x = cbind(1, matrix(runif(n*p, 0, 10), n, p))
true.beta = seq(-5, 5, length = ncol(x))
true.sig2 = 1
y = x %*% true.beta + rnorm(n, 0, sqrt(true.sig2))

# hyperparameter specifications
eps.a = 2
eps.b = 10
mu_0 = 0
sig_0 = 100

# indices for parameters
ind.beta = 1:ncol(x)
ind.eps = ncol(x)+1
nparams = ncol(x)+1

# support bounds for each parameter
lower = double(nparams)-Inf
upper = double(nparams)+Inf
lower[ind.eps] = 0

### calculatle log posterior
calc.post = function(params){
    beta = params[ind.beta]
    sig2 = params[ind.eps]

    # log-likelihood
    out = -n/2 * log(sig2) - 1/(2*sig2) *
        t(y - x %*% beta) %*% (y - x %*% beta)

    # priors
    # gamma on sigma^2
    out = out + (eps.a-1)*log(sig2)-sig2/eps.b
    # normal on betas (iid)
    out = out + sum(-1/2*log(sig_0)-1/(2*sig_0)*(beta-mu_0)^2)
    return (out)
    }


### MCMC settings
nmcmc = 100000
params = matrix(0, nmcmc, nparams)
accept = matrix(0, nmcmc, nparams)

# candidate sigmas
sigs = double(nparams)
sigs[ind.beta] = 0.05
sigs[ind.eps] = 0.5

# starting values (starting at true values)
params[1, ind.beta] = true.beta
params[1, ind.eps] = true.sig2

# initialize log posterior value
post = calc.post(params[1,])

# vector used to testing proposals
cand.param = params[1,]


### MCMC loop
system.time({
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
            cand.post = calc.post(cand.param)
            # check whether to accept draw or not
            if (log(runif(1)) < cand.post - post){
                post = cand.post        # store post
                params[i,j] = cand      # store parameter draws
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
})
