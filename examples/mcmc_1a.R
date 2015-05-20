### "Full" Method that calculates full posterior each update

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
    # normal on b0 and b1
    out = out-1/2*log(sig_0)-1/(2*sig_0)*(beta[1]-mu_0)^2
    out = out-1/2*log(sig_1)-1/(2*sig_1)*(beta[2]-mu_1)^2
    return (out)
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

# initialize log posterior value
post = calc.post(params[1,])

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
