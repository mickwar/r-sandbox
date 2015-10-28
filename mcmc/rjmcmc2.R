### Libraries
library(MASS)

### The data
dat = read.table("~/files/data/MPG_saturn.txt", header = TRUE)
y = dat$miles / dat$gallons
n = length(y)

plot(density(y))


### Probability of selecting model j from model i: P[i,j]
### In this case, there is 1/3 probability of staying in the same model,
### 1/3 of going down one dimension, and 1/3 of going up one dimension
P = diag(n)                     # add 1 to diagonal
P[(which(P == 1) + 1)[-n]] = 2  # get index of diagonal make one lower a 2
P[(which(P == 1) - 1)[-1]] = 3  # make one greater than diagonal a 3
P = ifelse(P > 0, 1/3, 0)       # now set all non-zero elements to 1/3
P[1,1] = P[1,2] = P[n,n] = P[n-1,n] = 1/2

P[1:6, 1:6]

# Note: since I use the sample() function to propose a possible dimension jump,
# I don't actually need to make sure these rows sum to 1 since sample() takes
# care of normalizing. I could just as well keep all the non-zero entries to
# be the same value and it would still correspond to have 1/3 probability for
# birth, stay, and death steps


### Function for computing the (proportional) posterior
# modelk denotes the number of dimensions to compute the posterior at
# paramsk are the parameters associated with modelk
calc.post = function(modelk, paramsk){
    out = 0

    # Fixing the weights, based on number of dimensions
    weights = ((1 - prior.alpha)^(1:modelk - 1) * prior.alpha) / (1 - (1 - prior.alpha)^modelk)
    # Consider using an actual prior for the weights (say Dirichlet distribution)
    # Note as modelk -> infinity, weights -> geometric(alpha), so weights
    # come from the "truncated geometric"

    # Likelihood
    # Odd indices are the means, even are the variances, so elements 3 and 4 is the mean
    # and variance of model 2
    temp = matrix(paramsk, ncol = 2, byrow = TRUE)
    p.means = temp[,1]
    p.vars = temp[,2]
    for (i in 1:n)
        out = out + log(sum(weights * dnorm(y[i], p.means, sqrt(p.vars))))
#       out = out + log(sum(dnorm(y[i], head(paramsk, modelk), sqrt(tail(paramsk, modelk)))))
#       # under the old setup


    # Priors
    # Geometric prior with parameter prior.k
    out = out + dgeom(modelk, prior.k, log = TRUE)

    # iid Normals for the mu vector parameter
    out = out + sum(dnorm(p.means, prior.mu.mean, prior.mu.sd, log = TRUE))

    # iid Gamms for the sigma^2 vector parameter
    out = out + sum(dgamma(p.vars, prior.sig2.a, prior.sig2.b, log = TRUE))

    return (out)
    }

prior.alpha = 0.5 # Prior for weights (though weights is fixed)

prior.k = 0.5 # Prior for K

prior.mu.mean = 28 # about mean(y), so same mean for -all- dimensions
prior.mu.sd = 10 # about 2*sd(y)

prior.sig2.a = 5 # mean of 5, variance of 5
prior.sig2.b = 1

### MCMC settings
nburn = 100
nmcmc = 100

# List of matrices containing the parameters. Since I know the maximum number of dimensions,
# I can pre-allocate space for the data, though most will remain unused. This could probably
# be improved
params = list()
for (i in 1:n)
    params[[i]] = matrix(0, nburn + nmcmc, i*2)

# Updating in blocks so only looking at one acceptance rate
accept = double(nburn + nmcmc)

# Initial values
kcurr = 1       # Starting at model 1 (i.e. a single univariate normal distribution)

# Set the means and variances for all mixture components in all models to data mean and variance
for (i in 1:n)
    params[[i]][1,] = rep(c(mean(y), var(y)), times = i)

# To keep track of which model I'm in at each iteration
model.path = double(nburn + nmcmc)
model.path[1] = kcurr

# The vector of indices marking which iteration was drawn last each model
model.iter = rep(1, n)


# Make sure this actually evalutes
calc.post(kcurr, params[[kcurr]][1,])


# MCMC loop
for (i in 2:(nburn + nmcmc)){
    knext = sample(1:n, 1, prob = P[kcurr,])

    # Birth
    if (knext == kcurr + 1){
        cand = c(params[[kcurr]][i-1,],
            mvrnorm(1, tail(params[[knext]][model.iter[knext],], 2), 1/knext*diag(2)))
        if (tail(cand, 1) > 0){
            if (log(runif(1)) < calc.post(knext, cand) - calc.post(kcurr, params[[kcurr]][i-1,]) +
                log(P[knext, kcurr]) - log(P[kcurr, knext] ) -
                sum(dnorm(tail(cand, 2), tail(params[[knext]][model.iter[knext],], 2),
                    1/knext, log = TRUE))){
                kcurr = knext
                accept[i] = 1
                params[[kcurr]][i,] = cand
                }
            }
        knext = -Inf
        }

    # Stay (standard Metropolis-Hastings step)
    if (knext == kcurr){
        params[[kcurr]][i,] = params[[kcurr]][i-1,]
        cand = mvrnorm(1, params[[kcurr]][i,], 1/kcurr*diag(kcurr*2))   
        if (all(cand[seq(2, 2*kcurr, by = 2)] > 0)){
            if (log(runif(1)) < calc.post(kcurr, cand) - calc.post(kcurr, params[[kcurr]][i,])){
                params[[kcurr]][i,] = cand
                accept[i] = 1
                }
            }
        knext = -Inf
        }

    # Death
    if (knext == kcurr - 1){
        }

    model.path[i] = kcurr
    model.iter[kcurr] = i
    }
