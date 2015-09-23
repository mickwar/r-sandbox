# Robert and Casella Example 11.4. A linear Jacobian.
# This is Green (1995)'s toy example (with my made up data)
set.seed(1)

n = 15
K = 2
sig = 1
x = runif(n, 0, 2)
y = rnorm(15, 1 + 0.2 * x, sig)

# Number of models
plot(x, y, pch = 20)

calc.post = function(k, theta){
    if (k == 1){
        # likelihood
        out = sum(dnorm(y, theta[1], sig, log = TRUE))
        # priors
        out = out + log(1/K) # prior probability that the model is K
        out = out + dnorm(theta[1], 0, 10, log = TRUE)
        }
    if (k == 2){
        # likelihood
        out = sum(dnorm(y, theta[1] + x*theta[2], sig, log = TRUE))
        # priors
        out = out + log(1/K) # prior probability that the model is K
        out = out + dnorm(theta[1], 0, 10, log = TRUE)
        out = out + dnorm(theta[2], 0, 10, log = TRUE)
        }
    return (out)
    }

# Transformation function for moving from model kcurr to kprime
# (kcurr is known implicitly in this case)
Tr = function(kprime, theta, u){
    if (kprime == 1)
        return (c((theta[1] + theta[2]) / 2, (theta[1] - theta[2]) / 2))
    if (kprime == 2)
        return (c(theta - u, theta + u))
    }

Pi = matrix(0.5, K, K)
# {Pi}_ij = the probability of choosing a jump to model j while in model i
# The rows must sum to 1

nburn = 0
nmcmc = 5000

params = list()
params[[1]] = matrix(0, nburn + nmcmc, 1)
params[[2]] = matrix(0, nburn + nmcmc, 2)
accept = matrix(0, nburn + nmcmc, K)
count = rep(1, K)

params[[1]][1,] = mean(y)
params[[2]][1,] = coef(lm(y ~ x))

kcurr = 1
calc.post(1, params[[1]][1,])
calc.post(2, params[[2]][1,])

for (i in 2:(nburn + nmcmc)){
    # Algorithm A.48 Green's Algorithm

    # 1. Select model M_kprime with probability {Pi}_{kcurr,kprime}
    kprime = sample(K, 1, prob = Pi[kcurr,])
    count[kprime] = count[kprime] + 1

    # 2. Generate u ~ g(u), conditioned on kcurr and kprime
    # In this case, only need to generate u if jumping from 
    # model 1 to model 2
    if ((kcurr == 1) && (kprime == 2)){
        u = rnorm(1, 0, 0.1)
    } else {
        u = NULL
        }


    # Proposing candidates
    if (kcurr == kprime){ # Not jumping between models, then this is just regular M-H
        cand = rnorm(NCOL(params[[kcurr]]), params[[kcurr]][count[kcurr]-1,], 0.1)
#       calc.post(kprime, cand)
#       calc.post(kcurr, params[[kcurr]][count[kcurr],])
        if (log(runif(1)) <= calc.post(kprime, cand) -
            calc.post(kcurr, params[[kcurr]][count[kcurr]-1,])){
            params[[kprime]][count[kcurr],] = cand
            accept[count[kcurr],kprime] = 1
        } else {
            params[[kprime]][count[kcurr],] = params[[kprime]][count[kcurr]-1,]
            }
    } else { # Jumping between models (steps 3 and 4)
        # 3. Set candidate parameter vector to the transformation of current
        # parameter and augmented u
        cand = Tr(kprime, params[[kcurr]][count[kcurr],], u)
        if (log(runif(1)) <= calc.post(kprime, cand) - calc.post(kcurr, params[[kcurr]][count[kcurr]+1,]) -
            ifelse(is.null(u), 0, dnorm(u, 0, 0.1, log = TRUE)) + log(kprime / kcurr)){
            params[[kprime]][count[kprime],] = cand[1:NCOL(params[[kprime]])]
            accept[count[kprime]+1,kprime] = 1
            kcurr = kprime
        } else {
            params[[kprime]][count[kprime],] = params[[kprime]][count[kprime]-1,]
            }
        }
    }

plot(1:count[1], params[[1]][1:count[1]], type = 'l')
plot(params[[2]][1:count[2],], pch = 20)

plot(1:count[1], accept[1:count[1],1], pch = 20)
plot(1:count[2], accept[1:count[2],2], pch = 20)


### There are some issues with this, I'm thinking it's mostly because of
### the choice of T_ij (the bijection). If k == 1, the model is only an
### intercept, if k == 2, the model is an intercept and slope. The
### MLEs for each model don't necessarily agree well with T_ij. By that
### I mean, theta_1 = theta - u and theta_2 = theta + u might be a bad choice.
### The problem is much more noticeable if x is generated from runif(n, 0, 15)
