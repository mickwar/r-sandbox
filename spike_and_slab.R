# Remarks on logistic spike and slab (mixture of two normals on the betas)
#
# The estimates for the betas tended to be shrunk toward zero. Say beta_1 = 8,
# then this would usually be estimated at something lower, say around 4.
# This is possibly because of a small sample size, increasing n might correct this
# If a large beta was included, say beta_2 = 60, then the new estimate for
# beta_1 would be even less than 4, and beta_2 itself would be much smaller than
# the truth. Is this because of all the superfluous parameters that when taking
# a linear combination still yields values away from 0 and the logistic function
# is sensitive to even moderate values? exp(3) is nearly 3 times as large as exp(2).
#
# The model did seem to do well in selecting the variables. Without any very large values,
# the true positive and true negative rates were close to 1. In the presense of the
# large values, we see false positives at about the same rate, but false negatives
# increase as those non-zero coefficients that are close to zero (in comparison to
# the large betas) are less likely to be selected.
#
# The model also seems pretty sensitive to the choice of g_j and tau_j^2 which specify
# the prior for the betas.
#
# Since the bernoulli likelihood produces a non-conjugate posterior for the betas,
# this means we require Metropolis-Hasting to obtain posterior samples for beta.
# Traceplots for beta showed very high autocorrelation. Some posterior correlation
# was also present (around 0.4-0.6 for the non-zero betas).

# Logistic spike and slab
spike_slab = function(dat, nmcmc, nburn){
    require(MASS)
    require(mwBASE)
    y = dat$y
#   X = cbind(1, dat$X)
    X = as.matrix(dat$X)
    n = NROW(X)
    p = NCOL(X)
    gj = dat$g.vec      # g_j
    tauj = dat$tau.vec  # tau_j^2

    XtX = t(X) %*% X
    Xty = t(X) %*% y

    # what are good values to choose? similarly for gj and tauj?
    # hyperpriors
#   a.w = rep(1, p+1)
#   b.w = rep(1, p+1)
    a.w = rep(1, p)
    b.w = rep(1, p)
    a.sig = 1
    b.sig = 1
#   a.w = dat$a.w
#   b.w = dat$b.w
#   a.sig = dat$a.sig
#   b.sig = dat$b.sig

    # +1 is for intercept
#   param.beta = matrix(0, nburn + nmcmc, p+1)
#   param.gamma = matrix(0, nburn + nmcmc, p+1)
#   param.w = matrix(0, nburn + nmcmc, p+1)
    param.beta = matrix(0, nburn + nmcmc, p)
    param.gamma = matrix(0, nburn + nmcmc, p)
    param.w = matrix(0, nburn + nmcmc, p)
#   param.sig2 = double(nburn + nmcmc)

    # logistic regression
    cand_sig = diag(1/p, p)
    cand_chol = chol(cand_sig)
    accept = double(nburn + nmcmc)
    calc.beta = function(beta, gamma, y, X){
        prob = 1/(1+exp(-X %*% beta))
        out = sum(log(prob[y == 1])) + sum(log(1-prob[y == 0]))
        D.vec = (gamma*(gj-1)+1) * tauj
        out = out -sum(0.5*log(D.vec)) -0.5*sum(D.vec * beta^2)
        return (out)
        }

    # starting values
    param.w[1,] = 0.5
#   param.sig2[1] = 1
    
    # sample (all are known closed-form posteriors)
    window = 100
    for (i in 2:(nburn + nmcmc)){
        if (floor(i / window) == i/window)
            cat(i, "/", nburn + nmcmc, "\r")
        # update beta
#       Dinv = diag(1/((param.gamma[i-1,]*(gj-1)+1) * tauj))
#       sigma = solve(XtX / param.sig2[i-1] + Dinv)
#       mu = sigma %*% Xty / param.sig2[i-1]
#       param.beta[i,] = mvrnorm(1, mu, sigma)
        param.beta[i,] = param.beta[i-1,]
        cand = t(param.beta[i-1,] + rnorm(p) %*% cand_chol)
        if (log(runif(1)) <= calc.beta(cand, param.gamma[i-1,], y, X) -
            calc.beta(param.beta[i-1,], param.gamma[i-1,], y, X)){
            param.beta[i,] = cand
            accept[i] = 1
            }
        if ((floor(i/window) == i/window) && (i <= nburn)){
            cand_sig = autotune(mean(accept[(i-window+1):i]), target = 0.25, k = window / 50) *
                (cand_sig + var(param.beta[(i-window+1):i,])/i)
            cand_chol = chol(cand_sig)
            }

        # update gamma
        h1 = param.w[i-1,] * dnorm(param.beta[i,], 0, gj*sqrt(tauj))
        h2 = (1-param.w[i-1,]) * dnorm(param.beta[i,], 0, sqrt(tauj))
#       param.gamma[i,] = rbinom(p+1, 1, h1 / (h1 + h2))
        param.gamma[i,] = rbinom(p, 1, h1 / (h1 + h2))

        # update w
#       param.w[i,] = rbeta(p+1, a.w + param.gamma[i,], b.w + 1 - param.gamma[i,]) 
        param.w[i,] = rbeta(p, a.w + param.gamma[i,], b.w + 1 - param.gamma[i,]) 

        # update sigma^2
#       param.sig2[i] = 1/rgamma(1, a.sig + n/2, b.sig + 0.5*sum((y - X %*% param.beta[i,])^2))

        }
    cat("\n", mean(accept), "\n")

    # removed burn-in
    param.beta = tail(param.beta, nmcmc)
    param.gamma = tail(param.gamma, nmcmc)
    param.w = tail(param.w, nmcmc)
#   param.sig2 = tail(param.sig2, nmcmc)

    # output
#   out = list("beta" = param.beta, "gamma" = param.gamma,
#       "w" = param.w, "sig2" = param.sig2)
    out = list("beta" = param.beta, "gamma" = param.gamma,
        "w" = param.w)
    return (out)

    }


library(MASS)
make.sigma = function(p, rho){
    out = matrix(0, p, p)
    inds = expand.grid(1:p, 1:p)
    fill = rho^apply(inds, 1, function(x) abs(diff(x)))
    out[1:p^2] = fill
    return (out)
    }
logit = function(p)
    log(p / (1-p))
ilogit = function(x)
    exp(x) / (1 + exp(x))

# set factors
n = 5000
p = 20
rho = 0.0
sig2 = 10

sigma = make.sigma(p, rho) * sig2
beta = c(rep(8, 5), rep(-4, 5),
    rep(2.5, 5), rep(0, p-15))
# if (i1 == 1){
#     beta = c(rep(3, 5), rep(0, p-5))   
# } else {
#     }
# beta[25] = 60
# beta[26] = -40


# simulate data
set.seed(1)
X = mvrnorm(n, rep(0, p), sigma)

X = cbind(1, matrix(rnorm(n*2), n, 2))
beta = c(sig2, 2, 0)
probs = ilogit(X %*% beta)
y = rbinom(n, size = 1, prob = probs)


dat = list("y" = y, "X" = X, "g.vec" = rep(10, p+1), "tau.vec" = rep(0.1, p+1))

mcmc = spike_slab(dat, nmcmc = 5000, nburn = 10000)

colMeans(mcmc$beta)
colMeans(mcmc$gamma)
colMeans(mcmc$w)
# mean(mcmc$sig2)
colMeans(mcmc$w) > 0.5
ind = which(colMeans(mcmc$w) > 0.5)

plot(mcmc$beta[,1], type = 'l')
plot(mcmc$w[,1], type = 'l')
# plot(mcmc$sig2, type = 'l')

mod = glm(y ~ 0+X, family = binomial(link = "logit"))
# 
summary(mod)
preds = predict(mod)
plot(probs, ilogit(preds))
abline(0, 1)

plot(X[,2], y)
points(X[,2], probs, col = 'red')
points(X[,2], ilogit(preds), col = 'blue')
points(X[,1], (ilogit(X %*% beta) >= 0.5) * 1, col = 'blue')

points(X[,1], y, col = 'green')

length(y[which(probs <= 0.10)]) / n
mean(y[which(probs <= 0.20)]
mean(y[which(probs <= 0.50)])

cbind(y, probs)
