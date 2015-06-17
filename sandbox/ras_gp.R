set.seed(1)
x = c(-4.6, -4.2, -4.0, -3.8, -3.7,
      -3.1, -2.0, -1.9, -1.6, -0.5,
      -0.1,  0.1,  0.4,  0.5,  0.8,
       1.2,  1.5,  1.7,  2.5,  3.1)

y = c( 5.1,  4.8,  3.9,  3.5,  3.9,
       2.0,  1.3,  1.7,  1.9, -2.3,
      -3.5, -2.8, -2.4, -2.2, -0.4,
       1.0,  1.1,  1.3,  1.0,  2.1)

# x = x[rep(1:length(x), each=2)]
# y = y[rep(1:length(y), each=2)]
# x = x + rnorm(length(x), 0, 0.25)
# y = y + rnorm(length(y), 0, 0.25)

plot(x, y, pch=3, xlim=c(-5, 5), ylim=c(-4, 8), lwd=2, cex=2)
abline(v=(-5:5), h=seq(-4, 8, by=2), col="gray50", lty=2)

m = function(x, theta_m){
    a = theta_m[1]
    b = theta_m[2]
    c = theta_m[3]
    return (a*x^2 + b*x + c)
    }
library(fields) # for Matern
k = function(x, xprime, theta_k){
    # if xprime = x, leave xprime missing (i.e. a covariance matrix
    # is created for all possible pairs in x)
    n = NROW(x)
    
    sig_y = theta_k[1]
    sig_n = theta_k[2]
    l = theta_k[3]

    if (missing(xprime)){
        d = as.matrix(dist(x, method = "manhattan"))
        # squared exponential
        out = sig_y * exp(-(d^2) / (2*l^2)) + diag(sig_n, n)

        # Matern
#       out = sig_y * Matern(d, range=l, nu=1, phi=1) + diag(sig_n, n)
    } else {
        m = NROW(xprime)
        d = as.matrix(dist(c(x, xprime), method = "manhattan"))
        d = d[1:n, (n+1):(n+m)]
        # squared exponential
        out = sig_y * exp(-(d^2) / (2*l^2)) 

        # Matern
#       out = sig_y * Matern(d, range=l, nu=1, phi=1)
        }
    return (out)
    }

# log posterior
calc.post = function(params){
    mu = m(x, theta_m = params[1:3])
    sigma = k(x, theta_k = params[4:6])
    # likelihood
    out = -0.5 * determinant(sigma)$modulus[1] +
        -0.5 * t(y-mu) %*% (solve(sigma) %*% (y-mu))
    # priors
    out = out + dnorm(params[1], 0, 1, log = TRUE)
    out = out + dnorm(params[2], 0, 1, log = TRUE)
    out = out + dnorm(params[3], 0, 1, log = TRUE)
    out = out + dgamma(params[4], 1, 1, log = TRUE)
    out = out + dgamma(params[5], 1, 1, log = TRUE)
    out = out + dgamma(params[6], 1, 1, log = TRUE)
    return (out)
    }

# mcmc setting
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

nburn = 10000
nmcmc = 10000
nparams = 6
params = matrix(0, nburn+nmcmc, nparams)
# starting values
params[1,] = 1
sigs = rep(1, nparams)

# initialize log posterior value
post = calc.post(params[1,])
cand.param = params[1,]

#
lower = c(-Inf, -Inf, -Inf, 0, 0, 0)
upper = c(Inf, Inf, Inf, Inf, Inf, Inf)

# confidence intervals on acceptance rates?
accept = matrix(0, nburn+nmcmc, nparams)
window = 100

# mcmc loop
for (i in 2:(nburn+nmcmc)){
    cat("\rIteration",i,"/",nburn+nmcmc)
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
    if (i == (nburn+nmcmc))
        cat("\n")
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

par(mfrow=c(3,2))
plot(params[,1], type='l')
plot(params[,2], type='l')
plot(params[,3], type='l')
plot(params[,4], type='l')
plot(params[,5], type='l')
plot(params[,6], type='l')

par(mfrow=c(1,1))

apply(params, 2, mean)
apply(accept, 2, mean)

sigs

library(MASS)
pred.x = seq(-6, 6, length = 100)
pred.y = matrix(0, nmcmc, length(pred.x))
for (i in 1:nmcmc){
    cat("\rIteration",i,"/",nmcmc)
    t_m = params[i,1:3]
    t_k = params[i,4:6]

    kstar = k(x, pred.x, t_k)
    sigma_inv = solve(k(x, theta_k = t_k))
    sigstar = k(pred.x, theta_k = t_k)

    post.mean = m(pred.x, t_m) + t(kstar) %*% sigma_inv %*%
        (y - m(x, t_m))
    post.var = sigstar - t(kstar) %*% sigma_inv %*% kstar

    pred.y[i,] = mvrnorm(1, post.mean, post.var)
    if (i == nmcmc)
        cat("\n")
    }

pred.upper = apply(pred.y, 2, quantile, 0.975)
pred.mean = apply(pred.y, 2, mean)
pred.lower = apply(pred.y, 2, quantile, 0.025)

# plot(x, y, pch=3, xlim=c(-5, 5), ylim=c(-4, 8), cex=2, lwd=2)
# abline(v=(-5:5), h=seq(-4, 8, by=2), col="gray50", lty=2)
# polygon(c(pred.x, rev(pred.x)), c(pred.upper, rev(pred.lower)),
#     col='gray85')
# lines(pred.x, pred.mean, col='red', lwd=2)
# lines(pred.x, pred.y[1,], lty=2)
# lines(pred.x, pred.y[2,], lty=2)
# lines(pred.x, pred.y[3,], lty=2)
# points(x, y, pch=3, xlim=c(-5, 5), ylim=c(-4, 8), cex=2, lwd=2)
# legend(-2.5, 7.5, legend = c("Data", "Posterior predictive mean",
#     "Single predictive draws", "95% predictive region"),
#     lty=c(NA,1,2,1), pch=c(3, NA, NA, NA), col=c(1,2,1,"gray85"),
#     lwd=c(2,2,1,8), bg="white", pt.cex=c(2,1,1,1))
# box()


### maximum likelihood estimates (much faster)
# produces similar predictions (smaller variance though)
g = function(p)
    0.5*determinant(k(x, theta_k=p[4:6]))$modulus[1] +
        0.5*t(y-m(x, p[1:3])) %*% (solve(k(x, theta_k=p[4:6])) %*%
        (y-m(x, p[1:3])))

(mle.p = optim(rep(1, 6), g, lower = c(rep(-Inf, 3), rep(0, 3)),
    method = "L-BFGS-B")$par)

MU = m(x, mle.p[1:3])
SIG = k(x, theta_k = mle.p[4:6])


t_m = mle.p[1:3]
t_k = mle.p[4:6]

kstar = k(x, pred.x, t_k)
sigma_inv = solve(k(x, theta_k = t_k))
sigstar = k(pred.x, theta_k = t_k)

post.mean = m(pred.x, t_m) + t(kstar) %*% sigma_inv %*%
    (y - m(x, t_m))
post.var = sigstar - t(kstar) %*% sigma_inv %*% kstar

mle.y = mvrnorm(10000, post.mean, post.var)

mle.upper = apply(mle.y, 2, quantile, 0.975)
mle.mean = apply(mle.y, 2, mean)
mle.lower = apply(mle.y, 2, quantile, 0.025)

plot(x, y, pch=3, xlim=c(-5, 5), ylim=c(-4, 8), cex=2, lwd=2)
abline(v=(-11:11), h=seq(-4, 8, by=2), col="gray50", lty=2)
polygon(c(pred.x, rev(pred.x)), c(pred.upper, rev(pred.lower)),
    col='lightgreen')
polygon(c(pred.x, rev(pred.x)), c(mle.upper, rev(mle.lower)),
    col='pink')
lines(pred.x, pred.mean, col='green', lwd=2)
lines(pred.x, mle.mean, col='red', lwd=2)
# lines(pred.x, mle.y[1,], lty=2)
# lines(pred.x, mle.y[2,], lty=2)
# lines(pred.x, mle.y[3,], lty=2)
points(x, y, pch=3, xlim=c(-5, 5), ylim=c(-4, 8), cex=2, lwd=2)
legend(-2.5, 7.5, legend = c("Data", "Posterior predictive mean",
    "MLE predictive mean", "95% posterior predictive region",
    "95% MLE predictive region"),
    lty=c(NA,1,1,1,1), pch=c(3, NA, NA, NA, NA), col=c("black", "green", "red",
        "lightgreen", "pink"),
    lwd=c(2,2,2,8,8), bg="white", pt.cex=c(2,1,1,1,1))
box()


