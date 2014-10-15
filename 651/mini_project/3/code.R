### bootstrap importance
library(truncnorm)
dat = as.matrix(read.table("~/files/R/651/data/faculty.dat"))
colnames(dat) <- rownames(dat) <- NULL

y = as.vector(dat)

# inverse gamma density
igpdf=function(x, a, b)
	(b^a)/gamma(a)*x^(-a-1)*exp(-b/x)

g = function(x, m = 5, s2 = 100, a = 2.5, b = 1.5){
    require(truncnorm)
    mu = x[1]
    sig2 = x[2]
    if (sig2 <= 0)
        return (0)
    prod(dtruncnorm(y, 1, 7, mu, sqrt(sig2)))*dnorm(mu, m, sqrt(s2)) *
        igpdf(sig2, a, b)
    }
modes = c(5.78, 0.31)

### multivariate t
dW = function(x)
    gamma((nu + 2)/2) / (gamma(nu/2)*nu*pi*sqrt(det(sigma)) *
        (1+1/nu*t(x-modes) %*% solve(sigma) %*% (x-modes))^((nu+2)/2))
rW = function(n){
    chi = rchisq(n, nu)
    L = chol(sigma)
    z = matrix(rnorm(n * nrow(L)), n, 2)
    (z %*% L) * sqrt(nu/chi) + matrix(rep(modes, n), n ,2,
        byrow = TRUE)
    }

nu = 4
sigma = matrix(c(0.24, 0.17, 0.17, 0.23), 2, 2)
sigma = sigma/4.589

theta = rW(1000000)

g.val = apply(theta, 1, g)
I.val = apply(theta, 1, dW)

# estimate of c
(const = mean(g.val / I.val))

weights = (g.val / I.val) / const

# bootstrap sample indices
b.ind = sample(nrow(theta), nrow(theta), replace = TRUE,
    prob = weights)

# posterior draws
b.draws = theta[b.ind,]

dens = kde2d(b.draws[,1], b.draws[,2], n = 100)
dens2 = kde2d(b.draws[,1], b.draws[,2], n = 25)
dens3 = kde2d(b.draws[,1], b.draws[,2], n = 200)

filled.contour(dens$x, dens$y, log(dens$z))
filled.contour(dens2$x, dens2$y, log(dens2$z))
filled.contour(dens3$x, dens3$y, log(dens3$z))

plot(b.draws[1:10000,], pch=20)

# expectations

### metropolis-hastings
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
window = 200

f.star = function(x)
    (1 + (x - 10)^2/3)^(-2)

xx = seq(0, 20, length=100)
plot(xx, f.star(xx), type='l')

nmcmc = 10000
nburn = 5000
params = double(nmcmc+nburn)
accept = double(nmcmc+nburn)
sig = 5.90

curr.post = log(f.star(params[1]))
cand.post = curr.post

for (i in 2:(nmcmc+nburn)){
    params[i] = params[i-1]
    cand = rnorm(1, params[i-1], sig)
    cand.post = log(f.star(cand))
    if (log(runif(1)) < cand.post - curr.post){
        curr.post = cand.post
        accept[i] = 1
        params[i] = cand
        }
    if (i <= nburn && floor(i/window) == i/window)
        sig = sig * autotune(mean(accept[(i-window+1):i]),
            target = 0.25, k = window/50)
    }

params = params[(nburn+1):(nburn+nmcmc)]
accept = accept[(nburn+1):(nburn+nmcmc)]

plot(params, type='l')
mean(accept)
sig

xx = seq(0, 20, length=100)
plot(xx, f.star(xx), type='l', lwd=3)
points(density(params,bw =((4*sd(params)^5)/(3*length(params)))^(1/5)),
    type='l', col='red', lwd=3)

### problem 3 (gibbs)
dat = scan("~/files/R/651/data/ballbearing2.dat")

autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

# see p.12 in study journal X
dweib = function(x, a, b, mu = 10)
    log(b) - log(a) + (b-1)*log(x - mu) - ((x-mu)^b)/a
rweib = function(n, a, b, mu = 10)
    rweibull(n, b, a^(1/b)) + mu

# some plots
hist(dat, col='gray', freq=FALSE)
curve(exp(dweib(x,5000,2.00)), from = 10, to = 200, add = TRUE, lwd=2)
abline(v=c(10, 50, 70), lty=c(2,1,1), lwd=2)

# dweib(dat, 0.5, 0.3)
# dweibull(dat, 0.3, 0.5^(1/0.3), log = TRUE)

calc.post = function(params){
    a = params[1]
    b = params[2]
    # likelihood
    out = sum(dweib(dat, a, b))
    # priors
    # inverse gamma on a
    out = out - (a.ig + 1) * log(a) - (b.ig / a)
    # gamma? on b
    out = out + (a.ga - 1) * log(b) - (b.ga * b)
    return (out)
    }

# mean of inverse gamma to be around 5000
a.ig = 4
b.ig = 20000

b.ig/(a.ig-1)
sqrt(b.ig^2/((a.ig-1)^2*(a.ig-2)))
plot(density(1/rgamma(10000, a.ig, b.ig)), type='l')

# mean of gamma to be around 2
a.ga = 8.0
b.ga = 6.0
a.ga/b.ga
a.ga/b.ga^2
plot(density(rgamma(10000, a.ga, b.ga)), type='l')

nburn = 10000
nmcmc = 50000
params = matrix(0, nburn+nmcmc, 2)
accept = double(nburn+nmcmc)
params[1,] = c(5000, 2.0)
sig = 0.15
window = 500

curr.post = calc.post(params[1,])
cand.post = curr.post

for (i in 2:(nburn+nmcmc)){
    # sample from the inverse gamma
    params[i, 1] = 1/rgamma(1, a.ig + length(dat),
        b.ig + sum(dat ^ params[i-1, 2]))

    # sample from the gamma (using Metropolis)
    params[i, 2] = params[i-1, 2]
    cand = rnorm(1, params[i-1, 2], sig)
    if (cand > 0){
        cand.post = calc.post(c(params[i,1], cand))
        if (log(runif(1)) < cand.post - curr.post){
            curr.post = cand.post
            accept[i] = 1
            params[i, 2] = cand
            }
        }
    if (floor(i/window) == i/window && i <= nburn)
        sig = sig * autotune(mean(accept[(i-window+1):i]),
            target = 0.25, k = window/50)
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc)]

mean(accept)

# trace plots
plot(params[,1], type='l'); abline(h = 10000)
plot(params[,2], type='l'); abline(h = 2)

# posterior and prior marginals
plot(density(params[,1]), type='l', lwd=2, xlim=c(0, 1e5))
points(density(1/rgamma(100000, a.ig, b.ig)), type='l',
    col='red', lwd=2)

plot(density(params[,2]), type='l', lwd=2, xlim=c(0, 5))
points(density(rgamma(10000, a.ga, b.ga)), type='l',
    col='red', lwd=2)

apply(params, 2, mean)
apply(params, 2, quantile, c(0.025, 0.5, 0.975))

plot(params, pch=20)

z = kde2d(params[,1], params[,2], n = 150)
.filled.contour(z$x, z$y, z$z, levels=seq(min(z$z), max(z$z),
    length=15), col=rainbow(14))
points(params, pch=20, cex = 0.01)

preds = double(nmcmc)
for (i in 1:nmcmc)
    preds[i] = rweib(1, params[i,1], params[i,2])

hist(dat, col='gray', freq = FALSE)
points(density(dat), lwd=2, type='l', col='red')
points(density(preds), lwd=2, type='l')

mean(preds > 120)
