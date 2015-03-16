### bootstrap importance
library(truncnorm)
y = scan("~/files/data/651/faculty.dat")

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

# subset to be within decent range
keep = which(b.draws[,1] <= quantile(b.draws[,1], 0.997) &
    b.draws[,2] <= quantile(b.draws[,2], 0.997))
b.sub = b.draws[keep,]

plot(b.draws, pch=20)

library(MASS)
library(fields)
system.time(dens <- kde2d(b.sub[,1], b.sub[,2], n = 100))

# default mar = c(5.1, 4.1, 4.1, 2.1)
pdf("./figs/boot_joint.pdf", width = 7, height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))
filled.contour(dens, color.palette = tim.colors, xlab = expression(mu),
    ylab = expression(sigma^2), cex.lab = 2, cex.main = 2,
    main = expression(paste("Joint Posterior for (", mu,", ", sigma^2,")")))
dev.off()

#filled.contour(dens$x, dens$y, dens$z)
#
#plot(b.sub, type='n')
#.filled.contour(dens$x, dens$y, dens$z,
#    levels = pretty(range(dens$z, finite = TRUE), 15),
#    col = tim.colors(14))
#points(b.sub[1:10000,], pch=20, cex = 0.15)

#plot(b.draws, pch=20, cex = 1.00)

# expectations
apply(b.draws, 2, mean)
var(b.draws)
sqrt(diag(var(b.draws)))

# posterior predictive
preds = rtruncnorm(1000000, 1, 7, b.draws[,1], sqrt(b.draws[,2]))

pdf("./figs/boot_pred.pdf")
plot(density(y), col='blue', lwd = 2, xlab = "Teacher Rating", ylab = "",
    main = "Posterior Predictive", cex.lab = 2, cex.main = 2, lty = 2)
points(density(preds), lwd=4, type = 'l')
legend("topleft", legend = c("Data", "Predictive"), lwd = c(2, 4),
    col = c("blue", "black"), lty = c(2, 1), box.lty = 0, cex = 1.5)
lines(c(5, 5), c(0, 0.28), lwd=3, col = 'red')
dev.off()

mean(preds >= 5)

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

pdf("./figs/mh_trace.pdf")
plot(params, type='l', ylab="x", cex.lab = 1.5, main = "Trace Plot", cex.main = 1.5)
dev.off()

mean(accept)
sig

xx = seq(0, 20, length=500)
dens = density(params,bw =((4*sd(params)^5)/(3*length(params)))^(1/5))

pdf("./figs/mh_post.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(xx, f.star(xx), type='l', lwd=2, xlab = "x", ylab = "density",
    col = "blue", lty = 2)
points(dens, type='l', lwd=4)
legend("topleft", box.lty = 0, legend = c("f*(x)", "Posterior"),
    col = c("blue", "black"), lwd = c(2, 4), cex = 1.5, lty = c(2, 1))

plot(xx, f.star(xx) * 0.35, type='l', lwd=2, xlab = "x", ylab = "density",
    col = "blue", lty = 2)
points(dens, type='l', col='black', lwd=4)
legend("topleft", box.lty = 0, legend = c("f*(x) * 0.35", "Posterior"),
    col = c("blue", "black"), lwd = c(2, 4), cex = 1.5, lty = c(2, 1))
dev.off()


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

# HPD
hpd.uni = function(x, prob = 0.95, precision = 1000){
    range = seq(0, 1-prob, length=precision)
    range = cbind(range, range+prob)
    best = range[which.min(apply(range, 1, function(y)
        diff(quantile(x, y)))),]
    return (quantile(x, best))
    }


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

# some plots (get an idea for the prior specification)
hist(dat, col='gray', freq=FALSE)
curve(exp(dweib(x, 800, 1.50)), from = 10, to = 200, add = TRUE, lwd=2)
abline(v=c(10, 50, 70), lty=c(2,1,1), lwd=2)

# mean of inverse gamma to be around 800
a.ig = 4
b.ig = 2500

b.ig/(a.ig-1)
sqrt(b.ig^2/((a.ig-1)^2*(a.ig-2)))
plot(density(1/rgamma(10000, a.ig, b.ig)), type='l')

# mean of gamma to be around 1.5 (+/- 1)
a.ga = 3.0
b.ga = 2.0
a.ga/b.ga
a.ga/b.ga^2
plot(density(rgamma(10000, a.ga, b.ga)), type='l')

nburn = 20000
nmcmc = 100000
params = matrix(0, nburn+nmcmc, 2)
accept = double(nburn+nmcmc)
params[1,] = c(800, 1.5)
sig = 0.22
window = 1000

for (i in 2:(nburn+nmcmc)){
    # sample from the inverse gamma
    params[i, 1] = 1/rgamma(1, a.ig + length(dat),
        b.ig + sum((dat-10) ^ params[i-1, 2]))

    # sample from the gamma (using Metropolis)
    params[i, 2] = params[i-1, 2]
    cand = rnorm(1, params[i-1, 2], sig)
    if (cand > 0){
        curr.post = calc.post(c(params[i, 1], params[i-1, 2]))
        cand.post = calc.post(c(params[i, 1], cand))
        if (log(runif(1)) < cand.post - curr.post){
            accept[i] = 1
            params[i, 2] = cand
            }
        }
    if (floor(i/window) == i/window && i <= nburn)
        sig = sig * autotune(mean(accept[(i-window+1):i]),
            target = 0.25, k = 2)
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc)]

mean(accept)

# trace plots
plot(params[,1], type='l'); abline(h = 800)
plot(params[,2], type='l'); abline(h = 1.6)

# posterior and prior marginals
pdf("./figs/gibb_marginal.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))
plot(density(1/rgamma(100000, a.ig, b.ig)), type='l', cex.lab = 2, ylab = "",
    col='red', lwd=2, xlim = c(0, 5000), xlab = expression(alpha), main = "")
points(density(params[,1]), type='l', lwd=4)
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "black"),
    lwd = c(2, 4), box.lty = 0, cex = 1.5)

plot(density(params[,2]), type='l', lwd=4, xlim=c(0, 5), cex.lab = 2,
    ylab = "", xlab = expression(beta), main = "")
points(density(rgamma(10000, a.ga, b.ga)), type='l',
    col='red', lwd=2)
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "black"),
    lwd = c(2, 4), box.lty = 0, cex = 1.5)
title("Marginals", outer = TRUE, line = -2, cex.main = 1.5)
dev.off()

apply(params, 2, mean)
var(params)
apply(params, 2, quantile, c(0.025, 0.5, 0.975))

plot(params, pch=20)

library(MASS)
library(fields)
system.time(z <- kde2d(params[,1], params[,2], n = 250))

pdf("./figs/gibb_joint.pdf", width = 7, height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))
filled.contour(z, xlim=c(0, 6000), color.palette = tim.colors, cex.main = 1.5,
    xlab = expression(alpha), ylab = expression(beta), cex.lab = 1.5,
    main = expression(paste("Joint Posterior for (",alpha,", ",beta,")")))
dev.off()

#plot(params, type='n')
#.filled.contour(z$x, z$y, z$z, levels=seq(min(z$z), max(z$z),
#    length=15), col=tim.colors(14))
#points(params, pch=20, cex = 0.01)

# posterior predictive
preds = double(nmcmc)
for (i in 1:nmcmc)
    preds[i] = rweib(1, params[i,1], params[i,2])

pdf("./figs/gibb_pred.pdf")
hist(dat, col='gray', freq = FALSE, main = "Posterior Predictive", ylab = "",
    xlab = "Revolutions before failure", cex.lab = 2, cex.main = 2)
points(density(dat), lwd=2, type='l', col='blue', lty=2)
points(density(preds), lwd=4, type='l')
legend("topright", legend = c("Data", "Density estimate of data", "Predictive"),
    col = c("gray", "blue", "black"), lwd = c(6, 2, 4), lty = c(1, 2, 1),
    box.lty = 0, cex = 1.5)
dev.off()

# 97% hpd
(hpd.alpha = hpd.uni(params[,1], prob = 0.97))
(hpd.beta = hpd.uni(params[,2], prob = 0.97))

# joint 97% hpd
plot(params, pch=20, cex=0.5)
# evaluate the joint posterior density for each draws
z = apply(params, 1, calc.post)
h.ind = which(z > quantile(z, 1.00 - 0.97))
p2 = params[h.ind,]
points(p2, pch=20, cex=0.5, col='yellow')
# marginals (notice that using the joint results in a bigger interval)
c(hpd.alpha, diff(hpd.alpha))
c(range(p2[,1]), diff(range(p2[,1])))
c(hpd.beta, diff(hpd.beta))
c(range(p2[,2]), diff(range(p2[,2])))
# how to get marginal hpds using this density method? just use density function?
# use the full conditional for each parameter? i would think either should
# give about the same answer, but using the density function is less accurate,
# but the full conditional might be harder to get?


# probability of a superior ball bearing
mean(preds > 120)

# final settings
nburn
nmcmc
sig
mean(accept)
