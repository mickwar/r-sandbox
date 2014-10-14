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
