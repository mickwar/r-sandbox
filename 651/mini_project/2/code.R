library(rgl)
library(truncnorm)
dat = as.matrix(read.table("~/files/R/651/data/faculty.dat"))
colnames(dat) <- rownames(dat) <- NULL

y = as.vector(dat)

# Likelihood: Truncated Normal, unknown mean and variance
# Prior: mu ~ Normal, sigma^2 ~ I.G, indepdent of each other
# Note: the mu and sigma^2 parameters for the truncated normal
#       are not the mean variance, but do have the same
#       support as mu and sigma^2 in a regular normal

# inverse gamma density
igpdf=function(x, a, b)
    (b^a)/gamma(a)*x^(-a-1)*exp(-b/x)

cross = function(...){
    vec = list(...)
    d = length(vec)
    N = double(d+2) + 1
    for (i in 1:d)
        N[i+1] = length(vec[[i]])
    out = matrix(0, prod(N), d)
    for (i in 1:d){
        out[,i] = rep(vec[[i]], times=prod(N[1:i]),
            each=prod(N[(i+2):(d+2)]))
        }
    return(out)
    }

calc.mode = function(dens, method="loess"){
    dx = dens$x
    dy = dens$y
    if (method == "loess"){
        l = loess(dy ~ dx)
        return (l$x[which.max(l$y)])
        }
    if (method == "density")
        return (dx[which.max(dy)])
    }

# unnormalize posterior
# m, s2, a, b are hyperparameters
g = function(x, m = 5, s2 = 100, a = 2.5, b = 1.5){
    require(truncnorm)
    mu = x[1]
    sig2 = x[2]
    if (sig2 <= 0)
        return (0)
    prod(dtruncnorm(y, 1, 7, mu, sqrt(sig2)))*dnorm(mu, m, sqrt(s2)) *
        igpdf(sig2, a, b)
    }
# log of adjusted unnormalize posterior
g.star = function(x)
    log(g(x)) - G 

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

# estimating G via a grid (use for plotting)
mu.vec = seq(4.0, 10.0, length=500)
sig.vec = seq(0.05, 4.0, length=500)

# fine grid around mode (not nice for plotting)
mu.vec = seq(5.6, 6.0, length=500)
sig.vec = seq(0.25, 0.30, length=500)
xy = cross(mu.vec, sig.vec)

nu = 4
sigma = matrix(c(0.24, 0.17, 0.17, 0.23), 2, 2)
sigma = sigma/4.589
z = log(apply(xy, 1, g))
w = log(apply(xy, 1, dW))
(G = max(z - w))

G = -20.613136

xy[which.max(z-w),]

niter = 10000
system.time(X <- reject(niter))
mean(X[,5]); max(X[,4])
# remove NA
#X = X[-which(is.na(X[,5])),]
# porportion of acceptances (about 0.35 so far)
mean(X[,5])

z.star = z - G
plot3d(cbind(xy[w > -15,], (w[w > -15])), col=rgb(0.0589,0.4980,0.2812),
    xlab = "mu", ylab = "sigma^2", zlab="density")
points3d(cbind(xy[z.star > -15,], (z.star[z.star > -15])),
    col='dodgerblue')

# take a screenshot of the rgl device
#rgl.postscript("figs/env2.pdf", "pdf")

plot3d(cbind(xy, exp(w)), col='blue')
points3d(cbind(xy, exp(z.star)))


reject = function(M = 100){
    # initialize
    out = matrix(0, M, 5)

    # envelope draws
    out[,1:2] = rW(M)

    # uniform draws
    out[,3] = log(runif(M))
    
    # calculate ratios
    for (i in 1:M)
        out[i,4] = g.star(out[i,1:2]) - log(dW(out[i,1:2]))

    # accept or not
    for (i in 1:M){
        out[i,5] = ifelse(out[i,3] <= out[i,4], 1, 0)
        }

    return (out)
    }

niter = 3500000
system.time(X <- reject(niter))
# porportion of acceptances (about 0.35 so far)
mean(X[,5])

# count of acceptances
sum(X[,5])

# ratios should not exceed 0, otherwise not a correct envelope
max(X[,4])

#hist(X[,4], breaks=1000, xlim=c(0, 2), ylim=c(0, 100), col='gray')

X[X[,4] == max(X[,4]),]

at = which.max(X[,4])

# accepted draws
Y = X[X[,5] == 1, 1:2]
#calc.mode(density(Y[,1]))
#calc.mode(density(Y[,2]))

apply(Y, 2, mean)
var(Y)
apply(Y, 2, sd)
cor(Y)

preds = rtruncnorm(nrow(Y), 1, 7, Y[,1], sqrt(Y[,2]))
dens = density(preds)

pdf("figs/pred.pdf")
plot(dens, xlab = "Score", ylab = "", cex.lab = 1.5, main = "")
polygon(dens, col='gray')
lines(c(5,5),c(0,0.28))
dev.off()

mean(preds >= 5)

# joint posterior
pdf("figs/draws.pdf")
par(mar=c(4,4.9,2,1))
set.seed(18)
plot(Y[sample(nrow(Y), 10000),], pch=20, xlab = expression(mu),
    ylab = expression(sigma^2), cex.lab = 2.0)
contour(mu.vec, sig.vec, z.mat, col='dodgerblue', lwd=2, add=TRUE)
dev.off()


plot(density(Y[,2]))

# marginal posterior distributions
plot(density(Y[sample(nrow(Y), 100000),1]))
plot(density(Y[sample(nrow(Y), 100000),2]))

plot(density(y))
curve(dnorm(x, mean(Y[,1]), sqrt(mean(Y[,2]))), add=TRUE, col='blue')

### contour plot
mu.vec = seq(5.25, 6.5, length=200)
sig.vec = seq(0, 1.25, length=200)
z.mat = matrix(0, length(mu.vec), length(sig.vec))
for (i in 1:length(mu.vec))
    for (j in 1:length(sig.vec))
        z.mat[i,j] = exp(g.star(c(mu.vec[i], sig.vec[j])))
pdf("figs/contour.pdf")
par(mar=c(4,4.9,2,1))
contour(mu.vec, sig.vec, z.mat, xlim=c(5.25,6.5), ylim=c(0,1.25),
    col='dodgerblue', xlab = expression(mu), cex.main=1.5, lwd=2,
    ylab = expression(sigma^2), cex.lab = 2.0)
dev.off()

points(Y, pch=20)
contour(mu.vec, sig.vec, z.mat, add=TRUE, col='blue')

### importance sampling
theta = rW(1000000)

# takes a moment
g.val = apply(theta, 1, g)
I.val = apply(theta, 1, dW)

# estimate of c
(const = mean(g.val / I.val))

# expectations
(est.mean = apply(theta * g.val / I.val, 2, mean)/const)
(est.var = apply(theta^2 * g.val / I.val, 2, mean)/const - est.mean^2)
(est.sd = sqrt(est.var))

# bootstrap importance
weights = (g.val / I.val) / const

b.ind = sample(nrow(theta), nrow(theta), replace = TRUE,
    prob = weights)

# posterior draws
b.draws = theta[b.ind,]

plot(b.draws[1:10000,], pch=20)
