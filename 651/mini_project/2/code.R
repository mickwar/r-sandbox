library(rgl)
library(truncnorm)
dat = as.matrix(read.table("./faculty.dat"))
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
    log(g(x)) -G 


modes = c(5.78, 0.31)

### t and ig prior
# # about the mode of the posterior for mean (after trial and error)
# t.adj = 16
# t.move = 4.3
# t.df = 8
# t.ncp = 25
# 
# # this fixes the mode at modes[2], increasing alpha decrease variance
# ig.alpha = 8
# ig.beta = modes[2]*(ig.alpha+1)
# 
# # envelope function
# dW = function(x){
#     mu = x[1]
#     sig2 = x[2]
#     dt((mu-t.move)*t.adj, t.df, t.ncp) * igpdf(sig2, ig.alpha, ig.beta)
#     }
# rW = function(n)
#     x = matrix(c(rt(n, t.df, t.ncp), 1/rgamma(n, ig.alpha, rate=ig.beta)),n,2)

### envelope -- multivariate normal
dW = function(x)
    ((2*pi)^2*det(sigma))^(-1/2) * exp(-0.5* t(x-modes) %*%
        solve(sigma) %*% (x-modes))
rW = function(n){
    L = chol(sigma)
    z = matrix(rnorm(n * nrow(L)), n, 2)
    z %*% L + matrix(rep(modes, n), n ,2, byrow = TRUE)
    }

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

x = rW(5000)
plot(x, pch=20)

# estimating G via a grid
mu.vec = seq(4.0, 20.0, length=100)
sig.vec = seq(0.15, 20.0, length=100)
xy = cross(mu.vec, sig.vec)

nu = 4
sigma = matrix(c(0.24, 0.17, 0.17, 0.23), 2, 2)
sigma = sigma/3.5
z = log(apply(xy, 1, g))
w = log(apply(xy, 1, dW))
(G = max(z - w))

z.star = z - G
plot3d(cbind(xy[w > -25,], (w[w > -25])), col='blue')
points3d(cbind(xy[z.star > -25,], (z.star[z.star > -25])))

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
calc.mode(density(Y[,1]))
calc.mode(density(Y[,2]))

# joint posterior
plot(Y[sample(nrow(Y), 100000),], pch=20)

plot(density(Y[,2]))

# marginal posterior distributions
plot(density(Y[sample(nrow(Y), 100000),1]))
plot(density(Y[sample(nrow(Y), 100000),2]))

plot(density(y))
curve(dnorm(x, mean(Y[,1]), sqrt(mean(Y[,2]))), add=TRUE, col='blue')

mu.vec = seq(3.0, 9.0, length=100)
sig.vec = seq(0.1, 4.0, length=100)
z.mat = matrix(0, length(mu.vec), length(sig.vec))
for (i in 1:length(mu.vec))
    for (j in 1:length(sig.vec))
        z.mat[i,j] = exp(g.star(c(mu.vec[i], sig.vec[j])))


contour(mu.vec, sig.vec, z.mat, xlim=c(5,7), ylim=c(0,1), col='blue')
points(Y, pch=20)
contour(mu.vec, sig.vec, z.mat, add=TRUE, col='blue')

