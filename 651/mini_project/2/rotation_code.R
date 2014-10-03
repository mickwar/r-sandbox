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

make.R = function(theta)
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)

# unnormalize posterior
# m, s2, a, b are hyperparameters
g = function(x, m = 5, s2 = 100, a = 2.5, b = 1.5){
    require(truncnorm)
    mu = x[1]
    sig2 = x[2]
    prod(dtruncnorm(y, 1, 7, mu, sqrt(sig2)))*dnorm(mu, m, sqrt(s2)) *
        igpdf(sig2, a, b)
    }

# envelope function
#dW = function(x){
#    mu = x[1]
#    sig2 = x[2]
#    dnorm(mu, mean(y), sqrt(0.04)) * dnorm(sig2, 0.3, sqrt(0.04))
#    }
#rW = function(n)
#    matrix(c(rnorm(n,mean(y),sqrt(0.04)),rgamma(n,5.5,scale=1/16)),n,2)

modes = c(5.78, 0.31)
# about the mode of the posterior for mean (after trial and error)
cauchy.adj = modes[1]

cauchy.scale = 1/32

# this fixes the mode at modes[2], increasing alpha decrease variance
ig.alpha = 6
ig.beta = modes[2]*(ig.alpha+1)

dW = function(uv, R, off){
    xy = (uv + off) %*% t(R)
    if (xy[2] > 0){
        return (dcauchy(xy[1]-cauchy.adj, scale=cauchy.scale) *
            igpdf(xy[2], ig.alpha, ig.beta))
    } else {
        return (0)
        }
    }
rW = function(n, R, off){
    x = matrix(c(rcauchy(n, scale=cauchy.scale)+cauchy.adj,
        1/rgamma(n, ig.alpha, rate=ig.beta)),n,2)
    x %*% R - matrix(rep(off, n), n, 2, byrow=TRUE)
    }

X = rW(50000, R, offset)
plot(X[,2:1], pch=20, xlim=c(-8, 15), ylim=c(-5, 10))

mu.vec = seq(3.0, 9.0, length=100)
sig.vec = seq(0.15, 4.0, length=100)
xy = cross(mu.vec, sig.vec)

# rotations
theta = pi/4.5
R = make.R(theta)
offset = (modes %*% R - modes)
offset.mat = matrix(rep(offset, nrow(xy)), nrow(xy), 2, byrow=TRUE)
xy.rotate = xy %*% R - offset.mat

z = double(nrow(xy))
w = double(nrow(xy))
for (i in 1:length(z)){
    z[i] = log(g(xy[i,]))
    w[i] = log(dW(xy[i,], R, offset))
    }
keep = which(w != -Inf)
w = w[keep]
z = z[keep]
xy = xy[keep,]


# estimate of the constant to bring g down
(G = max(z - w))
g.star = function(x)
    log(g(x)) - G 

z.star = z - G
plot3d(cbind(xy, z.star))
points3d(cbind(xy, w), col='red')

plot3d(cbind(xy[w > -25,], w[w > -25]), col='blue')
points3d(cbind(xy[w > -25,], w[w > -25]), col='red')
points3d(cbind(xy[z.star > -25,], z.star[z.star > -25]), col='green')
#points3d(cbind(xy, h), col='blue')

plot3d(cbind(xy.rotate[w > -15,] - offset, w[w > -15]), col='red')
points3d(cbind(xy[w > -15,], w[w > -15]), col='blue')


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

niter = 100000
system.time(X <- reject(niter))
# remove the NAs and Infs
X = X[!is.na(X[,4]),]
X = X[X[,4] != Inf,]
# proportion removed
(niter - nrow(X)) / niter

# porportion of acceptances (about 0.35 so far)
mean(X[,5])

# count of acceptances
sum(X[,5])

# ratios should not exceed 1, otherwise not a correct envelope
max(X[,4])

X[X[,4] == max(X[,4]),]

at = which.max(X[,4])

# accepted draws
Y = X[X[,5] == 1, 1:2]
calc.mode(density(Y[,1]))
calc.mode(density(Y[,2]))

# joint posterior
plot(Y[sample(nrow(Y), nrow(Y)),], pch=20)

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

