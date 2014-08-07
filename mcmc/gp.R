library(MASS)


# "computer experiment data"
dat.sim = function(x)
    -2.5+5.5*x+rnorm(length(x), 0, 0.4)

true = function(x)
    -2-16*x^3+19.2*x^2
true.sim = function(x)
    -2-16*x^3+19.2*x^2+rnorm(length(x), 0, 0.5)

X = c(0, 0.1, 0.3, 0.35, 0.5, 0.6, 0.65, 0.8, 0.95, 1)
n = 30
comp.dat = matrix(dat.sim(rep(X, n)), n, length(X), byrow=T)
field.dat = true.sim(X)

xx = seq(0,1,length=100)
plot(xx, true(xx),ylim=c(-3.5,4))
for (i in 1:n)
    lines(X, comp.dat[i, ], col='red')
meanf = apply(comp.dat, 2, mean)
lines(X, meanf, col='green', lwd=4)
lines(X, field.dat, col='blue', lwd=4)


cov.func = function(x, y=x, phi=0.5){
    out = matrix(0, length(x), length(y))
    for (i in 1:length(x))
        out[i, ] = exp(-phi*sqrt((x[i]-y)^2))
    return (out)
    }


### comparison tests for computing the log determinant
### of a positive definite symmetric matrix
n = 500
xx = seq(0, 1, length=n)

A = cov.func(xx)
determinant(A, log=TRUE)$modulus[1]

B = chol(A)
2*sum(log(diag(B)))

C = eigen(A)
sum(log(C$values))

D = svd(A)
sum(log(D$d))

system.time(determinant(A, log=TRUE)$modulus[1])
system.time(B <- chol(A)) + system.time(2*sum(log(diag(B))))
system.time(C <- eigen(A)) + system.time(sum(log(C$values)))
system.time(D <- svd(A)) + system.time(sum(log(D$d)))
###
###


# GP on field data
phi = 0.5
xx = seq(0, 1, length=100)
Sigxx = cov.func(X, X, phi)
Sigyy = cov.func(xx, xx, phi)
Sigxy = cov.func(X, xx, phi)

cov.inv = solve(Sigxx)
mu.mat = t(Sigxy) %*% cov.inv %*% field.dat
sig.mat = Sigyy - t(Sigxy) %*% cov.inv %*% Sigxy + diag(length(xx))*1e-1

yfit = mvrnorm(10000, mu.mat, sig.mat)
gpmean = apply(yfit, 2, mean)
gpup95 = apply(yfit, 2, quantile, 0.975)
gplo95 = apply(yfit, 2, quantile, 0.025)

plot(X, field.dat, col='blue', lwd=4)
lines(xx, gpmean)
lines(xx, gpup95)
lines(xx, gplo95)

# GP on computer data
phi = 0.5
xx = seq(0, 1, length=100)
X = c(0, 0.1, 0.3, 0.35, 0.5, 0.6, 0.65, 0.8, 0.95, 1)
#comp.dat = as.vector(comp.dat)
Sigxx = cov.func(X, X, phi)
Sigyy = cov.func(xx, xx, phi)
Sigxy = cov.func(X, xx, phi)

cov.inv = solve(Sigxx)
mu.mat = t(Sigxy) %*% cov.inv %*% t(comp.dat)
sig.mat = Sigyy - t(Sigxy) %*% cov.inv %*% Sigxy + diag(length(xx))*1e-1

yfit = mvrnorm(1, mu.mat[,1], sig.mat)
points(xx, yfit, col='red')
