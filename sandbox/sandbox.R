f1 = function(x)
    1/384*(1/6*x^3-2*x+8/3)
f2 = function(x)
    1/384*(8*x-64/3)
f3 = function(x)
    1/384*(-1/6*x^3+58*x-1064/3)

xx1 = seq(2,6,length=100)
xx2 = seq(6,10,length=100)
xx3 = seq(10,14,length=100)

int = function(FUN, range, blocks){
    out = 0
    x = seq(range[1], range[2], length=(blocks+1))
    size = x[2]-x[1]
    for (i in 1:blocks)
        out = out + FUN(x[i])*size   # Left-hand integration
    return (out)
    }

plot(xx1, f1(xx1),xlim=c(2,14),ylim=c(0,0.75))
points(xx2, f2(xx2),col='blue')
points(xx3, f3(xx3),col='red')

int(f1, c(2,6), 10000)+
int(f2, c(6,10), 10000)+
int(f3, c(10,14), 10000)


# basic mcmc example
p = function(x){
    out = double(length(x))
    for (i in 1:length(x)){
        if ( x[i] == 0) out[i] = (1/6)
        if ( x[i] == 1) out[i] = (2/6)
        if ( x[i] == 2) out[i] = (3/6)
        }
    return (out)
    }

run = function(A, x, nreps){
    for (i in 1:nreps){
        cat(i, ": ", paste(round(x,5),collapse=" "),"\n")
        x = x %*% A
        }
    }

A = matrix(0, 3, 3)
A[1, ] = c(1/2, 1/2, 0)
A[2, ] = c(1/4, 1/4, 1/2)
A[3, ] = c(0, 1/3, 2/3)

# yA^k -> (p(0), p(1), p(2)) as k->Inf
run(A, c(1,0,0), 15)
p(c(0, 1, 2))

mcmc = function(x, nreps){
    accept = 0
    out = double(nreps)
    out[1] = x
    for (i in 2:nreps){
        out[i] = out[i-1]
        cand.x = out[i] + sample(c(-1,1), 1)
        if (cand.x >= 0 && cand.x <= 2){
            ratio = min(1, p(cand.x)/p(out[i-1]))
            prob.check = runif(1)
            if (prob.check <= ratio){
                out[i] = cand.x
                accept = accept + 1
                }
            }
        }
    return (list(accept/nreps, out))
    }

out = mcmc(0, 100000)
ap = out[[1]]
draws = out[[2]]
mean(draws==0)
mean(draws==1)
mean(draws==2)

# lag and auto-correlation
x = draws
ar0 = cor(x, x) # lag 0

x = draws[-length(draws)]
y = draws[-1]
ar1 = cor(x, y) # lag 1

x = draws[-c(length(draws)-1, length(draws))]
z = draws[-c(1,2)]
ar2 = cor(x, z) # lag 2

# effective sample size
a = acf(draws)
factor = 2*sum(a$acf)-1
ess = length(draws)/factor


# max of three m exponentials
n = 100000
m = 14
lambda = 1/17
draws = matrix(rexp(n*m, lambda), n, m)
draws = cbind(draws, 0)
for (i in 1:n)
    draws[i, 4] = max(draws[i, ])
f = function(x, lambda)
    (m/lambda*(1-exp(-x/lambda))^(m-1))*(exp(-x/lambda))

plot(density(draws[, 1]), col=1)
for (i in 2:4)
    lines(density(draws[, i]), col=i)
xx = seq(0,max(draws[,4]),length=100)
points(xx, f(xx, 1/lambda))



# choleski and random multivariate normal
library(MASS)
mu = c(1,2)
A = matrix(c(1,0.9,0.9,1),2,2)
B = chol(A)

n = 10000
draws = matrix(0, n, 2)
for (i in 1:n)
    draws[i, ] = mu + t(B) %*% rnorm(2)

draws2 = mvrnorm(n, mu, A)
plot(draws)
points(draws2, col='red')


# integration
# riemann: sum [ f(x_i) * prod(a_i, b_i) ]
# value of function times the "area", all summed up
dmvnorm = function(X, mu, Sigma)
    1/(2*pi)^(length(X)/2) * det(Sigma)^(-1/2) * exp(-1/2*t(X - mu) %*% solve(Sigma) %*% (X - mu))
dmvnorm = function(X)
    1/(2*pi)^(length(X)/2) * exp(-1/2*t(X) %*% (X))
cross = function(...){ # also see expand.grid()
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

# one-dimensional
n = 10000
d = 1
xx = seq(0, 5, length=n)
yx1 = matrix(0, n, 2)
yx1[, 2] = xx
yx1[, 1] = apply(as.matrix(yx1[, 2]), 1, dmvnorm)
sum(yx1[,1] * (xx[2] - xx[1])) # 0.5 as n -> Inf, and max(x) -> Inf

# two-dimensions
n = 100
d = 2
xx = seq(0, 5, length=n)
yx2 = matrix(0, n^d, d+1)
yx2[, 2:(d+1)] = cross(xx, xx)
yx2[, 1] = apply(yx2[, 2:(d+1)], 1, dmvnorm)
sum(yx2[,1] * (xx[2]-xx[1])^d) # 0.25

# three-dimensions
n = 25
d = 3
xx = seq(0, 4, length=n)
yx3 = matrix(0, n^d, d+1)
yx3[, 2:(d+1)] = cross(xx, xx, xx)
yx3[, 1] = apply(yx3[, 2:(d+1)], 1, dmvnorm)
sum(yx3[,1] * (xx[2]-xx[1])^d) # 0.125
# plot
mu = c(1,2,3)
sig = matrix(c(1, 0.9, 0.6, 0.9, 1, 0.7, 0.6, 0.7, 1), 3, 3)
ddtest = function(x)
    dmvnorm(x, mu, sig)
n = 50
d = 3
aa = seq(mu[1]-3, mu[1]+3, length=n)
bb = seq(mu[2]-3, mu[2]+3, length=n)
cc = seq(mu[3]-3, mu[3]+3, length=n)
hh = matrix(0, n^d, d+1)
hh[, 1:d] = cross(aa, bb, cc)
hh[, 4] = apply(hh[,1:d], 1, ddtest)
cols = (hh[,4] - min(hh[,4]))/(max(hh[,4])-min(hh[,4]))
plot3d(hh[,1:d], col=gray(1-cols), type='s', size=cols) #slow

# four-dimensions
n = 25
d = 4
xx = seq(0, 4, length=n)
yx4 = matrix(0, n^d, d+1)
yx4[, 2:(d+1)] = cross(xx, xx, xx, xx)
yx4[, 1] = apply(yx4[, 2:(d+1)], 1, dmvnorm)
sum(yx4[,1] * (xx[2]-xx[1])^d) # 0.0625

# try composite trapezoidal rule


# middle-square random numbers
msr = function(n, seed=123456){
    # count digits
    digit = function(x){
        d = 1
        while (x / 10^d >= 1)
            d = d + 1
        return (d)
        }
    m = digit(seed)
    out = numeric(n)
    out[1] = seed
    for (i in 2:(n+1)){
        new = out[i-1]^2 
        m2 = digit(new)
        start = ceiling(m/2)
        end = floor(m/2) + ((m2+1) %% 2)
        if (!(m %% 2) && !(m2 %% 2)){    # both even
            start = start + 1
            end = end - 1
            }
        out[i] = floor((new %% 10^(m2 - start + 1)) / 10^end)
        }
    return (out[-1])
    }

msr(100, seed=9187988)


# t-distr
tdist = function(x, df)
    (gamma((df+1)/2) / gamma(df/2)) * (1 / sqrt(df * pi)) * (1 / (1+(x^2)/df)^((df+1)/2))

xx = seq(-5, 5, length=100)
plot(xx, tdist(xx, 0.5), type='l', col='red', ylim=c(0,0.35))
lines(xx, tdist(xx, 1), col='blue')
lines(xx, tdist(xx, 2), col='green')
lines(xx, tdist(xx, 0.25))
lines(xx, tdist(xx, 1/8))
lines(xx, tdist(xx, 1/16))

d1 = rcauchy(100000)
d2 = rt(100000, 1/4)
d3 = rt(100000, 1/16)
d4 = rt(100000, 1/64)

mean(d1)
mean(d2)
mean(d3)
var(d1)
var(d2)
var(d3)
var(d4)

ndraws = 10000
mu = matrix(c(5,1), 2, 1)
cov = diag(2)*3
cov[1,2]=0.5
cov[2,1]=0.5
A=svd(cov)
U=A$u
S=A$d


draw = (U%*%sqrt(S))%*%(rnorm(ndraws))
draw = t(draw)

plot(draw)


# sum of normals to make normal(0, 1)
n = 1000000
means = c(-2,-3,0,9,2) # on test it was c(-2,-1,0,1,2)
sds = c(9,2,13,4,5)    # and c(1,2,3,4,5)
draws1 = rnorm(n, means[1], sds[1])
draws2 = rnorm(n, means[2], sds[2])
draws3 = rnorm(n, means[3], sds[3])
draws4 = rnorm(n, means[4], sds[4])
draws5 = rnorm(n, means[5], sds[5])

out = (draws1 - means[1])/(sqrt(5)*sds[1]) +
      (draws2 - means[2])/(sqrt(5)*sds[2]) + 
      (draws3 - means[3])/(sqrt(5)*sds[3]) + 
      (draws4 - means[4])/(sqrt(5)*sds[4]) + 
      (draws5 - means[5])/(sqrt(5)*sds[5])

mean(out)
sd(out)

plot(density(out))
curve(dnorm(x,0,1), add=T, from=-4, to=4, col='red')


# qnorm and pnorm discrepancy
x = seq(6.0, 8.292, by=0.001) # no problems less than around 7
y = qnorm(pnorm(x))
z = x - y
plot(z)
abline(0, 0, col='red')
plot(z[x>7.8])
abline(0, 0, col='red')


##########
# qq-plot
y = rnorm(100)
p = (1:length(y) - 0.5)/length(y) # continuity correction
qqnorm(y, xlim=c(-3.5,3.5), ylim=c(-3.5,3.5))
pars = par("usr")
abline(0,1)
sy = sort(y)
qy = qnorm(p)
points(qy, sy, col='red')

y = rgamma(100, 2, 3)
p = (1:length(y) - 0.5)/length(y)
qqnorm(y, xlim=c(-3.5,3.5), ylim=c(-3.5,3.5))
pars = par("usr")
abline(0,1)
sy = sort(y)
qy = qnorm(p)
points(qy, sy, col='red')
##########

b = 0.02998598

# GP example from rasmussen (process -> distribution)
library(MASS)
k = function(x, y=x, phi=0.5){
    out = matrix(0, length(x), length(y))
    for (i in 1:length(x))
        out[i, ] = exp(-phi*(x[i]-y)^2)
    return(out)
    }
m = function(x) 0.25 * x^2
good.dat = function(x) m(x) + rnorm(length(x), 0, 0.25)
bad.dat = function(x) -m(x) + rnorm(length(x), 0, 0.25)

# consider looking at different lengths for xx
x.data = seq(-5, 5, length=10)
x.test = seq(-5, 5, length=100)

mu.data = as.matrix(m(x.data))
mu.test = as.matrix(m(x.test))
sig.dd = k(x.data)
sig.dt = k(x.data, x.test)
sig.tt = k(x.test)

mm = m(x.data)
kk = k(x.data)
xx = x.data

yy = mvrnorm(10000, mm, kk)


up = apply(yy, 2, quantile, 0.975)
do = apply(yy, 2, quantile, 0.025)
plot(0, type='n', xlim=c(-5, 5), ylim=c(-3, 9))
# 95% probability interval
polygon(c(xx, xx[length(xx):1]), c(up, do[length(do):1]), col='skyblue')
# three draws
for (i in 1:3)
    lines(xx, yy[i,], lty=i)


# cholesky and inverses (of symmetric, positive-definite matrices)
A = matrix(c(
    10.7, 6.4, 1,
    6.4, 5.2, 4,
    1, 4, 9.1), 3, 3)
eigen(A)
B = chol(A) # B is upper triangle
# A^(-1) = (B'B)^(-1) = B^(-1)(B')^(-1) = B^(-1)(B^(-1))'
solve(A)
solve(t(B) %*% B)
solve(B) %*% solve(t(B))
solve(B) %*% t(solve(B)) # notice we only need to compute inverse of B


# singular value decomposition and PCA (see abdi) (incomplete)
f = function(x)
    c(x+rnorm(1), x^2+rnorm(1), x-3+rnorm(1), sqrt(abs(x))+rnorm(1), x^3+rnorm(1))
n = 100
A = matrix(f(rcauchy(n)), n, 5)
S = svd(A)
U = S$u
D = diag(S$d)
V = S$v
# variability portions
S$d/sum(S$d)
U %*% D %*% t(V) # equals A, but with some numerical fluctuations
V %*% t(V)  # equals I, V is orthogonal

F. = U %*% D # factor scores also = A %*% V

p = 2 # number of components
Am = U[,1:p] %*% D[1:p,1:p]/sqrt(n) # correct?

# making a covariance matrix
CS = apply(A, 2, function(x){ x-mean(x) })      # center
CS = apply(A, 2, function(x){ x/sqrt(n-1) })    # divide by root(n-1)
COV = t(CS) %*% CS

# making a correlation matrix
CC = apply(A, 2, function(x){ x-mean(x) })      # center
CC = apply(A, 2, function(x){ x/sqrt(sum(x^2)) })   # divide by norm of each variable
COR = t(CC) %*% CC


# cholesky and determinant relationship (symmetric positive-definite matrix)
A = matrix(c(
    10.7, 6.4, 1,
    6.4, 5.2, 4,
    1, 4, 9.1), 3, 3)
B = chol(A)
log(sqrt(det(A)))
sum(diag(log(B)))


# time comparison with cholesky decomposition
cov.func = function(x, y=x){
    out = matrix(0, length(x), length(y))
    for (i in 1:length(x))
        out[i, ] = exp(-sqrt((x[i]-y)^2))
    return (out)
    }

n = 10000
x = seq(0, 1, length=n)
A = cov.func(x)

# consider testing without assigning all variables with in system.time
# determinant
system.time(B <- chol(A))+
system.time(C <- sum(diag(log(B)))) # 196.542
system.time(D <- log(sqrt(det(A)))) # 283.960
# C and D should be equivalent (with high n, D might be -Inf)

# inverse (cholesky)
system.time(B <- chol(A))+
system.time(E <- solve(B))+
system.time(G <- E %*% t(E))
# 1484.692

# built-in R function
system.time(K <- chol2inv(A, n))
# 334.696?

# inverse (spectral decomposition)
system.time(I <- eigen(A, TRUE)) + 
system.time(J <- I$vectors %*% solve(diag(I$values)) %*% t(I$vectors))
# 

# straight solve
system.time(H <- solve(A))   # 1099.668
# 
# G, H, J, and K should be equivalent, perhaps with some numerical fluctuations


# outer and inner products (for one vector)
n = 100
V = runif(n, -10,10)
A = t(V) %*% V  # inner product
B = V %*% t(V)  # outer
# the next two are equivalent
sum(eigen(B)$values)
A


# what solves the equation x - exp(1) = log(x) ?
n = 1000
xx = seq(0, 5, length=n)
plot(xx, xx-exp(1), type='l', col='red')
lines(xx, log(xx), col='blue')
f1 = function(x) x - exp(1)
f2 = function(x) log(x)
f1(xx) - f2(xx)
b1 = 4.138652
f1(b1)
f2(b1)
b2 = 0.0708316
f1(b2)
f2(b2)


# welch, buck, sacks, wynn, mitchell, morris equation 6 example
y = function(x)
    (30 + x[1]*sin(x[1]))*(4+exp(-x[2]))
cov.func = function(w, x=w, theta, p=c(2,2)){
    out1 = matrix(0, length(x), length(x))
    out2 = matrix(0, length(x), length(x))
    for (i in 1:length(w)){
        out1[i, ] = exp(-theta[1]*abs(w[i]-x)^p[1])
        out2[i, ] = exp(-theta[2]*abs(w[i]-x)^p[2])
        }
    return (out1 * out2)
    }
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
xx = seq(0, 5, length=100)
xx = cross(xx, xx)
zz = apply(xx, 1, y)
library(rgl)
plot3d(xx[,1], xx[,2], zz)


# absolute cauchy distribution
n = 100000
X = abs(rcauchy(n))
hist(X, breaks=c(seq(0, 10, length=100), seq(10.1, max(X), length=2)), xlim=c(0, 10), freq=F)
curve(dexp(x), add=T, col='blue')


# schur complement method (incomplete)
A = matrix(rnorm(9), 3, 3)
B = cbind(A, rnorm(3))
B = rbind(B, rnorm(4))
C = solve(A)
D = solve(B)


# simple latin hypercubes
lh = function(n, dimensions, domain=c(0, 1)){
    d = dimensions
    r = domain
    out = matrix(seq(r[1], r[2], length=n), n, d)
    for (i in 1:d)
        out[,i] = sample(out[,i])
    return (out)
    }
pairs(out<-lh(11, 5), xlim=c(0, 1), ylim=c(0, 1))
plot(lh(10,2))

# minimax (johnson 1990), doesn't cover boundaries like maximin
# d = 1
n = 13
S.star = seq(1/(2*n), (2*n-1)/(2*n), length=n)
d.star = 1/(2*n)
# d = 2
n = 3
S.star = matrix(0, 3, 2)
S.star[1,] = c(1/16, 9/16)
S.star[2,] = c(1/2, 1/2-1/4)
S.star[3,] = c(1/2, 1/2+1/4)
d.star = sqrt(65)/16


# cross-validation
# data generation
f = function(x)
    3 + 2*x[1] - 6*x[2]
n = 1000
design = matrix(c(runif(n, 2, 6), runif(n, -1, 3)), n, 2)
y = apply(design, 1, f) + rnorm(n)

# subset
val.size = floor(n/2)
reps = 1000
mse.vec = double(reps)
for (i in 1:reps){
    val.samp = sort(sample(n, val.size))
    predict = (1:n)[-val.samp]
    mod = lm(y[val.samp] ~ design[val.samp,])
    mse.vec[i] = mean((cbind(1,design[predict,])%*%mod$coef-y[predict])^2)
    }
mean(mse.vec)

# LOOCV (leave one out cross-validation)
mse.vec2 = double(n)
for (i in 1:n){
    mod = lm(y[-i] ~ design[-i,])
    mse.vec2[i] = (c(1,design[i,])%*%mod$coef-y[i])^2
    }
mean(mse.vec2)


########## polynomial fitting
design = seq(0, 2, length=12)
y = design+rnorm(length(design), sd=0.25)
plot(design, y, ylim=c(-1, 3))

xx = seq(0, 2, length=100)

mod1 = lm(y ~ design)
pred1 = cbind(1, xx) %*% mod1$coef
lines(xx, pred1, lwd=2)

newD = design
newX = cbind(1, xx)
for (i in 2:length(design)){
    readline(paste(i,"/",length(design)))
    newD = cbind(newD, design^i)
    newX = cbind(newX, xx^i)
    mod = lm(y ~ newD)
    pred = newX %*% mod$coef
    lines(xx, pred, lwd=2, col=i)
    }
# another example with not really a polynomial
x = seq(-5, 10, length=1000)
y = x/2 + sin(2*x) + rnorm(length(x), sd=0.1)
plot(x, y, pch=20)

max.P = 20
X = rep(1, length(x))
for (i in 1:max.P){
    readline(paste(i,"/",max.P))
    X = cbind(X, x^i)
    mod = lm(y ~ -1 + X)
    lines(x, X %*% coef(mod), lwd=3, col=rainbow(max.P)[i])
    }
##########


# kronecker product (o times)
A = matrix(rnorm(9),3,3)
kronecker(A, diag(2))
kronecker(diag(2), A)



xx = seq(0, 0.050, length=100000)
theta = seq(0.001, 0.03, length=100)
plot(xx, dnorm(xx, theta[1], theta[1]), type='l', xlim=c(0,0.04))
for (i in 2:length(theta))
    points(xx, dnorm(xx, theta[i], theta[i]), type='l')


X = matrix(rnorm(20), 10, 2)
for (i in 1:2){
    a = min(X[,i])
    b = max(X[,i])
    X[,i] = X[,i] - a - 1


X[,2] = X[,2] - 
X[,1] %*% solve(t(X[,1]) %*% X[,1]) %*% t(X[,1]) %*% X[,2]



# colors - 21 ^j
dat = read.table("~/files/R/rgb.txt")
index = NULL
for (i in 1:(nrow(dat)-1)){
    if (dat[i,1] == dat[i+1,1] && dat[i,2] == dat[i+1,2] && dat[i,3] == dat[i+1,3])
        index = c(index, i)
    }
dat = dat[-index,]
label = as.character(dat[,4])
n = length(label)

k = 5 # number of columns in par
pdf("./colors.pdf", height=ceiling(n/k)/9)
par(mfrow=c(ceiling(n/k), k), oma=double(4), mar=double(4))
for (i in 1:n){
    plot(0, type='n', xlab="", ylab="", xaxt='n', yaxt='n')
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
        col = rgb(dat[i,1], dat[i,2], dat[i,3], maxColorValue=255))
    text(0, labels=label[i])
    }
dev.off()

# alternative
all = colors()
n = length(all)
k = 5 # number of columns in par
pdf("./colors2.pdf", height=ceiling(n/k)/8)
par(mfrow=c(ceiling(n/k), k), oma=double(4), mar=double(4))
for (i in 1:n){
    plot(0, type='n', xlab="", ylab="", xaxt='n', yaxt='n')
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
        col = all[i])
    text(0, labels=all[i])
    }
dev.off()



# simple gaussian process for plotting
library(MASS)
n = 15
dat.func = function(x)
    -1*(x-0.5)^2+3*(x-0.5)^4 + rnorm(length(x), 0, 1/128)
dat.x = c(0,1,runif(n-2))
dat.y = dat.func(dat.x)
dat.x = c(0.49479955, 0.54777443, 0.48692351, 0.56489071, 0.08003327,
    0.73539188, 0.69787756, 0.02739255, 0.40758957, 0.21963842,
    0.92037507, 0.50282820, 0.51918529, 0.73909302, 0.15778479)
dat.y = c(-0.0005333700, -0.0039474049, -0.0006762263, -0.0056543233,
    -0.0868673788, -0.0517171887, -0.0361714292, -0.0753287717,
    -0.0074365065, -0.0593913496, -0.0842185453, 0.0003879083,
    -0.0023916517, -0.0485891627, -0.0755961365)
plot(dat.x, dat.y)

cov.func = function(x, y=x, phi){
    out = matrix(0, length(x), length(y))
    for (i in 1:length(x))
        out[i, ] = exp(-phi*(x[i]-y)^2)
    return (out)
    }

phi = 2
eps = 0.00025
m = 100
test.x = seq(min(dat.x)-0.05, max(dat.x)+0.05, length=m)

sig.xx = cov.func(dat.x, dat.x, phi)
sig.yy = cov.func(test.x, test.x, phi)
sig.xy = cov.func(dat.x, test.x, phi)

cov.inv = solve(sig.xx + diag(eps, n))
mu.mat = t(sig.xy) %*% cov.inv %*% dat.y
sig.mat = sig.yy - t(sig.xy) %*% cov.inv %*% sig.xy

yfit = mvrnorm(10000, mu.mat, sig.mat)
gpmean = apply(yfit, 2, mean)
gpup95 = apply(yfit, 2, quantile, 0.975)
gplo95 = apply(yfit, 2, quantile, 0.025)

pdf("~/files/afrl_project/slides/figs/gp_example.pdf", height=5, width=7)
plot(dat.x, dat.y, col='forestgreen', lwd=4, xlim=c(min(test.x)+0.04,
    max(test.x)-0.04), ylim=c(min(gplo95), max(gpup95)), pch=20,
    xlab="", ylab="")
grid(lty=2)
polygon(c(test.x, test.x[m:1]), c(gpup95, gplo95[m:1]),
    col='cadetblue1', border='gray')
#for (i in 1:100)
#    points(test.x, yfit[i,], col='gray', type='l')
points(dat.x, dat.y, col='yellow1', lwd=12, pch=20)
points(dat.x, dat.y, col='yellow3', lwd=2, cex=2.3)
lines(test.x, gpmean, lwd=2)
dev.off()


### rescale function
transform = function(x, new.scale, old.scale=range(x)){
    (x - old.scale[1])*(new.scale[2]-new.scale[1])/
        (old.scale[2]-old.scale[1])+new.scale[1]
    }


### hexagonal lattice in two dimensions
n1 = 6
n2 = 4
n = n1 * n2

design = matrix(0, n, 2)
for (i in 1:n){
    design[i, 1] = ( (i-1) %% n1  ) + 1/2*(floor((i-1)/n1) %% 2)
    design[i, 2] = sqrt(3)/2*floor( (i-1)/n1 )
    }
Xrange = c(0, (n1-1) + 0.5)
Yrange = c(0, sqrt(3)/2*(n2-1))

plot(design)

# rescale (d1: x in [0,1]; d2: y in [0,1])
new.d1 = design %*% diag(1/Xrange[2], 2)
new.d2 = design %*% diag(1/Yrange[2], 2)
points(new.d1, col='red')
points(new.d2, col='blue')

# reposition
new.d3 = new.d1 + rep(c(4, 1), each=n)
new.d4 = new.d2 + rep(c(2, 1.5), each=n)
points(new.d3, col='green')
points(new.d4, col='orange')



### alpha levels, p-values
y.gen = function(x)
    0 + 0*x + rnorm(length(x))

n = 100
m = 10000
xx = seq(0, 1, length=n)
results = matrix(0, m, 2)
for (i in 1:m){
    yy = y.gen(xx)
    results[i,] = summary(lm(yy ~ xx))[[4]][,4]
    }
mean(results[,1] < 0.05)
mean(results[,2] < 0.05)
hist(results[,1])


########## log-scale and very large values
#
# bi-symmetric log transformation
trans = function(x, center=median(x), base=exp(1))
    sign(x-center)*log(abs(x-center)+1, base)

# some views of the function
xx = seq(-5, 5, length=1000,)
plot(xx, trans(xx, base=exp(1)), type='l')
points(xx, trans(xx, base=2), type='l', col='green')
points(xx, trans(xx, base=10), type='l', col='red')
points(xx, xx, lty=2, col='gray', type='l')

# large positive values, not as a extreme negative, many close to 0
n.good = 100
n.high = 5
n.low = 5
Z = c(rnorm(n.good, 2), rnorm(n.high, 150, 25), rnorm(n.low, -40, 5))
ord = sample(n.good+n.high+n.low)
pchs = c(rep(1, n.good), rep(20, n.high+n.low))[ord]
cols = ifelse(Z>=median(Z), "red", "blue")[ord]
Z = Z[ord]

par(mfrow=c(3,1), mar=c(2,4.1,2,2.1))
plot(Z, pch=pchs, col=cols)

# zoomed in
# plot(Z, pch=pchs, col=cols, ylim=range(Z[ord<=n.good]))
# abline(h=0, col='gray', lty=2)

# log-ish transformation centered around 0
D = trans(Z, center=0)
plot(D, pch=pchs, col=cols)
abline(h=0, col='gray', lty=2)

# centered around the median (gives a better spread of middle points)
E = trans(Z, median(Z))
plot(E, pch=pchs, col=cols)
abline(h=0, col='gray', lty=2)



# only positive values (comparison of bi-symmetric log transformation
# and a regular log transform, plots 2 and 3, respectively)
n.good = 100
n.high = 10
Z = c(rgamma(n.good,1,1), rgamma(n.high, 5, 0.1))
ord = sample(n.good+n.high)
pchs = c(rep(1, n.good), rep(20, n.high))[ord]
cols = ifelse(Z>=median(Z), "red", "blue")[ord]
Z = Z[ord]

par(mfrow=c(3,1), mar=c(2,4.1,2,2.1))
plot(Z, pch=pchs, col=cols)
abline(h=median(Z), col='gray', lty=2)

D = trans(Z)
plot(D, pch=pchs, col=cols)
abline(h=trans(median(Z)), col='gray', lty=2)

E = log(Z)
plot(E, pch=pchs, col=cols)
abline(h=log(median(Z)), col='gray', lty=2)
#
##########



########## arcsin(sqrt(unif(0, 1))) approx normal
#
par(mfrow=c(1,1))
S = asin(sqrt(runif(1000000)))
plot(density(S))
curve(dnorm(x, pi/4, pi/8), add=TRUE, col='red')
#
##########



##########
#
# 642 9.2
theta = 10
n = 100
probs = double(1000)
for (i in 1:length(probs)){
    X = rnorm(n, theta, 1)
    int = mean(X)+c(-1,1)*qnorm(0.975)/sqrt(n)
    temps = rnorm(1000, theta, 1)
    probs[i] = mean(temps >= int[1] & temps <= int[2])
    }

#
##########


##########
#
# some intuition behind likelihood and M-H
true.alpha = 4
true.beta = 18

n = 1000
dat = rbeta(n, true.alpha, true.beta) # log posterior: 1155.545
dat = rnorm(n, 0.5, 0.1) # bad likelihood, log post: 892.0784
hist(dat, breaks=25, col="gray", freq=FALSE, xlim=c(0,1))

calc.post = function(params)
    sum(dbeta(dat, params[1], params[2], log=TRUE))+
        dgamma(params[1], a.alpha, b.alpha, log=TRUE)+
        dgamma(params[2], a.beta, b.beta, log=TRUE)


a.alpha = 1
b.alpha = 0.075
a.beta = 1
b.beta = 0.075

sig = c(0.7,1.1)

nburn = 100
nmcmc = 500
params = matrix(0, nburn+nmcmc, 2)
params[1,] = c(55, 0.1) # starting values far from true

accept = c(0,0)

for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:2){
        hist(dat, breaks=25, col="gray", freq=FALSE, xlim=c(0,1))
        legend("topright", c("Current", "Proposal"), col=c('blue','green'),
            lty=1, lwd=2)
        curve(dbeta(x, params[i,1], params[i,2]), add=TRUE, col='blue', lwd=2)
        post = calc.post(params[i,])
        cand = rnorm(1, params[i,j], sig[j])
        if (cand > 0){
            cand.param = params[i,]
            cand.param[j] = cand
            curve(dbeta(x, cand.param[1], cand.param[2]), add=TRUE, col='green', lwd=2)
            cand.post = calc.post(cand.param)
            if (log(runif(1)) < cand.post - post){
                legend("bottomright", "ACCEPT", box.lwd=0, cex=1.5)
                params[i,j] = cand
                if (i > nburn)
                    accept[j] = accept[j] + 1
                }
            }
        }
#   Sys.sleep(1)
    }
post
cand.post

accept/nmcmc
#
##########


##########
#
# logit-beta
f = function(y, a, b)
    1/beta(a,b)*1/(1+exp(-y))^a*(exp(-y)/(1+exp(-y)))^b
a <- b <- 2.234756201
a = 5
b = 2
f(0, a, b) > dnorm(0)

yy = seq(-4, 4, length=1000)
plot(yy, f(yy, a, b), type='l', ylim=c(0,0.4))
#points(yy, dbeta(yy, a, b), type='l', col='red')
points(yy, dnorm(yy), type='l', col='green')
#points(yy, dcauchy(yy), type='l', col='blue')
#

for (a in c(0.1, 0.5, 1, 2, 5, 10, 50, 100)){
    for (b in c(0.1, 0.5, 1, 2, 5, 10, 50, 100)){
        yy = seq(-50, 50, length=1000)
        plot(yy, f(yy, a, b), type='l', main=paste0("a=", a, " -- ",
            "b=",b))
#       points(yy, dcauchy(yy), type='l', col='blue')
        readline()
        }
    }
##########


##########
#
# moving average
n = 250
xx = seq(0, 1, length=n)
yy = double(n)
yy[1] = rnorm(1, 0, 1)
for (i in 2:n)
    yy[i] = rnorm(1, yy[i-1], 1)
plot(xx, yy, type='l')

prop = c(0.05, 0.10, 0.15)
cols = c("red", "blue", "green")
for (i in 1:3){
    sub.n = ceiling(n * prop[i])
    ma.x = seq(min(xx), max(xx), length=n-sub.n)
    ma.y = double(length(ma.x))
    for (j in 1:length(ma.x))
        ma.y[j] = mean(yy[j:(j+sub.n-1)])
    points(ma.x, ma.y, col=cols[i], lwd=2, type='l')
    }
legend("topleft", legend=prop, col=cols, lty=1, lwd=2,
    bg="white")
#
##########


##########
#
# moving average-based mcmc auto tuning
# just a first look at random bernoulli trials
# and moving average
n = 1000
xx = seq(0, 1, length=n)
yy = rbinom(n, 1, 0.5)
plot(xx, yy, pch=20)

# prop is proportional to window sizes
prop = c(0.05, 0.10, 0.15)
cols = c("red", "blue", "green")
for (i in 1:3){
    sub.n = ceiling(n * prop[i])
    ma.x = seq(min(xx), max(xx), length=n-sub.n)
    ma.y = double(length(ma.x))
    for (j in 1:length(ma.x))
        ma.y[j] = mean(yy[j:(j+sub.n-1)])
    points(ma.x, ma.y, col=cols[i], lwd=2, type='l')
    }
legend("topleft", legend=prop, col=cols, lty=1, lwd=2,
    bg="white")

# the auto tuning part, goal: a function only of
# the acceptance rate

n = 1000
# true values
a.alpha = 5
b.alpha = 3
a.beta = 2
b.beta = 4
dat = rgamma(n, rgamma(n, a.alpha, b.alpha),
    rgamma(n, a.beta, b.beta))
hist(dat, breaks=100, col='gray', freq=FALSE)

# estimating alpha, beta, a.alpha, b.alpha, a.beta, b.beta
# hyperparameters
a = 1
b = 0.01
calc.post = function(params)
    sum(dgamma(dat, dgamma(params[1], params[3], params[4]),  # likelihood
        dgamma(params[2], params[5], params[6]), log=TRUE)) + 
        dgamma(params[1], params[3], params[4], log=TRUE) +   # alpha prior
        dgamma(params[2], params[5], params[6], log=TRUE) +   # beta prior
        sum(dgamma(params[3:6], a, b, log=TRUE))              # a.alpha, etc prior

autotune = function(accept, target, k){
    k = c((1-1/k)/(cosh(target)-1), (k-1)/(cosh(target-1)-1))
    1+sign(accept-target)*(cosh(accept-target)-1)*
        k[(sign(accept-target)+1)/2+1]
    }

 m = 0.25
 k = 2
 xx = seq(0, 1, length=100)
 plot(xx, autotune(xx, m, k), type='l')
# points(xx, 1/f(xx, m, k), type='l', col='red')
# range(f(xx, m, k))

new.seq = function(x, back){
    out = seq(x-back+1, x)
    return (out[out > 0])
    }

# a = 0.2
# b = 0.35
# L = 0.01
# U = 2.9
# ll = seq(0, a, length=100)
# uu = seq(b, 1, length=100)
# lower = function(x)
#     L+1-exp(x/a*log(L))
# upper = function(x)
#     exp((x-b)/(1-b)*log(U))
# plot(ll, lower(ll), type='l', ylim=c(L, U), col='blue', xlim=c(0,1))
# points(uu, upper(uu), type='l', col='red')
# points(seq(a, b, length=10), double(10)+1, type='l')

nburn = 50000
window = 500 # how often the autotune function is used
nmcmc = 10000
params = matrix(0, nburn+nmcmc, 6)
params[1,] = c(10, 10, 10, 10, 10, 10)

sig = c(100, 100, 100, 0.0001, 0.0001, 0.0001)
accept = matrix(0, nburn+nmcmc, 6)
plot.me = matrix(0, ceiling((nburn+nmcmc)/window), 6)

keep.sig = matrix(0, 1+nrow(plot.me), 6)
keep.sig[1,] = sig

par(mfrow=c(6,2), mar=double(4)+1.5, oma=double(4)+0.5)
for (i in 2:(nburn+nmcmc)){
    params[i,] = params[i-1,]
    for (j in 1:6){
        post = calc.post(params[i,])
        cand = rnorm(1, params[i,j], sig[j])
        if (cand > 0){
            cand.param = params[i,]
            cand.param[j] = cand
            cand.post = calc.post(cand.param)
            if (log(runif(1)) < cand.post - post){
                params[i,j] = cand
                accept[i,j] = 1 #accept[i,j] + 1
                }
            }
        if (floor(i/window) == i/window && i < nburn) # every window'th iteration run autotune function
            sig[j] = autotune(sig[j], mean(accept[(i-window+1):i,j]))
        }
    if (floor(i/window) == i/window){
        plot.me[i/window, ] = apply(accept[(i-window+1):i,], 2, mean)
        keep.sig[i/window+1,] = sig
        for (j in 1:6){
            plot(plot.me[1:(i/window),j], type='l', ylim=c(0,1))
            abline(h=0.25, lty=2, col='gray')
            plot(keep.sig[new.seq(i/window+1, 50),j], type='l', col='blue')
            }
        }
    }
apply(accept[(nburn+1):(nburn+nmcmc),], 2, mean)

par(mfrow=c(2,3))
xx = seq(1, nburn+nmcmc)
for (h in 1:6){
    plot(xx, accept[,h], pch=20)
    prop = c(0.01, 0.05, 0.10)
    cols = c("red", "blue", "green")
    for (i in 1:3){
        sub.n = ceiling((nburn+nmcmc) * prop[i])
        ma.x = seq(min(xx), max(xx), length=nburn+nmcmc-sub.n)
        ma.y = double(length(ma.x))
        for (j in 1:length(ma.x))
            ma.y[j] = mean(accept[j:(j+sub.n-1), h])
        points(ma.x, ma.y, col=cols[i], lwd=2, type='l')
        }
    }

burn = params[(nburn+1):(nburn+nmcmc),]
par(mfrow=c(2,3))
for (j in 1:6) 
    plot(burn[,j], type='l')
#
##########


##########
# recursion
multiply = function(x, i){
    if (i > 1)
        return (x + multiply(x, i-1))
    return (x)
    }

multiply(4, 3)
multiply(5, 2)
#
##########

##########
# congruential generator on {0, 1, ..., M}
# Note 2.6.1 of Robert and Casella "Monte Carlo Statistical Methods"
f = function(x) (69069*x) %% 1
xx = seq(0, 1, length=2000)
yy = f(xx)
plot(yy, xx, pch=20, col='gray')
for (i in 1:length(xx))
    points(yy[i], xx[i], pch=20)

# the Kiss Algorithm
k = 30
L = matrix(0, k, k)
R = matrix(0, k, k)
I = diag(k)
for (i in 2:k){
    L[i-1, i] = 1
    R[i, i-1] = 1
    }
L13 = L %*% L %*% L %*% L %*% L %*%
        L %*% L %*% L %*% L %*% L %*%
        L %*% L %*% L
L15 = L %*% L %*% L %*% L %*% L %*%
        L %*% L %*% L %*% L %*% L %*%
        L %*% L %*% L %*% L %*% L
R17 = R %*% R %*% R %*% R %*% R %*%
        R %*% R %*% R %*% R %*% R %*%
        R %*% R %*% R %*% R %*% R %*%
        R %*% R
R18 = R %*% R %*% R %*% R %*% R %*%
        R %*% R %*% R %*% R %*% R %*%
        R %*% R %*% R %*% R %*% R %*%
        R %*% R %*% R
two.32 = 2^32
two.31 = 2^31
kiss = function(i, j, k){
    i = (69069 * i + 23606797) %% (two.32)
    j = ((I + L15) %*% (I + R17) * j) %% (two.32)
    k = ((I + L13) %*% (I + R18) * k) %% (two.31)
    x = (i + j + k) %% (two.32)
    return (c(x, i, j, k))
    }
kiss(1,1,1)
kiss(23675866, 71, 71)
#
##########

##########
#
# plots of casella and berger note 2.6.3

f1 = function(x) dnorm(x)
f2 = function(x) f1(x)*(1+0.5*sin(2*pi*x))
K1 = function(t) (t^2)/2
K2 = function(t) K1(t) + log(1+0.5*exp(-2*pi^2)*sin(2*pi*t))

par(mfrow=c(2,1), mar=c(3.1,3.1,2.1,2.1))
xx = seq(-4, 4, length=10000)
tt = xx
plot(xx, f2(xx), type='l', col='blue', lwd=2)
points(xx, f1(xx), type='l', col='red', lwd=2)
plot(xx, K1(tt), type='l', col='red', lwd=2, lty=2)
points(xx, K2(tt), type='l', col='blue', lwd=2, lty=3)
#
##########

##########
#
# cauchy and normal ratios
n = 10000
X = rnorm(n)
Y = rnorm(n)
C = X/Y
plot(density(C, n=5*n), xlim=c(-5,5))
curve(dcauchy(x), add=TRUE, col='red')

n = 100000
X = rnorm(n, 0, sqrt(6))
Y = rnorm(n, 0, sqrt(2))
C = X/Y
plot(density(C, n=10*n), xlim=c(-5,5))
curve(dcauchy(x, 0, sqrt(3)), add=TRUE, col='red')
#
##########

##########
#
x = rnorm(10000, 3, sqrt(3))
mean(x)
yy = seq(-20, 20, length=10000)
mm = sum(yy*dnorm(yy, 3, sqrt(3))) * (yy[2] - yy[1])
vv = sum(yy*dnorm(yy, 3, sqrt(3))) * (yy[2] - yy[1]) - mm^2
var(x)

n = 100
m = 10000
out = matrix(0, m, n)
for (i in 1:n){
    for (j in 1:m){
        x = rbeta(i, 1, 1)
        out[j,i] = mean(x^(1/x))
        }
    }
up = apply(out, 2, quantile, 0.975)
do = apply(out, 2, quantile, 0.025)
me = apply(out, 2, mean)
maxes = apply(out, 2, max)
mines = apply(out, 2, min)
plot(me, pch=20, ylim=c(0,1))
points(up, type='l')
points(do, type='l')
points(maxes, type='l', col='red')
points(mines, type='l', col='red')
#
##########

##########
#
# regression on a function of the y's
set.seed(1)
predx = cbind(1, seq(0, 10, length=100))
x = runif(15, 0, 10)
b0 = 2
b1 = 1.2
sig2 = 0.8
y = b0 + b1*x + rnorm(length(x), 0, sqrt(sig2))
plot(x, y, pch=20)

xstar = cbind(1, x)
est = solve(t(xstar) %*% xstar) %*% t(xstar) %*% y
predy = predx %*% est
points(predx[,2], predy, type='l', lty=2, col='gray30')

# log transformation
fy = log(y)
plot(x, fy, pch=20)
fest = solve(t(xstar) %*% xstar) %*% t(xstar) %*% fy
predy = predx %*% fest
points(predx[,2], predy, type='l', lty=2, col='gray30')
# notice we get different estimates

plot(x, y, pch=20)
points(predx[,2], exp(predy), type='l', lty=2, col='gray30')
# and inaccurate predictions on the original scale, but
# this may because of our model choice: the linear model
# could not be appropriate in this case

# with a gaussian process
# no transformation
library(geoR)

set.seed(1)
x = runif(15, 0, 10)
b0 = 2
b1 = 1.2
sig2 = 0.8
y = b0 + b1*x + rnorm(length(x), 0, sqrt(sig2))
orig.y = y
plot(x, y, pch=20)

nu = 2
field.dat = as.geodata(cbind(y, x, 0),  data.col=1, coords.col=c(2,3))
gp.field = likfit(field.dat, cov.model="matern", kappa=nu,
    fix.kappa=TRUE, ini.cov.pars=c(10, 10), trend="cte",
    messages=FALSE)
phi = 1/gp.field$phi
s2 = gp.field$sigmasq
beta = as.matrix(gp.field$beta)
tau2 = gp.field$tausq

predx = seq(0, 10, length=100)
N = length(x)
K = length(predx)
D = rdist(c(predx, x))
V = s2*Matern(D, alpha=phi, nu=nu)
inv.V = solve(V[K+(1:N),K+(1:N)] + diag(tau2, N))
pred.EV = predx %*% beta + V[1:K,K+(1:N)] %*% inv.V %*%
    (y - x %*% beta)
plot(x, y, pch=20)
points(predx, pred.EV, type='l')

# log transformation
y = log(y)
nu = 2
field.dat = as.geodata(cbind(y, x, 0),  data.col=1, coords.col=c(2,3))
gp.field = likfit(field.dat, cov.model="matern", kappa=nu,
    fix.kappa=TRUE, ini.cov.pars=c(10, 10), trend="cte",
    messages=FALSE)
phi = 1/gp.field$phi
s2 = gp.field$sigmasq
beta = as.matrix(gp.field$beta)
tau2 = gp.field$tausq

predx = seq(0, 10, length=100)
N = length(x)
K = length(predx)
D = rdist(c(predx, x))
V = s2*Matern(D, alpha=phi, nu=nu)
inv.V = solve(V[K+(1:N),K+(1:N)] + diag(tau2, N))
pred.log.EV = predx %*% beta + V[1:K,K+(1:N)] %*% inv.V %*%
    (y - x %*% beta)
plot(x, y, pch=20)
points(predx, pred.log.EV, type='l')

plot(x, orig.y, pch=20)
points(predx, pred.EV, type='l')
points(predx, exp(pred.log.EV), type='l', col='red')

# there is some discrepancy between the predictions made on the
# untransformed response with the transformed response, but when
# the predictions from the transformed response are brought
# back to the original scale, the model still does pretty good

##########


##########
#
# computing the modified bessel function of the second kind

f = function(x, nu, maxk = 20){
    calc.i = function(x, nu){
        out = 0
        for (k in 0:maxk)
            out = out + ((1/4 * x^2)^k) / (factorial(k) * 
                gamma(nu + k + 1))
        out = (1/2 * x)^(nu) * out
        return (out)
        }
    return ((pi/2) * (calc.i(x, -nu) - calc.i(x, nu)) / sin(nu * pi))
    }

f(3, 1.5, 10)
# built-in R
besselK(3, 1.5)

pi/2*(besselI(3, -1.5) - besselI(3, 1.5))/sin(pi*1.5)
##########

##########
# ranjan (2011) example (the discrepancy function doesn't
# match with Figure 4.1, so I use other values)
eta = function(x){
    x1 = x[1]
    x2 = x[2]
    t = x[3]
    (30 + 5*x1*sin(5*x1))*(6*t+1+exp(-5*x2))
    }
disc = function(x){
    x1 = x[1]
    x2 = x[2]
    -50*exp(-1.0*x1 - 0.5*x2)
    }

x1 = seq(0, 1, length=100)
x2 = seq(0, 1, length=100)
t = 0.5

Zeta = matrix(0, length(x1), length(x2))
Zdisc = matrix(0, length(x1), length(x2))
for (i in 1:length(x1)){
    for (j in 1:length(x2)){
        Zeta[i,j] = eta(c(x1[i], x2[j], t))
        Zdisc[i,j] = disc(c(x1[i], x2[j]))
        }
    }

par(mfrow=c(1,2))
contour(x1, x2, Zeta, col=topo.colors(20))
contour(x1, x2, Zeta+Zdisc, col=topo.colors(20))
##########

##########
smooth = function(x, y, m, d){
    # m = number of points to do the smoothing,
    # should be greater than length(x) to look nice
    if (missing(m))
        m = max(100, 8*length(x))
    if (missing(d))
        d = diff(range(x)) / length(x)

    # gaussian kernel, delta is the bandwidth parameter
    kern = function(x, y)
        exp(-1/(2*d^2)*(x-y)^2)

    # go outside the range of the data by some
    w = diff(range(x))*0.25

    outx = seq(min(x) - w, max(x) + w, length=m)
    outy = double(m)
    # loop to compute the weighted value at each new x
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i]))/
            sum(kern(x, outx[i]))
    return (list("x"=outx, "y"=outy))
    }

# 1-d
n = 30
m = max(100, n*8)
daty = rgamma(n, 2, 1)

datx = seq(0, 100, length=n)
r = diff(range(datx))
w = diff(range(datx))*0.25

plot(datx, daty, type='h', xlim=c(min(datx)-w, max(datx)+w))
ss = smooth(datx, daty)
points(ss$x, ss$y, type='l')

# higdon (2001) 4.1 (the function, not the multiresolution stuff)
set.seed(2)
n = 30
gen.y = function(x)
    sin(2*pi*x/10) + 0.2*sin(2*pi*x/2.5) + rnorm(length(x), 0, 0.1)
xx = seq(1, 10, length=n)
yy = gen.y(xx)
plot(xx, yy, pch=20)
ss = smooth(xx, yy)
points(ss$x, ss$y, type='l')
##########

##########
# choleksy decomposition algorithm
# from p.27 of rencher and christensen
cov.func = function(x, y=x){
    out = matrix(0, length(x), length(y))
    for (i in 1:length(x))
        out[i, ] = exp(-sqrt((x[i]-y)^2))
    return (out)
    }
xx = seq(0, 1, length=6)
A = cov.func(xx)
B = chol(A)

mw.chol = function(x){
    n = nrow(x)
    out = matrix(0, n, n)
    out[1,1] = sqrt(x[1,1])
    for (j in 2:n)
        out[1,j] = x[1,j] / out[1,1]
    for (i in 2:n){
        for (j in i:n){
            if (i == j){
                out[i,i] = sqrt(x[i,i] - sum(out[1:(i-1),i]^2))
            } else {
                out[i,j] = (x[i,j] - sum(out[1:(i-1),i] * out[1:(i-1),j])) / out[i,i]
                }
            }
        }
    return (out)
    }

xx = seq(2, 5, length=6)
A = cov.func(xx)
(B = chol(A))
mw.chol(A)

##########
n = 2000
d = 2000
x = matrix(runif(n*d), n, d)
j = rep(1, n)

# overall sum
system.time(sum(x))
system.time(t(j) %*% x %*% j)
# overall mean
system.time(mean(x))
system.time(t(j/n) %*% x %*% (j/n))

# row sum
system.time(x %*% j)
system.time(apply(x, 1, sum))
# row mean
system.time(x %*% (j/n))
system.time(apply(x, 1, mean))

# column sum
system.time(t(j) %*% x)
system.time(apply(x, 2, sum))
# column mean
system.time(t(j/n) %*% x)
system.time(apply(x, 2, mean))


##########
# gamma regression? (epsilon distributed gamma(a, b))

n = 100000
true.beta = c(2, 0.5, -1)
true.a = 1.2
true.b = 2.5
X = matrix(runif(n*2, -5, 10), n, 2)
Y = cbind(1, X) %*% true.beta + (rgamma(n, true.a, scale=true.b) - true.a * true.b)
#Y = cbind(1, X) %*% true.beta + rnorm(n, 0, sqrt(true.a*true.b^2))

m = floor(n * 0.5)
train = sort(sample(n, m))
test = sort((1:n)[-train])

dat = data.frame(cbind(Y, X))
names(dat) = c("Y", "X1", "X2")

mod = lm(Y ~ 1 + X1 + X2, data = dat[train,])
# try glm(..., family = gamma)?
plot(fitted(mod), rstudent(mod), pch=20)

coef(mod) # the betas are estimated fine, but how to estimate a and b?
preds = predict(mod, newdata = dat[test,], interval="prediction")

mean(preds[,2] <= Y[test,]) # should be about 0.975 for both
mean(Y[test,] <= preds[,3])

mean(preds[,2] <= Y[test,] & Y[test,] <= preds[,3]) # should be about 0.95




##########

##########
# logistic, sum to one
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
logit = function(x)
    log(x/(1-x))

xx = seq(0.01, 0.99, by=0.01)
yy.cross = cross(xx, xx)
yy.sum = yy.cross
for (i in 1:nrow(yy.sum))
    yy.sum[i,] = yy.sum[i,] / sum(yy.sum[i,])

zz.cross = logit(yy.cross)
zz.sum = logit(yy.sum)

plot(zz.cross[1:99,1], type='l', col='red', ylim=c(-5,5))
points(zz.cross[1:99,2], type='l', col='blue')

plot(zz.sum[1:99,1], type='l', col='red', ylim=c(-5,5))
points(zz.sum[1:99,2], type='l', col='blue')
##########

##########
# multiple tests
# getting "significant" results when that is not true
set.seed(1)
m = 1000
alpha = double(m)
minp = double(m)
maxp = double(m)

for (i in 1:m){
    n = 5000
    d = 100
    X = cbind(1, matrix(runif(n*(d-1)), n, d-1))
    # Y is distribution XB + epsilon, B = 0 and
    # epsilon ~ N(0, 1)
    Y = rnorm(n)

    beta = solve(t(X) %*% X) %*% t(X) %*% Y
    s2 = as.vector(1/(n-d) * t(Y - X %*% beta) %*% (Y - X %*% beta))
    stderr = sqrt(s2 * diag(solve(t(X) %*% X)))
    t.stat = beta/stderr
    p.vals = 2*pt(abs(t.stat), n-d, lower.tail=FALSE)

    alpha[i] = mean(p.vals <= 0.05)
    minp[i] = min(p.vals)
    maxp[i] = max(p.vals)
    }
alpha
range(alpha)
mean(alpha)
hist(alpha)
names(which.max(table(alpha)))
###########

###########
set.seed(1)
n = 10
p = 5
Y = matrix(rnorm((n+1)*p), n+1, p)

X = Y[1:n,]
z = matrix(Y[n+1,], p, 1)

X.cov = cov(X)
# current inverse
X.inv = solve(X.cov)

Y.cov = cov(Y)
# target
Y.inv = solve(Y.cov)

# B + zz'
z = sqrt(n+1)*(z - apply(Y, 2, mean)) / n

(z - apply(Y, 2, mean))/d

# a check, z should equal d, which produces correct result
D = Y.cov - (n-1)/n*X.cov
d = matrix(sqrt(diag(D)), p, 1) * sign(D[1,])


(n-1)/n*X.cov + d %*% t(d)
Y.cov

new.inv = n/(n-1)*X.inv

# should all be equal
Updated = new.inv - (new.inv %*% z %*% t(z) %*% new.inv) / as.vector(1 + t(z) %*% new.inv %*% z)
new.inv - (new.inv %*% d %*% t(d) %*% new.inv) / as.vector(1 + t(d) %*% new.inv %*% d)
Y.inv

# add new data point
Y = rbind(Y, rnorm(p))
n = n + 1

z = Y[n+1,]
z = sqrt(n+1)*(z - apply(Y, 2, mean)) / n

new.inv = n/(n-1)*Updated
Updated = new.inv - (new.inv %*% z %*% t(z) %*% new.inv) / as.vector(1 + t(z) %*% new.inv %*% z)
solve(cov(Y))

### speed test
set.seed(1)
n = 10000
p = 2000
Y = matrix(rnorm((n+1)*p), n+1, p)

i = 5000
A.cov = cov(Y[1:i,])
A.inv = solve(A.cov)

system.time({z <- Y[i+1,];
z <- sqrt(i+1)*(z - apply(Y[1:(i+1),], 2, mean)) / i;

A.inv <- i/(i-1) * A.inv;
A.inv <- A.inv - (A.inv %*% z %*% t(z) %*% A.inv) / as.vector(1 + t(z) %*% A.inv %*% z)})

system.time(solve(cov(Y[1:(i+1),])))

###########

###########
# normal mixture model
source("~/files/R/mcmc/bayes_functions.R")
rho = 0.7
mu1 = 2
mu2 = 6
sig1 = 0.8
sig2 = 1

gen = function(n){
    r = runif(n)
    l = sum(r < rho)
    out = double(n)
    out[r < rho] = rnorm(l, mu1, sig1)
    out[r >= rho] = rnorm(n-l, mu2, sig2)
    return (out)
    }

calc.post = function(params){
    rho = params[1]
    mu1 = params[2]
    mu2 = params[3]
    sig1 = params[4]
    sig2 = params[5]
    out = 0
    # likelihood
    out = sum(log(rho*dnorm(y, mu1, sig1) +
        (1-rho)*dnorm(y, mu2, sig2)))
    # priors
#   out = out + dbeta(rho, 7, 3, log = TRUE)
#   out = out + dnorm(mu1, 3, 10, log = TRUE)
#   out = out + dnorm(mu2, 3, 10, log = TRUE)
#   out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
#   out = out + dgamma(sig2, 1.5, 0.5, log = TRUE)
    out = out + dbeta(rho, 10, 15, log = TRUE)
    out = out + dnorm(mu1, 24, 1, log = TRUE)
    out = out + dnorm(mu2, 31, 1, log = TRUE)
    out = out + dgamma(sig1, 1.5, 0.5, log = TRUE)
    out = out + dgamma(sig2, 1.75, 0.5, log = TRUE)
    return (out)
    }

y = gen(500)
# mpg data
y=c(209.5/6.171, 296.0/9.664, 240.5/9.867, 221.8/10.231,
        319.3/10.404, 287.5/9.543, 307.3/9.911, 227.9/9.405,
        318.1/9.812, 309.4/9.634, 269.2/9.271, 297.9/10.334,
        163.3/9.913, 300.0/10.12, 355.1/10.384, 286.6/9.611,
        330.6/10.390, 301.0/9.956, 316.2/10.020, 366.9/9.942,
        262.1/10.055, 215.8/10.048, 283.1/8.690, 357.6/9.879,
        242.5/10.121, 188.9/8.485, 311.3/9.870, 225.3/10.236,
        259.8/10.304, 264.8/9.904, 277.3/10.465, 272.5/11.277,
        223.3/10.018)
y = y[-c(13, 20)] # removing the most extreme values
plot(density(y))

nparams = 5
sigs = rep(1, nparams)
nburn = 15000
nmcmc = 25000
upper = c(1, Inf, Inf, Inf, Inf)
lower = c(0, -Inf, -Inf, 0, 0)
window = 500

params = matrix(0, nburn + nmcmc, nparams)
params[1,] = c(0.5, 0, 0, 1, 1)
cand.param = params[1,]
accept = matrix(0, nburn + nmcmc, nparams)

post = calc.post(params[1,])
cand.post = post

for (i in 2:(nburn+nmcmc)){
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
    if (floor(i/window) == i/window)
        cat(i, "/", nburn+nmcmc, "\n")
    }

params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

for (i in 1:nparams){
    plot(params[,i], type='l')
    if (i != nparams)
        readline()
    }

apply(params, 2, mean)
apply(accept, 2, mean)
# it's okay to estimate rho as 1-rho, just make the
# changes with the other parameters and it all
# works out

m = nrow(params)
pred.y = double(m)
mix = rbinom(m, 1, params[,1])
for (i in 1:m){
    if (mix[i] == 1){
        pred.y[i] = rnorm(1, params[i,2], sqrt(params[i,4]))
    } else {
        pred.y[i] = rnorm(1, params[i,3], sqrt(params[i,5]))
        }
    }

plot(density(pred.y))
points(density(y), col='red', type='l')

##########


##########
# random draw? from the function? doubtful
# this is called griddy gibbs i think
# not so useful in high dimensional spaces
f = function(x)
    dnorm(x)
xx = seq(-5, 5, length=100000)

plot(xx, f(xx), type='l')

yy = f(xx)
new.x = sample(xx, 100000, replace = TRUE, prob = yy/sum(yy))
new.y = rnorm(100000)

hist(new.x, freq=FALSE, breaks=100, col='gray')
curve(f(x), add=TRUE, col='red', lwd=2)
points(density(new.y), col='blue', type='l', lwd=2)
points(density(new.x), col='green', type='l', lwd=2)

# some cosine density
f = function(x)
    (cos(x) + 1)/(2*pi)
xx = seq(0, 2*pi, length=100000)

plot(xx, f(xx), type='l')

yy = f(xx)
new.x = sample(xx, 100000, replace = TRUE, prob = yy/sum(yy))

hist(new.x, freq=FALSE, breaks=100, col='gray')
curve(f(x), add=TRUE, col='red', lwd=2)
points(density(new.x), col='green', type='l', lwd=2)
##########

##########
# random draws from a Wishart distribution
# trace function
tr = function(x)
    sum(diag(x))

# multivariate gamma function
# p = 1 is univariate gamma
mgamma = function(a, p)
    pi^(p*(p-1)/4) * prod(gamma(a+(1-1:p)/2))
# log multi gamma
lmgamma = function(a, p)
    p*(p-1)/4*log(pi) + sum(lgamma(a+(1-1:p)/2))

# density of Wishart(V, n) at X
dwishart = function(X, V, df){
    # X is p x p, positive definite (the support)
    # V is p x p, positive definite
    # df > p - 1 is degrees of freedom
    p = nrow(V)
    det(X)^((df-p-1)/2) * exp(-0.5*tr(solve(V) %*% X)) /
        (2^(df*p/2) * det(V)^(df/2) * mgamma(df/2, p))
    }

# log density of Wishart
ldwishart = function(X, V, df){
    # X is p x p, positive definite (the support)
    # V is p x p, positive definite
    # df > p - 1 is degrees of freedom
    p = nrow(V)
    (df-p-1)/2*determinant(X)$modulus[1] - (1/2)*tr(solve(V) %*% X) - (df*p/2)*log(2) -
        (df/2)*determinant(V)$modulus[1] - lmgamma(df/2, p)
    }

V = matrix(c(1, 0.9, 0.9, 1), 2, 2)
X = diag(2)
p = nrow(V)
df = 10

dwishart(X, V, df)
ldwishart(X, V, df)

# Bartlett decomposition (for random draws)
# get one random draw
rwishart = function(V, df){
    p = nrow(V)
    # function to randomize lower triangle of a matrix with independent N(0,1)s
    lower.rand = function(x){
        n = nrow(x)
        z = sequence(n)
        index = cbind(
            row = unlist(lapply(2:n, function(x) x:n), use.names = FALSE),
            col = rep(z[-length(z)], times = rev(tail(z, -1))-1))
        x[index] = rnorm(nrow(index))
        return (x)
        }
    L = t(chol(V))
    A = diag(sqrt(rchisq(p, df-(1:p)+1)))
    A = lower.rand(A)
    return (L %*% A %*% t(A) %*% t(L))
    }

n = 1000
out = rep(list(matrix(0, p, p)), n)
for (i in 1:n)
    out[[i]] = rwishart(V, df)
mean.out = matrix(0, 2, 2)
for (i in 1:n)
    mean.out[1,1] = mean.out[1,1] + out[[i]][1,1]
for (i in 1:n)
    mean.out[1,2] = mean.out[1,2] + out[[i]][1,2]
for (i in 1:n)
    mean.out[2,1] = mean.out[2,1] + out[[i]][2,1]
for (i in 1:n)
    mean.out[2,2] = mean.out[2,2] + out[[i]][2,2]
mean.out = mean.out / n

# theoretical mean: df*V
df * V
mean.out
##########

##########
# rejection algorithm
reject = function(dg, denv, renv, niter){
    # dg is target density
    # denv is envelope density
    # renv a function to produce random draws
    out = matrix(0, niter, 4)
    out[,1] = renv(niter)
    out[,2] = log(runif(niter))
    out[,3] = log(dg(out[,1])) - log(denv(out[,1]))
    out[,4] = out[,2] < out[,3]
    return (out)
    }

a = 2.7
b = 6.3
g = function(x)
    dbeta(x, a, b)
g.star = function(x)
    g(x) / G

# The envelope need not be constrained to the support
# of the target distribution, but you do lose out on
# some efficiency. That cost isn't much if the
# envelope is easier to work with than others.
denv = function(x)
    dnorm(x, 0.2, 0.25)
renv = function(n)
    rnorm(n, 0.2, 0.25)

G = max(dbeta(xx, a, b) / denv(xx))

xx = seq(0, 1, length=1000)
plot(xx, g.star(xx), type='l')
points(xx, denv(xx), type='l', col='green')

out = reject(g.star, denv, renv, 100000)
plot(xx, dbeta(xx, 2.7, 6.3), type='l')
points(density(out[out[,4] == 1, 1]), type='l', col='red')

mean(out[,4])
#########

########
# polynomial contrasts

poly.contrast = function(x, n){
    # default when given x = number of x's
    if (length(x) == 1)
        x = seq(1, x)

    # center x
    x = x - (max(x) - min(x)) / 2

    out = matrix(0, length(n), length(x))

    for (i in 1:length(n)){
        y = x^n[i]
        out[i,] = y - mean(y)
        }

    return (out)
    }

poly.contrast(5, 1:2)
#########
# simple gaussian process example
library(MASS)
dat.x = c(0.2, 0.3, 0.6, 0.8)
dat.y = c(0.5, 0.56, 0.52, 0.40)

l = 0.25
k = function(x){
    out = matrix(0, length(x), length(x))
    for (i in 1:length(x))
        for (j in 1:length(x))
            out[i, j] = exp(-((x[i] - x[j])^2)/(2*l^2))
    return (out)
    }
test.x = seq(0, 1, length=200)

## prior
R = k(test.x)
set.seed(29)
draws = mvrnorm(6, rep(0.5, length(test.x)), R)

#c.lu = apply(draws, 2, quantile, c(0.025, 0.975)) 
#c.mean = apply(draws, 2, mean)

pdf ("~/files/R/figs/gp_prior.pdf", bg = "gray86", width = 9, height = 9)
par(mar=c(5.1, 4.6, 4.1, 2.1))
plot(0, type='n',xlim=c(0,1), ylim=c(-3, 4), xlab = "x", ylab = "y",
    cex.lab = 2.0, main = "Six random draws from a Gaussian process prior", cex.main = 2)
polygon(c(-0.04, -0.04, 1.04, 1.04), c(-3.28, 4.28, 4.28, -3.28), col='white')
for (i in 1:nrow(draws))
    points(test.x, draws[i,], type='l', col='dodgerblue', lwd=2.5)
#points(test.x, c.mean, type='l', lwd=3)
#points(test.x, c.lu[1,], type='l', lwd=3)
#points(test.x, c.lu[2,], type='l', lwd=3)
dev.off()

## posterior
R = k(c(test.x, dat.x))
ind.test = 1:length(test.x)
ind.train = (1+length(test.x)):(length(test.x) + length(dat.x))
post.mu = rep(0.5, length(test.x)) + R[ind.test, ind.train] %*% (solve(R[ind.train, ind.train]) %*%
    (dat.y - rep(0.5, length(dat.y))))
post.sig = R[ind.test, ind.test] - R[ind.test, ind.train] %*%
    solve(R[ind.train, ind.train]) %*% R[ind.train, ind.test]

draws = mvrnorm(50, post.mu, post.sig)

# need a polygon over the plotting region since pdf() will make all the white space
# transparent. the bg option in pdf changes the plotting region and the margins. i
# change everything to gray86 and then use polygon to make the plotting region white
# instead of gray86
pdf("~/files/R/figs/gp_data.pdf", bg = "gray86", width = 9, height = 9)
par(mar=c(5.1, 4.6, 4.1, 2.1))
plot(0, type='n',xlim=c(0,1), ylim=c(0, 1), xlab = "x", ylab = "y",
    cex.lab = 2.0, main = "GP conditioned on noise free observations", cex.main = 2)
polygon(c(-0.04, -0.04, 1.04, 1.04), c(-0.04, 1.04, 1.04, -0.04), col='white')
for (i in 1:nrow(draws))
    points(test.x, draws[i,], type='l', col='dodgerblue', lwd=2.5)
points(dat.x, dat.y, pch=20, lwd=16)
dev.off()

##########
# kernel density estimation

k = function(x, x.obs, h){
    # x:     values at which to estimate density
    # x.obs: random draws
    # h:     bandwidth

    out = double(length(x))
    n = length(x.obs)
    
    # from wikipedia page on kde, when using gaussian kernel
    if (missing(h))
        h = (4/3)^(1/5) * sd(x.obs) * n^(-1/5)

    # standard gaussian kernel
    f = function(x)
        dnorm((x.obs - x)/h, 0, 1)

    # non vectorized form, but not much of a difference in time
    # still slower than just density()
#   for (i in 1:length(x)){
#       out[i] = 1/(n*h)*sum(dnorm((x.obs - x[i])/h, 0, 1))
#       }
    out = sapply(x, function(x) 1/(n*h)*sum(f(x)))

    return (out)
    }

x.draws = rgamma(10000, 1.5, 1.5)
x.draws = rbeta(10000, 1.5, 1.5)
x.draws = rbeta(10000, 0.5, 0.5)

x.seq = seq(min(x.draws), max(x.draws), length=512)
hist(x.draws, freq=FALSE, breaks=50)
points(density(x.draws, from=min(x.draws), to=max(x.draws)),
    type='l', col='red')
points(x.seq, k(x.seq, x.draws), type='l', col='blue')

##########

##########
# visualing importance sampling
g = function(x)
    (3*exp(-0.5*(x-2)^2) + 2*exp(-1/64*(x-8)^2)) / 40
I = function(x)
    dt((x-5)/4, 2)/4

xx = seq(-4, 22, length=100)
w = g(xx)/I(xx) 

par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(xx, g(xx), type='l')
points(xx, I(xx), type='l', col='red')
abline(v=(-5:25), col='gray', lty=2)
abline(v=seq(-5, 25, by=5), col='gray30', lty=2)
legend("topright", box.lty=0, legend=c("Unnormalized target density", "Importance function"),
    col=c("black","red"), lty=1, lwd=2, bg="white"); box()

plot(xx, w, type='l', col='blue'); abline(h=1)
abline(v=(-5:25), col='gray', lty=2)
abline(v=seq(-5, 25, by=5), col='gray30', lty=2)
legend("topright", box.lty=0, legend="g(x) / I(x)", col="blue", lty=1, lwd=2, bg="white"); box()
# Around x=5 we see a valley (below y=1), this means we are going to be getting
# more draws than we should be in that region, so we need to downweight
# those draws. Around x=15 (above y=1), we are getting less draws than we should,
# so upweight those draws. This is akin to the probabilities when using
# bootstrap importance. In the valleys we have low probability of keeping that
# draw in the bootstrap because we will have more than enough draws in that
# region because of the importance function.

# bootstrapping
n = 10000   # don't do too many
theta = 4*rt(n, 2) + 5

g.t = g(theta)
i.t = I(theta)
w.t = (g.t / i.t) / sum(g.t / i.t)

par(mfrow=c(1,1), mar=c(5.1,4.1,2.1,2.1))
plot(theta, w.t, pch=20, xlim=range(xx))

par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(density(w.t))
draws = sample(theta, replace = TRUE, prob = w.t)
plot(density(draws), col='red')
points(xx, g(xx), type='l')

##########

##########
# drawing from a Dirichlet distribution

# n is a scalar, alpha is a vector
rdirichlet = function(n, alpha){
    k = length(alpha)

    # get all the chisq draws
    x = matrix(rchisq(n*k, 2*alpha), n, k, byrow=TRUE)

    # alternate form (the scale parameter doesn't matter, each
    # draw just needs to have the same scale)
    # x = matrix(rgamma(n*k, alpha, 1), n, k, byrow=TRUE)

    # normalize
    x / apply(x, 1, sum)
    }

alpha = c(1.5, 2.8, 7.6, 0.18)
a0 = sum(alpha)
x = rdirichlet(10000, alpha)

# marginals (each is a beta)
for (i in 1:length(alpha)){
    plot(density(x[,i]))
    curve(dbeta(x, alpha[i], a0 - alpha[i]), add=TRUE, col='red')
    if (i < length(alpha))
        readline()
    }

# dots and plots!
plot(x[,c(1,2)], pch=20, cex = 0.1, col = 'black', xlim=c(0,1),
    ylim=c(0,1))
points(x[,c(1,3)], pch=20, cex = 0.1, col = 'blue')
points(x[,c(1,4)], pch=20, cex = 0.1, col = 'red')
points(x[,c(2,3)], pch=20, cex = 0.1, col = 'green')
points(x[,c(2,4)], pch=20, cex = 0.1, col = 'yellow')
points(x[,c(3,4)], pch=20, cex = 0.1, col = 'orange')
##########

##########
# random draws from a discrete
# incomplete

# log of sum identity
# x is a vector for which we compute log(sum(x))
lsum = function(x, method = 2){
    if (method == 0)
        return (log(sum(x)))
    if (method == 1)
        return (log(x[1]) + log(1 + sum(x[2:length(x)] / x[1]))) 
    if (method == 2)
        return (log(x[1]) + log(1 + sum(exp(log(x[2:length(x)]) - log(x[1])))))
    return (0)
    }

x = dpois(0:50, 10)
lx = ppois(0:50, 10, log = TRUE)
y = double(51)
y[1] = lsum(x[1], 0)
for (i in 2:51)
    y[i] = lsum(x[1:i], 2)
cbind(log(cumsum(x)), y, lx)

# n is number of draws
# f is the density (assumed to start at 0)
rdisc = function(n, f, use.log = TRUE){
    x = log(runif(n))
    y = cumsum( )
    }

##########

##########
# priors that have support outside of the parameters in the likelihood
# y ~ Gamma(alpha, beta), alpha > 0, beta > 0
# alpha ~ Normal(0, 1)
# beta ~ Normal(0, 1)

# new dgamma function to return 0 when given invalid parameter values
df = function(y, a, b){
    if (a > 0 && b > 0){
        return (dgamma(y, shape = a, scale = b))
    } else {
        return (0)
        }
    }
autotune = function(accept, target, k){
    k = c((1-1/k)/(cosh(target)-1), (k-1)/(cosh(target-1)-1))
    1+sign(accept-target)*(cosh(accept-target)-1)*
        k[(sign(accept-target)+1)/2+1]
    }
calc.post = function(params){
    a = params[1]
    b = params[2]
    # likelihood
    out = sum(log(df(y, a, b)))
    # priors
    out = out + dnorm(a, 0, 1)
    out = out + dnorm(b, 0, 1)
    return (out)
    }

n = 1000
set.seed(1)
y = rgamma(n, shape = 0.5, scale = 0.2)

# prior predictive
m = 10000
prior = rgamma(m, shape = rnorm(m), scale = rnorm(m))
prior = prior[!is.na(prior)]
hist(y, col='gray', freq = FALSE, breaks=50)
points(density(prior), lwd=3, type='l')

# calc.post(c(0,0))
# calc.post(c(3,2))

# metropolis settings
nburn = 5000
nmcmc = 10000
window = 500

params = matrix(0, nburn+nmcmc, 2)
accept = matrix(0, nburn+nmcmc, 2)
sig = c(1, 1)

params[1,] = c(1,1)
curr.post <- cand.post <- calc.post(params[1,])

for (iter in 2:(nburn+nmcmc)){
    params[iter,] = params[iter-1,]
    cand.vec = params[iter,]
    for (j in 1:2){
        cand.vec[j] = rnorm(1, params[iter-1, j], sig[j])
        cand.post = calc.post(cand.vec)
        if (log(runif(1)) < cand.post - curr.post){
            curr.post = cand.post
            params[iter, j] = cand.vec[j]
            accept[iter, j] = 1
        } else {
            cand.vec[j] = params[iter, j]
            }
        if (floor(iter/window) == iter/window && iter <= nburn)
            sig[j] = sig[j] * autotune(mean(accept[(iter-window+1):iter,j]), 0.25, k = window/50)
        }
    }

# remove the burn-in
params = params[(nburn+1):(nburn+nmcmc),]
accept = accept[(nburn+1):(nburn+nmcmc),]

# some diagnostics
rbind("sig"=sig, "accept"=apply(accept, 2, mean))
apply(params, 2, range)

# trace plots with mean (solid) and 2.5% and 97.5% quantiles (dotted)
plot(params[,1], type='l')
abline(h=mean(params[,1]))
abline(h=quantile(params[,1], c(0.025, 0.975)), lty=2)

plot(params[,2], type='l')
abline(h=mean(params[,2]))
abline(h=quantile(params[,2], c(0.025, 0.975)), lty=2)

# bivariate
plot(params, pch=20)

# posterior predictive
preds = double(nmcmc)
for (i in 1:nmcmc)
    preds[i] = rgamma(1, shape = params[i, 1], scale = params[i, 2])

hist(y, col='gray', freq = FALSE, breaks=50)
points(density(y), col='blue', lty=2, type='l', lwd=3)
points(density(preds), type='l', lwd=3)
##########

##########
# entropy
# normal dist
mu = 5
sig = 7.2
x = rnorm(100000, mu, sig)
mean(-log(dnorm(x, mu, sig))) # sample entropy
0.5*log(2*pi*exp(1)*sig^2) # theoretical

# cauchy
loc = 3
sca = 2.7
x = rcauchy(100000, loc, sca)
mean(-log(dcauchy(x, loc, sca))) # sample entropy
log(sca) + log(4*pi) # theoretical

# beta
a = 0.5
b = 1.2
x = rbeta(100000, a, b)
mean(-log(dbeta(x, a, b))) # sample
lbeta(a, b) - (a-1)*digamma(a) -(b-1)*digamma(b) +
    (a+b-2)*digamma(a+b) # theoretical

# uniform
a = -2
b = 7
x = runif(100000, a, b)
mean(-log(dunif(x, a, b))) # sample
# the uniform entropy is constant for all x and
# is equal to -log(1/(b-a))
# Unif(0,1) is the only continuous probability
# distribution with entropy 0. (i think anyway)
##########

##########
# time difference between runif and rnorm

# initialize
system.time(runif(1))
system.time(rnorm(1))

m = 10
nmcmc = 20000
nparams = floor(seq(1, 1000, length = m))
niter = nmcmc*nparams
tunif = double(m)
tnorm = double(m)

# n is for number of repititions to be averaged over
n = 10

for (i in 1:m){
    cat(nparams[i], "/", max(nparams), "\n")
    for (j in 1:n){
        tunif[i] = tunif[i] + system.time(runif(niter[i]))[3]/n
        tnorm[i] = tnorm[i] + system.time(rnorm(niter[i]))[3]/n
        }
    }


plot(nparams, tunif, type='l')
points(nparams, tnorm, type='l', col='red')
plot(nparams, tnorm, type='l')
points(nparams, tunif, type='l', col='red')

plot(tnorm/tunif, type='l')
##########

##########
# removing a variable is supposed to reduce variation in the betas?
n = 200
p = 1
q = 11

x = matrix(rnorm(n*q), n, q)
x2 = x[,-1]
y = matrix(rnorm(n*p), n, p)

B = solve(t(x) %*% x) %*% t(x) %*% y
s2 = as.vector(1/(n - q) * t(y - x %*% B) %*% (y - x %*% B))
zz1 = diag(s2*solve(t(x) %*% x))

q = 10
B = solve(t(x2) %*% x2) %*% t(x2) %*% y
s2 = as.vector(1/(n - q) * t(y - x2 %*% B) %*% (y - x2 %*% B))
zz2 = c(Inf, diag(s2*solve(t(x2) %*% x2)))

cbind(zz1, zz2, zz1 < zz2)
##########

##########
system.time(for (i in 1:5000000) a = 5)
system.time(for (i in 1:5000000) a <- 5)

niter = 5000000
nloop = 250

# for any preprocessing problems
nburn = 5

t1 = double(nloop+nburn)
t2 = double(nloop+nburn)

for (i in 1:(nburn+nloop)){
    t1[i] = system.time(for (j in 1:niter) a = 5)[3]
    t2[i] = system.time(for (j in 1:niter) a <- 5)[3]
    cat("\rIteration:", i)
    }
t1 = t1[-(1:nburn)]
t2 = t2[-(1:nburn)]
plot(t1, type='l')
points(t2, type='l', col='red')
c(mean(t1), mean(t2))


t3 = double(nloop+nburn)
t4 = double(nloop+nburn)

for (i in 1:(nburn+nloop)){
    t4[i] = system.time(for (j in 1:niter) a <- 5)[3]
    t3[i] = system.time(for (j in 1:niter) a = 5)[3]
    cat("\rIteration:", i)
    }
t3 = t3[-(1:nburn)]
t4 = t4[-(1:nburn)]
plot(t3, type='l')
points(t4, type='l', col='red')
c(mean(t3), mean(t4))

matrix(c(mean(t1), mean(t3), mean(t2), mean(t4)), 2, 2)

plot(t1, type='l', col='blue', lty=2)
points(t3, type='l', col='blue')

points(t2, type='l', col='red', lty=2)
points(t4, type='l', col='red')
##########

##########
# which.max of a matrix

mat.max = function(x, method = 1){
    index = double(2)
    if (method == 1){
        n = nrow(x)
        w = which.max(x) - 1
        index[1] = 1 + w %% n
        index[2] = 1 + w %/% n
        }
    return(index)
    }
##########

##########
# cdf of multivariate normal
library(MASS)

mu = c(2, 5)
sigma = matrix(c(1, 0.8, 0.8, 2), 2, 2)

n = 5000
x = mvrnorm(n, mu, sigma)

xbar = apply(x, 2, mean)
S = var(x)

plot(x, pch=20)

dis = x - matrix(mu, n, length(mu), byrow = TRUE)
inv = solve(S)

mal = -0.5*diag(dis %*% (inv %*% t(dis)))

# uniform?
plot(seq(0, 1, length=n), sort(exp(mal)), pch=20, cex = 0.5)
abline(0, 1)

pmin = which.min(mal)
pmax = which.max(mal)

x[pmin,]
x[pmax,]
plot(x, pch=20)
points(x[pmin, 1], x[pmin, 2], col='red', pch=20, cex=2)
points(x[pmax, 1], x[pmax, 2], col='green', pch=20, cex=2)
##########

##########
# comparison of mean, median, and mode
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
n = 100
niter = 10000

mn = double(niter)
md = double(niter)
#mo = double(niter)
for (i in 1:niter){
    cat("\rIteration:",i,"/",niter)
    x = rnorm(n)
#   x = rbeta(n, 3, 1)
#   x = rlnorm(n, 0, 1)
#   x = rgamma(n, 0.5, 0.1)
#   x = rcauchy(n)
    md[i] = median(x)
    mn[i] = mean(x)
    mo[i] = calc.mode(density(x, n = 2^12))
    if (i == niter)
        cat("\n")
    }

plot(density(x))
plot(density(mn, n=2^20), xlim=range(c(mn, md, mo)), ylim=c(0, max(density(mn)$y, density(md)$y,
    density(mo)$y)))
lines(density(md, n =2^20), col='red')
lines(density(mo, n =2^20), col='blue')

var(mn)
var(md)
var(mo)
# for non-heavy tail distributions, it seems var(mean) < var(median) < var(mode)
# for heavy-tail distributions, it seems var(mode) < var(median) < var(mean)

# for cauchy, it seems var(median) < var(mode) < var(mean), but this could be
# due to the way i'm calculating the mode
############

############
# factor analysis

set.seed(1)
p = 7
m = 2
lambda = matrix(rnorm(p*m), p, m)
psi = diag(sort(runif(p)))

sigma = lambda %*% t(lambda) + psi

# get random draws
n = 10000
library(MASS)
y = mvrnorm(n, rep(0, p), sigma)

R = cor(y)
eig = eigen(R)
eig$values / sum(eig$values)
cumsum(eig$values) / sum(eig$values)

# take the first k eigenvectors
k = 2

L = eig$vectors[,1:k] * sqrt(eig$values[1:k])

P = diag(1 - diag(L %*% t(L)))

(L %*% t(L) + P) - R
# these seem similar enough, but P and L are fairly
# different from psi and lambda

##########

##########
iter = 10000
n = 10
mins = double(iter)
for (i in 1:iter)
    mins[i] = min(runif(n, 0, 1))
plot(density(mins))
curve(dbeta(x, 2, 5), col='red')
##########

##########
make.mat = function(a, b, n)
    matrix(a, n, n) + diag(b, n)

invert.mat = function(a, b, n)
    diag(1/b, n) - a/(b*(b+a*n))
#########

#########
set.seed(1)
n = 20000

f = function(x, a, c, d)
    d^c / beta(a, c) * x^(a-1) / (x + d)^(a+c)

g = function(x, a, c, d)
    (x^a) / (x+d)^(a+c)
h = function(x)
    g(x, a, c, d)

plot(xx, h(xx), type='l')

int(h, c(0, 100000000000), 100000)


a = 1.5
c = 1.0
d = 1
xx = seq(0, 1000, length=5000)

y = rgamma(n, c, d)
z = rgamma(n, a, y)
mean(z); var(z)

plot(z, y, pch=20)

plot(density(z, n=10000), xlim=c(0, 50))
lines(xx, f(xx, a, c, d), col='red')
# mean:
1/(beta(a, c)) * 1/(a+c-2) * 1/(d^(a-2)) * 1/(a+c-1)
mean(z)

cor(cbind(z, y))

plot(density(z))
curve(dgamma(x, mean(z)^2/var(z), mean(z)/var(z)), add = TRUE, col='red')

plot(density(z))
lines(density(rgamma(n, mean(z)^2/var(z), mean(z)/var(z))), col='red')

ks.test(z, j)
ks.test(y, "pgamma", 1, 1)


plot(density(y))
curve(dgamma(x, mean(y)^2/var(y), mean(y)/var(y)), add = TRUE, col='red')


ss = seq(0, 1, length = 1000)
ll = seq(1, 100000, length = 1000)

g = function(x, a, c, d)
    -(x+d)^(-(a+c)+2) / (a+c-2)

a = 1.5
c = 1.0
d = 2
par(mfrow=c(2,1))
plot(ss, g(ss, a, c, d), pch=20)
plot(ll, g(ll, a, c, d), pch=20)

 -1/((a+c)-2) / (d^(a+c-2))

library(hypergeo)
hypergeo(a+1, a+c, a+2, -100000000000/d)
#######
