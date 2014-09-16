dat = read.table("~/files/R/666/data/oliver3b.txt", head=TRUE)
dat = as.matrix(dat)
n = nrow(dat)
p = ncol(dat)

mw.pairs = function(x){
    require(MASS)
    par(mfrow=rep(ncol(x), 2), mar=rep(0, 4))
    for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
            if (i == j){
                hist(x[,i], axes = FALSE, main = "", col = "gray")
                #plot(density(x[,i]), axes = FALSE, main = "")
                legend("topright", legend = i, box.lty = 0, cex = 1.5)
            } else {
                if (i > j){
                    z = kde2d(x[,i], x[,j])
                    plot(NA, xlim = range(z$x), ylim = range(z$y), axes = FALSE)
                    .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
                        col=gray(seq(0.0, 1.0, length=20)))
#                   points(dat[,4], dat[,1], pch=20)
                } else {
                    plot(x[,j], x[,i], pch=20, axes = FALSE)
                    points(x[40,j], x[40,i], pch=20, col='red')
                    }
                }
            box()
            }
        }
    }

mw.pairs(dat)

x.bar = apply(dat, 2, mean)
S = var(dat)
S.inv = solve(S)

# qq plotting
D = double(n)
for (i in 1:length(D))
    D[i] = as.vector(t(dat[i,]-x.bar) %*% S.inv %*% (dat[i,]-x.bar))

plot(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D),
    pch = 20)
abline(0,1)

### box-cox
# do the transformation
lam.func = function(x, lam){
    out = matrix(0, nrow(x), ncol(x))
    for (i in 1:ncol(x)){
        if (lam[i] == 0){
            out[,i] = log(x[,i])
        } else {
            out[,i] = (x[,i]^lam[i] - 1)/lam[i]
            }
        }
    return (out)
    }
# multivariate function to maximize
max.func = function(x, lam){
    y = lam.func(x, lam)
    -n/2*determinant((n-1)/n*cov(y))$modulus[1] +
        sum(apply(log(x), 2, sum) * (lam-1))
    }
# R function to do the random walk maximization
max.loop = function(x, lambda, niter = 1000, temperature = 0.99, width, print = FALSE){
    x = as.matrix(x)
    n = nrow(x)
    p = ncol(x)
    if (missing(lambda)){
        lam = runif(p, -2 ,2)
    } else {
        lam = lambda
        }
    if (missing(width))
        width = rep(1, p)
    temp = temperature
    current = max.func(x, lam)
    window = round(niter/20)
    for (i in 1:niter){
        cand.lam = lam
        for (j in 1:p){
            cand.lam[j] = runif(1, lam[j]-width[j], lam[j]+width[j])
            candidate = max.func(x, cand.lam)
            if (candidate > current){
                current = candidate
                lam[j] = cand.lam[j]
            } else {
                width[j] = width[j] * temp
                }
            }
        if (floor(i/window) == i/window && print)
            cat(i,"/",niter, "\n")
        }
    return (lam)
    }

start.lam = double(p)
for (i in 1:p)
    start.lam[i] = max.loop(dat[,i], print = TRUE)
mult.lam = max.loop(dat, lam = start.lam, print = TRUE)

max.func(dat, start.lam)
max.func(dat, mult.lam)

mw.pairs(dat)
mw.pairs(lam.func(dat, start.lam))
mw.pairs(lam.func(dat, mult.lam))

plot(dat[,4], dat[,1], pch=20)
#identify(dat[,4], dat[,1])
# observation 40 looks to be an outlier

plot(dat[,1], type='l')

new.dat = lam.func(dat, mult.lam)
x.bar = apply(new.dat, 2, mean)
S = var(new.dat)
S.inv = solve(S)

# qq plotting
D = double(n)
for (i in 1:length(D))
    D[i] = as.vector(t(new.dat[i,]-x.bar) %*% S.inv %*%
        (new.dat[i,]-x.bar))

plot(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D),
    pch = 20)
identify(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D))
abline(0,1)

### kurtosis and skewness
