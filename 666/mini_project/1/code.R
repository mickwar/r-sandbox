dat = read.table("~/files/R/666/data/oliver3b.txt", head=TRUE)
dat = as.matrix(dat)
n = nrow(dat)
p = ncol(dat)

mw.pairs = function(x){
    par(mfrow=rep(ncol(x), 2), mar=rep(0, 4))
    for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
            if (i == j){
                hist(x[,i], axes = FALSE, main = "", col = "gray")
                #plot(density(x[,i]), axes = FALSE, main = "")
                legend("topright", legend = i, box.lty = 0, cex = 1.5)
            } else {
                plot(x[,j], x[,i], pch=20, axes = FALSE)
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

# do it
lam = runif(p, -2, 2)
width = rep(2, p)
temp = 0.99
niter = 10000

current = max.func(dat, lam)
for (i in 1:niter){
    cand.lam = lam
    for (j in 1:p){
        cand.lam[j] = runif(1, lam[j]-width[j], lam[j]+width[j])
        candidate = max.func(dat, cand.lam)
        if (candidate > current){
            current = candidate
            lam[j] = cand.lam[j]
        } else {
            width[j] = width[j] * temp
            }
        }
    if (floor(i/500) == i/500)
        cat(i,"/",niter, "\n")
    }

new.dat = lam.func(dat, lam)
mw.pairs(new.dat)

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
abline(0,1)
