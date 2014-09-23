raw.dat = read.table("~/files/R/666/data/oliver3b.txt", head=TRUE)
raw.dat = as.matrix(raw.dat)
dat = raw.dat
#dat = dat[-c(40, 3, 165, 132, 176),]
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
                    z = kde2d(x[,j], x[,i])
                    plot(NA, xlim = range(z$x), ylim = range(z$y), axes = FALSE)
                    .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
                        col=gray(seq(0.0, 1.0, length=20)))
#                   .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20),
#                       col=two.colors(20, start="dodgerblue"))
#                   points(dat[,4], dat[,1], pch=20)
                } else {
                    plot(x[,j], x[,i], pch=20, axes = FALSE)
#                   points(x[40,j], x[40,i], pch=20, col='red')
#                   points(x[191,j], x[1,i], pch=20, col='blue')
                    }
                }
            box()
            }
        }
    }

#pairs(dat)
pdf("figs/pairs_raw.pdf", 8, 8)
mw.pairs(dat)
dev.off()


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
plot(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(E),
    pch = 20)
abline(0,1)
order(D, decreasing = TRUE)

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
                width[j] = width[j] * temperature
                }
            }
        if (floor(i/window) == i/window && print)
            cat(i,"/",niter, "\n")
        }
    return (lam)
    }

mult.lam = max.loop(dat, print = TRUE)
pdf("figs/pairs_trans.pdf", 8, 8)
mw.pairs(lam.func(dat, mult.lam))
dev.off()

box.dat = lam.func(dat, mult.lam)
box.x.bar = apply(box.dat, 2, mean)
box.S = var(box.dat)
box.S.inv = solve(box.S)

# qq plotting
D = double(n)
for (i in 1:length(D))
    D[i] = as.vector(t(box.dat[i,]-box.x.bar) %*% box.S.inv %*%
        (box.dat[i,]-box.x.bar))

# remove 3 highest D's (possible outliers)
dat.out = raw.dat[-order(D, decreasing = TRUE)[1:3],]
n.out = nrow(dat.out)
p.out = ncol(dat.out)
mult.lam.out = max.loop(dat.out, print = TRUE)
box.dat.out = lam.func(dat.out, mult.lam.out)
box.x.bar.out = apply(box.dat.out, 2, mean)
box.S.out = var(box.dat.out)
box.S.inv.out = solve(box.S.out)
E = double(n.out)
for (i in 1:length(E))
    E[i] = as.vector(t(box.dat.out[i,]-box.x.bar.out) %*% box.S.inv.out %*%
        (box.dat.out[i,]-box.x.bar.out))

pdf("figs/qq_trans.pdf", 8, 8)
u = n/(n-1)^2*sort(D)
v = qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2)
plot(u, v, pch = 20, xlab = "u", ylab = "v")
abline(0,1)
dev.off()
pdf("figs/qq_outlier.pdf", 8, 8)
u = n.out/(n.out-1)^2*sort(E)
v = qbeta((1:n.out - 0.5)/n.out, p.out/2, (n.out - p.out - 1)/2)
plot(u, v, pch = 20, xlab = "u", ylab = "v")
abline(0,1)
dev.off()
