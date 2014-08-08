# Find a near-maximin, near-orthogonal subset of design points.
# Motivated by a set of computer experiments containing a very large
# number of runs. We wish to find a subset of the runs which have
# properties similar to Latin hypercubes having maximin and
# orthogonality properties since use of the full dataset will be very
# costly in obtaining estimates for a Gaussian process model where
# matrix inversion is required.

# todo: run diagnostics to see how well the chosen subset
#       meets desired properties (maximin, orthogonal factors, uniform
#       spread over marginals, latin hypercube, no repitition of
#       factor levels, and others?)
msub = function(x, m, w1 = 1/3, w2 = 1/3, pb = FALSE){
    # x: the design matrix, number of points by number of dimensions
    # m: <= nrow(x), the size of the subset to find
    # note: there seems to be a tendency to choose around the
    #       edges before filling in the middle, but a high enough
    #       m this shouldn't be too much of an issue
    # perhaps other distance metrics may be used. i tried minkowski
    # distance to compute the 'norms' variable below and all choices
    # of p gave the same result, but maybe not so in the j for loop

    if (pb)
        source("~/files/R/sandbox/pb_mickey.R")

    euclid = function(x)
        sqrt(sum(x^2))
    n = nrow(x)
    d = ncol(x)
    xmin = apply(x, 2, min)
    xrange = apply(x, 2, function(x) diff(range(x)))
    # scale x to [-0.5, 0.5]^d (really just needs to be centered
    # around 0)
    for (i in 1:d)
        x[,i] = ((x[,i] - xmin[i]) / xrange[i]) - 0.5

    # a default value for m
    if (missing(m))
        m = floor(n^(1/d))
    out = double(m)

    # compute distance from the origin
    norms = apply(x, 1, euclid)
    # take the first point in the subset to be that which is furthest
    # from the origin. the point will be on the edge and more likely
    # to be in a "corner" than, say, on the wider edge of an ellipse
    out[1] = which.max(norms)

    # initialize the dists matrix. these values are the distances
    # from each point in the design and the ones currently selected
    # to be in the subset. the column of Inf's is included so the
    # for loop will always consider dists to be a matrix on not
    # automatically change it to a vector (so apply always works)
    t.vec = double(n)
    v.vec = Inf+double(n)
    r.vec = 1 + double(n)
    
    if (pb)
        pb.Init()

    for (j in 2:m){
        if (pb)
            pb.Update(j-1, m-1)
        if (w1 != 0){
            t.vec = apply(x - matrix(x[out[j-1],], n, d, byrow=TRUE),
                1, euclid)
            t.vec = t.vec / max(t.vec)
            v.vec = apply(cbind(v.vec, t.vec), 1, min)
            }

        if (w2 != 0){
            if (j > 2 && d > 1){
                r.vec = 1 + double(n)
                for (k in (1:n)[-out])
                    r.vec[k] = max(abs(cor(x[c(k, out),])) - diag(d))
                }
            }
        if (w1 + w2 != 1){
            z.vec = double(n)
            s.vec = Inf+double(n)
            for (l in 1:d){
                for (k in 1:n){
                    z.vec[k] = min(abs(x[out,l] - x[k, l]))
                    }
                z.vec = z.vec / max(z.vec)
                s.vec = apply(cbind(s.vec, z.vec), 1, min)
                }
        } else {
            s.vec = double(n)
            }
        out[j] = which.max(w1*v.vec + w2*(1-r.vec) + (1-w1-w2)*s.vec)
        }
    # returns the index of points in the subset
#   if (length(unique(out)) < length(out))
#       print("Duplicates selected in the subset. Removing.")
#   return (unique(out))
    return (out)
    }


# the correlation is good (low), but the spread is still not
# quite right. i need to put weights back on the factors with
# similar levels of previously chosen subset points. but where?
# perhaps on ( w*vec + (1-w)*(1-r.vec) ) ?





# example
n = 10000
d = 10
m = 500
x = matrix(runif(n*d), n, d)
# x = matrix(rbeta(n*d, 1, 0.25), n, d)
# x = matrix(mvrnorm(n,c(0.5,0.5),matrix(c(1,0.8,0.8,1),2,2)),n,d)

# takes a bit over a minute (n = 10000, d = 10, m = 500)
system.time(out1 <- msub(x1, m1, w1 = 0.5, w2 = 0.5, pb = TRUE))


cols = rep("gray30", n)
cols[out] = "red"
lwds = rep(1, n)
lwds[out] = 3
pairs(x, pch=20, col=cols, lwd=lwds)



pairs(x[out,], pch=20, col='red')
pairs(x[out2,], pch=20)
pairs(x[out3,], pch=20)

rbind(out1, out2, out3)

plot3d(x)
points3d(x[out,], col='red', size=5)

pairs(x, pch=20)
for (i in 2:m){
    readline()
    pairs(x[out[1:i],], col='red', xlim=c(0,1), ylim=c(0,1))
    }

# for d = 2
plot(x, pch=20)
points(matrix(x[out,], ncol=2),col='red',cex=1.5, lwd=2)

# see the choice of subset points one at a time
for (i in 1:m){
    points(matrix(x[out[i],], ncol=2),col='red',cex=1.5, lwd=2)
    readline()
    }

# diagnostics should be run on x[out,], it being the subset
# out (from out <- msub(x, m)) is the index of points

###
library(foreach)
library(doMC)
registerDoMC(2)

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

howmany = 11
m = 11
x = cbind(rep(seq(0, 1, length=howmany), each=howmany),
    rep(seq(0, 1, length=howmany), times=howmany))
n = nrow(x)
d = ncol(x)

where = cross(seq(0, 1, length=101), seq(0, 1, length=101))
for (i in 1:nrow(where)){
    if (sum(where[i,]) > 1){
        where[i,] = Inf
        }
    }
where = where[!duplicated(where),]
where = where[-which.max(where),]

doit = function(i)
    msub(x, m, w1 = where[i,1], w2 = where[i,2])

out = foreach (i = 1:nrow(where), .combine = rbind) %dopar% doit(i)

that = function(i)
    dist(where) 
testing = foreach (i = 1:100, .combine = rbind) %dopar% that(i)

p = 2
func.phi.row = function(DIST)
    apply(1/(DIST^p+diag(Inf, m)), 1, sum)^(1/p)
func.rho = function(COR)
    (sum(COR^2)-d)/(d*(d-1))
func.phi = function(DIST)
    (sum(func.phi.row(DIST)^p)/2)^(1/p)

rho = foreach (i = 1:nrow(where), .combine = rbind) %dopar% 
    func.rho(cor(x[out[i,],]))
phi = foreach (i = 1:nrow(where), .combine = rbind) %dopar% 
    func.phi(as.matrix(dist(x[out[i,],], method="manhattan")))

i = which.min(phi)
for (i in which(rho == min(rho))){
    readline()
    plot(x, pch=20, main=paste0("Dist: ", where[i,1], " -- Corr: ",
        where[i,2], " -- Comp: ", 1-where[i,1]-where[i,2]))
    points(x[out[i,],], cex=1.5, lwd=2, col='red')
    }


### Regression example

# true model
gen.y = function(x1, x2, x3)
    1.5 + 0.2*x1 + 2*x2 + 1.3*x1*x2 + 0.01*x3 + rnorm(length(x1), 0, 0.1)

n = 10000
ncr = floor(n*0.5)
x1 = runif(n, -5, 5)
x2 = runif(n, 2, 7)
x3 = runif(n, -1, 3)

y = gen.y(x1, x2, x3)

X = data.frame(x1, x2, x3)

ind.train = sample(n, ncr)
ind.test = (1:n)[-ind.train]
mod = lm(y ~ x1 + x2 + x1*x2 + x3, data = X)

y.pred = predict(mod, newdata = X)

mean((y.pred - y[ind.test])^2)

m = floor(n*0.05)
mcr = floor(m*0.5)
x = as.matrix(X)
system.time(sub.ind <- sort(msub(x, m, w1 = 1.0, w2 = 0.0, pb = TRUE)))

sub.train = sort(sample(sub.ind, mcr))
sub.test = sub.ind[which(table(c(sub.ind, sub.train)) == 1)]

sub.mod = lm(y ~ x1 + x2 + x1*x2 + x3, data = X, subset = sub.ind)
summary(sub.mod)

sub.train2 = sort(sample(sub.ind2, mcr))
sub.test2 = sub.ind[which(table(c(sub.ind2, sub.train2)) == 1)]
sub.mod2 = lm(y ~ x1 + x2 + x1*x2 + x3, data = X, subset = sub.ind2)

sub.ind3 = sample(n, m)
sub.mod3 = lm(y ~ x1 + x2 + x1*x2 + x3, data = X, subset = sub.ind3)

th = cbind(summary(mod)[[4]][,2], summary(sub.mod)[[4]][,2],
    summary(sub.mod2)[[4]][,2],
    summary(sub.mod3)[[4]][,2])

y.sub1 = predict(sub.mod, newdata = X[sub.ind,])
y.sub2 = predict(sub.mod, newdata = X[sub.ind2,])
y.sub3 = predict(sub.mod, newdata = X[sub.ind3,])
mean((y.pred - y)^2)
mean((y.sub1 - y[sub.ind])^2)
mean((y.sub2 - y[sub.ind2])^2)
mean((y.sub3 - y[sub.ind3])^2)

plot(x[sub.ind, 1], y[sub.ind], pch=20)
pairs(x[sub.ind,])

library(rgl)
plot3d(x[,1], x[,2], y)
points3d(x[sub.ind,1], x[sub.ind,2], y[sub.ind], col='red', size=8)


system.time(for (i in 1:10000){
    test = max(abs(cor(x[c(k, out),])) - diag(d))
    })

system.time(for (i in 1:100){
    test = 0
    for (j in 1:(d-1))
        for (k in (j+1):d)
            test = max(test, cor(x[c(k, out),j], x[c(k, out),k]))
    })

mod = lm(y ~ x1 + x2 + x1*x2 + x3, data = X)

x = as.matrix(X)
m = 100
sub.ind = msub(x, m, w1 = 0.5, w2 = 0.5, pb = TRUE)

sub.mod = lm(y ~ x1 + x2 + x1*x2 + x3, data = X, subset = sub.ind)

pred = predict(mod)
sub.pred = predict(sub.mod, newdata = X)

mean((y - pred)^2)
mean((y - sub.pred)^2)


### cmaq example
cmaq.dat = read.table("~/files/R/data/cmaq_cmaq_data.csv", header=TRUE,
    sep = ',')
cmaq.dat = cmaq.dat[,-3]

cmaq.rand.ind = sample(nrow(cmaq.dat), 500)
cmaq.msub.ind = msub(as.matrix(cmaq.dat), 500, w1 = 1.0, w2 = 0.0, T)
nu = 1
s2.start.val = 35
phi.start.val = 1.5

cmaq.mod = likfit(as.geodata(cmaq.dat[cmaq.rand.ind,], data.col = 3,
    coords.col = c(1,2)), cov.model = "matern", kappa = nu,
    fix.kappa = TRUE, ini.cov.pars = c(s2.start.val, phi.start.val),
    trend = "cte")

cmaq.mod2 = likfit(as.geodata(cmaq.dat[cmaq.msub.ind,], data.col = 3,
    coords.col = c(1,2)), cov.model = "matern", kappa = nu,
    fix.kappa = TRUE, ini.cov.pars = c(s2.start.val, phi.start.val),
    trend = "cte")

