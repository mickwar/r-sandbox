library(rgl)
library(MASS)
library(Matrix)
source("~/files/R/experimental_design/omlhd.R")
source("~/files/R/sandbox/pb_linux.R")
source("./construct_functions.R")
source("./mcmc_functions.R")
source("~/files/R/mcmc/bayes_functions.R")

# ranjan example
sim.f = function(x){
    x1 = x[1]
    x2 = x[2]
    t = x[3]
    (30+5*x1*sin(5*x1))*(6*t+1+exp(-5*x2))
    }

var.y = 0.5^2
eps.f = function(x)
    rnorm(length(x), 0, sqrt(var.y))

disc.f = function(x){
    x1 = x[1]
    x2 = x[2]
    -50*exp(-0.2*x1-0.1*x2)
    }

# for new data
n = 15
field.X = (maximin(n, 2, w=0.1, t=0.1, alpha=10,niter=10000 , pb=TRUE) - 1)/n
field.Y = apply(field.X, 1, function(x) sim.f(c(x, 0.5)) + eps.f(1) + disc.f(x))

m = 200
sim.X = (maximin(m, 3, w=0.5, t=0.1, alpha=20, niter=10000, pb=T) - 1)/m
sim.Y = apply(sim.X, 1, sim.f)

label = c("lamW", "lamV", "lamY", "lamEta", "rhoW1", "rhoW2",
    "rhoW3", "rhoV1", "rhoV2", "theta")

# load data
all = read.table("./out_all_draws.txt")
field = read.table("./out_field.txt")
field.X = field[,1:2]
field.Y = field[,3]
n = nrow(field.X)
px = ncol(field.X)
simul = read.table("./out_simul.txt")
sim.X = simul[,1:3]
sim.Y = simul[,4]
m = nrow(sim.X)
pt = ncol(sim.X) - px
lam.w = as.matrix(all[,1])
lam.v = as.matrix(all[,2])
lam.eta = as.matrix(all[,3])
lam.y = as.matrix(all[,4])
rho.w = as.matrix(all[,5:7])
rho.v = as.matrix(all[,8:9])
theta = as.matrix(all[,10])

# higdon model
# 2.2.1 design of simulator runs
# inputs are already in [0,1]^3
px = 2
pt = 1

# 2.2.2 simulator model

# set up principal components part
n.eta = 1
ysim.mean = mean(sim.Y)
ysim.sd = sd(sim.Y)
Ksi = t((sim.Y - ysim.mean)/ysim.sd)
sim.yStd = t(Ksi)
K.eta = as.matrix(svd(Ksi)$u)

p.eta = 1 # not to exceed n.eta
K.eta = as.matrix(K.eta[,1:p.eta])

# set up K matrix
K = matrix(0, m*n.eta, m*p.eta)
for (i in 1:p.eta)
    K[1:(m*n.eta),1:(m*n.eta)] = kronecker(diag(m), K.eta)


# 2.2.3 discrepancy model
Fgroups = 1
field.yStd = (field.Y - ysim.mean)/ysim.sd
p.delta = 1
D = diag(15)



# D_i,jk = d_k(x_ij, y_ij)
# Dgrid = seq(min(sim.Y), max(sim.Y), length=p.delta)
# Dwidth = Dgrid[2] - Dgrid[1]
# Dgrid = t(c(Dgrid[1] - Dwidth, Dgrid, Dgrid[p.delta] + Dwidth))
# p.delta = length(Dgrid)     # p.delta + 2
# Dsim = matrix(0, ncol(Ksi), p.delta)
# Dobs = matrix(0, length(field.yStd), p.delta)
# for (i in 1:p.delta){
#     Dsim[, i] = dnorm(t(Ksi), Dgrid[i], Dwidth)
#     Dobs[, i] = dnorm(field.yStd, Dgrid[i], Dwidth)
#     }

# matrices for likelihood
#W = R.x(t(field.yStd), double(length(field.yStd))) # assuming no correlation between obsverations
W = diag(var.y, 15)
P.D = matrix(0, n*p.delta, n*p.delta)
for (i in 1:p.delta){
    for (j in 1:n){
        P.D[j+n*(i-1), (j-1)*p.delta+i] = 1
        }
    }
P.K = matrix(0, n*p.eta, n*p.eta)
for (i in 1:p.eta){
    for (j in 1:n){
        P.K[j+n*(i-1), (j-1)*p.eta+i] = 1
        }
    }
P = diag2(t(P.D), t(P.K))
B = cbind(D, diag(n))  %*% P





#joint.out = rbind(as.matrix(field.yStd), t(Ksi)) # for regular
# using result 1
ridge = 1e-6
z.vec = solve(t(B) %*% W %*% B + diag(ridge, 2*n))%*%t(B)%*%W%*%field.yStd
z.vec = rbind(z.vec, solve(t(K)%*%K)%*%t(K)%*%sim.yStd)


# MCMC

# ends at line 384
niter = 1000 # 249 ^j

# initialize paramters, first rows are starting values
lam.w = matrix(NA, niter, p.eta)
lam.w[1,] = 0.4
acc.lam.w = double(p.eta)

rho.w = matrix(NA, niter, px+pt)
rho.w[1,] = c(0.9,0.9,0.9) 
acc.rho.w = double(px+pt)

lam.eta = matrix(NA, niter, 1)
lam.eta[1,] = 70
acc.lam.eta = double(1)

lam.v = matrix(NA, niter, Fgroups)
lam.v[1,] = 2
acc.lam.v = double(Fgroups)

rho.v = matrix(NA, niter, px)
rho.v[1,] = c(0.6, 0.6)
acc.rho.v = double(px)

lam.y = matrix(NA, niter, 1)
lam.y[1,] = 6
acc.lam.y = double(1)

theta = matrix(NA, niter, pt) # uniform prior
theta[1,] = 0.3
acc.theta = double(pt)

params = list("lamw"=c(1.1, 1.1), "rhow"=c(1.1, 1.1), "lameta"=c(1.1, 1.1),
    "lamv"=c(1.1, 1.1), "rhov"=c(1.1, 1.1), "lamy"=c(1.1, 1.1))
# using result 1
params$lameta = params$lameta + c(m*(n.eta-p.eta)/2,
    0.5*t(sim.yStd)%*%(diag(m)-K%*%solve(t(K)%*%K)%*%t(K))%*%sim.yStd)
params$lamy = params$lameta + c((n-rankMatrix(B))/2,
    0.5*t(field.yStd)%*%(W-W%*%B%*%solve(t(B)%*%W%*%B + diag(ridge, n*2))%*%t(B)%*%W)%*%field.yStd)

#sigmas = list("lamw"=0.4, "rhow"=c(0.13, 0.1, 0.1), "lameta"=15,
#    "lamv"=3.5, "rhov"=c(0.8, 0.85), "lamy"=5.5, "theta"=0.8)
sigmas = list("lamw"=1, "rhow"=c(1, 1, 1), "lameta"=1,
    "lamv"=1, "rhov"=c(1, 1), "lamy"=1, "theta"=1)

# to decrease acceptance rate, increase sigma
# to increase acceptance rate, decrease sigma

# only if it accepts, allow to increase sigma, not just
# if acceptance rate is too high
# or look at acceptance rate of previous 100 draws.

cov.mat = covmat()  # 176 ^j
post = compute.posterior(cov.mat)
for (iter in 2:niter){

    # f(y|x) = the height at y when the distribution is based off x
    # f(candidate | current) / f(current | candidate)

    # pb.linux(iter-1, niter-1)
    # lam.w
    for (i in 1:p.eta){
        LAMW = row.col(lam.w, apply(lam.w, 2, function(x) max(which(!is.na(x)))))
#   the commented section here is for an example of the changing sigmas idea
#       a = 1/10 # asympototically sigma = 1/(a * sqrt(2*pi))
#       b = 100  # at x = 0, sigma = 1/(b * sqrt(2*pi)), if a=b then sigma is always 1/(a*sqrt(2*pi))
#       d = 10   # higher d means approach to a faster
#       sig.curr = sig.calc(LAMW, a, b, d)
        lam.w[iter] = lam.w[iter-1]
        cand = rnorm(1, LAMW, sigmas[["lamw"]][i])
        if (cand > 0){
#           sig.cand = sig.calc(cand, a, b, d)
            lam.w[iter] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
#           hastings.top = dnorm(cand, LAMW, sig.curr, log=TRUE)
#           hastings.bot = dnorm(LAMW, cand, sig.cand, log=TRUE)
#           if (log(runif(1)) < cand.post - post + hastings.top - hastings.bot){
            if (log(runif(1)) < cand.post - post){
                acc.lam.w = acc.lam.w + 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                lam.w[iter] = lam.w[iter-1]
                }
            }
        }
    if (acc.lam.w/(iter-1) > 0.25)
        sigmas[["lamw"]] = sigmas[["lamw"]] * 1.01
    if (acc.lam.w/(iter-1) < 0.19)
        sigmas[["lamw"]] = sigmas[["lamw"]] * 0.99
    if (post == cand.post)
        cat("\nlam.w  ",trailing(acc.lam.w/(iter-1)),trailing(post),trailing(sigmas[["lamw"]]))

    # lam.v
    for (i in 1:Fgroups){
        LAMV = row.col(lam.v, apply(lam.v, 2, function(x) max(which(!is.na(x)))))
        lam.v[iter] = lam.v[iter-1]
        cand = rnorm(1, LAMV, sigmas[["lamv"]][i])
        if (cand > 0){
            lam.v[iter] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                acc.lam.v = acc.lam.v + 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                lam.v[iter] = lam.v[iter-1]
                }
            }
        }
    if (acc.lam.v/(iter-1) > 0.25)
        sigmas[["lamv"]] = sigmas[["lamv"]] * 1.01
    if (acc.lam.v/(iter-1) < 0.19)
        sigmas[["lamv"]] = sigmas[["lamv"]] * 0.99
    if (post == cand.post)
        cat("\nlam.v  ",trailing(acc.lam.v/(iter-1)),trailing(post),trailing(sigmas[["lamv"]]))

    # lam.y
    LAMY = row.col(lam.y, apply(lam.y, 2, function(x) max(which(!is.na(x)))))
    lam.y[iter] = lam.y[iter-1]
    cand = rnorm(1, LAMY, sigmas[["lamy"]])
    if (cand > 0){
        lam.y[iter] = cand
        cand.covmat = covmat()
        cand.post = compute.posterior(cand.covmat)
        if (log(runif(1)) < cand.post - post){
            acc.lam.y = acc.lam.y + 1
            cov.mat = cand.covmat
            post = cand.post
        } else {
            lam.y[iter] = lam.y[iter-1]
            }
        }
    if (acc.lam.y[1]/(iter-1) > 0.25)
        sigmas[["lamy"]] = sigmas[["lamy"]] * 1.01
    if (acc.lam.y[1]/(iter-1) < 0.19)
        sigmas[["lamy"]] = sigmas[["lamy"]] * 0.99
    if (post == cand.post)
        cat("\nlam.y  ",trailing(acc.lam.y/(iter-1)),trailing(post),trailing(sigmas[["lamy"]]))

    # lam.eta
    LAMETA = row.col(lam.eta, apply(lam.eta, 2, function(x) max(which(!is.na(x)))))
    lam.eta[iter] = lam.eta[iter-1] 
    cand = rnorm(1, LAMETA, sigmas[["lameta"]])
    if (cand > 0){
         lam.eta[iter] = cand
         cand.covmat = covmat()
         cand.post = compute.posterior(cand.covmat)
         if (log(runif(1)) < cand.post - post){
             acc.lam.eta = acc.lam.eta + 1
             cov.mat = cand.covmat
             post = cand.post
         } else {
             lam.eta[iter] = lam.eta[iter-1]
             }
         }
    if (acc.lam.eta/(iter-1) > 0.25)
        sigmas[["lameta"]] = sigmas[["lameta"]] * 1.01
    if (acc.lam.eta/(iter-1) < 0.19)
        sigmas[["lameta"]] = sigmas[["lameta"]] * 0.99
    if (post == cand.post)
        cat("\nlam.eta",trailing(acc.lam.eta/(iter-1)),trailing(post),trailing(sigmas[["lameta"]]))

    # rho.w
    for (i in 1:(px+pt)){
        RHOW = row.col(rho.w, apply(rho.w, 2, function(x) max(which(!is.na(x)))))
        rho.w[iter, i] = rho.w[iter-1, i]
        cand = rnorm(1, RHOW[i], sigmas[["rhow"]][i])
        if (cand >= 0 && cand <= 1){
            rho.w[iter, i] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                acc.rho.w[i] = acc.rho.w[i] + 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                rho.w[iter, i] = rho.w[iter-1, i]
                }
            }
        if (acc.rho.w[i]/(iter-1) > 0.25)
            sigmas[["rhow"]][i] = sigmas[["rhow"]][i] * 1.01
        if (acc.rho.w[i]/(iter-1) < 0.19)
            sigmas[["rhow"]][i] = sigmas[["rhow"]][i] * 0.99
        if (post == cand.post)
            cat("\nrho.w",i,trailing(acc.rho.w[i]/(iter-1)),trailing(post),trailing(sigmas[["rhow"]][i]))
        }

    # rho.v
    for (i in 1:px){
        RHOV = row.col(rho.v, apply(rho.v, 2, function(x) max(which(!is.na(x)))))
        rho.v[iter, i] = rho.v[iter-1, i]
        cand = rnorm(1, RHOV[i], sigmas[["rhov"]][i])
        if (cand >= 0 && cand <= 1){
            rho.v[iter, i] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                acc.rho.v[i] = acc.rho.v[i] + 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                rho.v[iter, i] = rho.v[iter-1, i]
                }
            }
        if (acc.rho.v[i]/(iter-1) > 0.25)
            sigmas[["rhov"]][i] = sigmas[["rhov"]][i] * 1.01
        if (acc.rho.v[i]/(iter-1) < 0.19)
            sigmas[["rhov"]][i] = sigmas[["rhov"]][i] * 0.99
        if (post == cand.post)
            cat("\nrho.v",i,trailing(acc.rho.v[i]/(iter-1)),trailing(post),trailing(sigmas[["rhov"]][i]))
        }

    # theta
    for (i in 1:pt){
        THETA = row.col(theta, apply(theta, 2, function(x) max(which(!is.na(x)))))
        theta[iter, i] = theta[iter-1, i]
        cand = rnorm(1, THETA[i], sigmas[["theta"]][i])
        if (cand >= 0 && cand <= 1){
            theta[iter, i] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                acc.theta[i] = acc.theta[i] + 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                theta[iter, i] = theta[iter-1, i]
                }
            }
        if (acc.theta[i]/(iter-1) > 0.25)
            sigmas[["theta"]][i] = sigmas[["theta"]][i] * 1.01
        if (acc.theta[i]/(iter-1) < 0.19)
            sigmas[["theta"]][i] = sigmas[["theta"]][i] * 0.99
        if (post == cand.post)
            cat("\ntheta",i,trailing(acc.theta[i]/(iter-1)),trailing(post),trailing(sigmas[["theta"]][i]))
        }
    par(mfrow=c(10,1), mar=double(4)+0.2, oma=c(3,5,0,0))
    plot(lam.w[1:iter], type='l', xaxt='n')
    plot(lam.v[1:iter], type='l', xaxt='n')
    plot(lam.y[1:iter], type='l', xaxt='n')
    plot(lam.eta[1:iter], type='l', xaxt='n')
    plot(rho.w[1:iter,1], type='l', xaxt='n')
    plot(rho.w[1:iter,2], type='l', xaxt='n')
    plot(rho.w[1:iter,3], type='l', xaxt='n')
    plot(rho.v[1:iter,1], type='l', xaxt='n')
    plot(rho.v[1:iter,2], type='l', xaxt='n')
    plot(theta[1:iter], type='l', xaxt='n')
    mtext(iter, 1, 2)
    mtext(label, 2, 3, at=seq(0.95, 0.05, length=length(label)),
        outer=TRUE, cex=0.8)
    }

acc.lam.w/niter
acc.lam.v/niter
acc.lam.eta/niter
acc.lam.y/niter
acc.rho.w/niter
acc.rho.v/niter
acc.theta/niter

# high acceptance rates: increase sigma
# i.e., increase sigma to decrease acceptance rate, optimal = 0.22
# i.e., decrease sigma to increase acceptance rate, optimal = 0.22

# save the draws and design matrices
# all = cbind(lam.w, lam.v, lam.eta, lam.y, rho.w, rho.v, theta)
# write.table(all, "./out_all_draws.txt", row.names=FALSE, col.names=FALSE)
# write.table(cbind(sim.X, sim.Y), "./out_simul.txt", row.names=FALSE, col.names=FALSE)
# write.table(cbind(field.X, field.Y), "./out_field.txt", row.names=FALSE, col.names=FALSE)



# make burn-in variables
burn.lam.w = lam.w[101:1000]
burn.lam.v = lam.v[101:1000]
burn.rho.w = rho.w[101:1000,]
burn.rho.v = rho.v[101:1000,]
burn.lam.eta = lam.eta[101:1000]
burn.lam.y = lam.y[101:1000]
burn.theta = theta[101:1000]

# trace
par(mfrow=c(10,1), mar=double(4)+0.2, oma=c(0,3,0,0))
plot(burn.lam.w, type='l', xaxt='n')
plot(burn.lam.v, type='l', xaxt='n')
plot(burn.lam.eta, type='l', xaxt='n')
plot(burn.lam.y, type='l', xaxt='n')
plot(burn.rho.w[,1], type='l', xaxt='n')
plot(burn.rho.w[,2], type='l', xaxt='n')
plot(burn.rho.w[,3], type='l', xaxt='n')
plot(burn.rho.v[,1], type='l', xaxt='n')
plot(burn.rho.v[,2], type='l', xaxt='n')
plot(burn.theta, type='l', xaxt='n')

# density
par(mfrow=c(10,1), mar=double(4)+0.2, oma=c(0,3,0,0))
plot(density(burn.lam.w), type='l', xaxt='n')
plot(density(burn.lam.v), type='l', xaxt='n')
plot(density(burn.lam.eta), type='l', xaxt='n')
plot(density(burn.lam.y), type='l', xaxt='n')
plot(density(burn.rho.w[,1]), type='l', xaxt='n')
plot(density(burn.rho.w[,2]), type='l', xaxt='n')
plot(density(burn.rho.w[,3]), type='l', xaxt='n')
plot(density(burn.rho.v[,1]), type='l', xaxt='n')
plot(density(burn.rho.v[,2]), type='l', xaxt='n')
plot(density(burn.theta), type='l', xaxt='n')

label = c("lam.w", "lam.v", "lam.eta", "lam.y", "rho.w1",
    "rho.w2", "rho.w3", "rho.v1", "rho.v2", "theta")
# highest probability density intervals
burn.all = all[101:1000,]
hpds = t(apply(burn.all, 2, hpd))
rownames(hpds) = label
colnames(hpds) = c("lower", "upper")

# predictions
x.star = matrix(0, 121, 2)
x.star[,1] = rep(seq(0, 1, length=11), each=11)
x.star[,2] = rep(seq(0, 1, length=11), times=11)
w = matrix(0, 900, nrow(x.star))
for (i in 1:900){
    pb.linux(i, 900)
    w[i,] = mvrnorm(1, double(nrow(x.star)), 1/burn.lam.w[i]*R.x(cbind(x.star, burn.theta[i]), burn.rho.w[i,]))
    }
v = matrix(0, 900, nrow(x.star))
for (i in 1:900){
    pb.linux(i, 900)
    v[i,] = mvrnorm(1, double(nrow(x.star)), 1/burn.lam.v[i]*R.x(x.star, burn.rho.v[i,]))
    }

x.star = t(c(0.5,0.5))
preds = 0
for (i in 1:1){
    sig.11 = solve(covmat(i))
    sig.22 = diag(c(1/burn.lam.v[i], 1/burn.lam.w[i]))
    sig.12 = cbind(1/burn.lam.v[i]*R.cross(z.vec, x.star, burn.rho.v[i,]),
        1/burn.lam.w[i]*R.cross(z.vec, x.star, burn.rho.w[i,]))
    preds = mvrnorm(1, t(sig.12)%*%sig.11%*%z.vec, sig.22-t(sig.12)%*%sig.11%*%sig.12)
    }


# [ 230 x 230   230 x 2 ]
# [   2 x 230     2 x 2 ]

w.star = t(apply(w, 1, function(x) x*ysim.sd + ysim.mean))
#v.star = t(apply(v, 1, function(x) x*ysim.sd + ysim.mean))
new.star = w.star*ysim.sd + ysim.mean + v.star
# back transform the simulator part, not the discrepancy part

meanit = apply(new.star, 2, mean)
one = cbind(x.star, meanit)

calibrated = apply(cbind(simul[,1:2], 0.5), 1, sim.f)
plot3d(cbind(simul[,1:2], calibrated))
points3d(one, col='red')
points3d(field, col='blue')
