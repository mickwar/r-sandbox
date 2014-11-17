library(rgl)
library(MASS)
#library(Matrix)
source("~/files/R/omlhd.R")
source("~/files/R/bayes_functions.R")
source("./construct_functions.R")
source("./group_functions.R")

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
field.X = (omlhd(n, 2) - 1)/n
field.Y = apply(field.X, 1, function(x) sim.f(c(x, 0.5)) + eps.f(1) + disc.f(x))

m = 200
sim.X = (omlhd(m, 3) - 1)/m
sim.Y = apply(sim.X, 1, sim.f)

label = c("lamW", "lamV", "lamY", "lamEta", "rhoW1", "rhoW2",
    "rhoW3", "rhoV1", "rhoV2", "theta")

# write X
write.table(field.X, "./out_field_x.txt", quote=FALSE,
    row.names=FALSE, col.names=FALSE)
write.table(sim.X, "./out_simul_x.txt", quote=FALSE,
    row.names=FALSE, col.names=FALSE)

# loading previous data
n = 15
m = 200
burn.draws = read.table("./out_burn.txt", header=TRUE)
sims = read.table("./out_simul.txt")
sim.X = sims[,1:3]
sim.Y = sims[,4]
fields = read.table("./out_field.txt")
field.X = fields[,1:2]
field.Y = fields[,3]
ysim.mean = mean(sim.Y)
ysim.sd = sd(sim.Y)
label = c("lamW", "lamV", "lamY", "lamEta", "rhoW1", "rhoW2",
    "rhoW3", "rhoV1", "rhoV2", "theta")

# higdon model
# 2.2.1 design of simulator runs
# inputs are already in [0,1]^3
px = 2
pt = 1

####### load chemistry data
field = matrix(scan("~/files/R/651/final_project/data/field.txt"), 9, 4, byrow=TRUE)
simul = matrix(scan("~/files/R/651/final_project/data/simul.txt"), 693, 13, byrow=TRUE)

field.x = field[,1]
field.y = field[,2:4]

simul.x = simul[,1]
simul.t = simul[,2:10]
simul.y = simul[,11:13]

# set constant variables
n = nrow(field)
m = nrow(simul)
px = 1
pt = ncol(simul.t)
p = px + pt
q = ncol(field.y)

### transform the data
fdat.y = logit(field.y)
sdat.y = logit(simul.y)
y.mean = apply(sdat.y, 2, mean)
y.sd = apply(sdat.y, 2, sd)

# standardize output based on simulation mean and sd
fdat.y = (fdat.y - matrix(rep(y.mean, n), n, q, byrow=TRUE)) /
    matrix(rep(y.sd, n), n, q, byrow=TRUE)
sdat.y = (sdat.y - matrix(rep(y.mean, m), m, q, byrow=TRUE)) /
    matrix(rep(y.sd, m), m, q, byrow=TRUE)
#######

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

##### EXPAND THIS SECTION
# 2.2.3 discrepancy model
Fgroups = 1
field.yStd = (field.Y - ysim.mean)/ysim.sd
p.delta = 1
D = diag(15)


#########################

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
nmcmc = 1000
nburn = 1000

# initialize paramters, first rows are starting values
draws = NULL
accepts = NULL
sigmas = NULL
for (var.name in label){
    draws[[var.name]] = NA*double(nmcmc)
    draws[[var.name]][1] = 0.5
    accepts[[var.name]] = double(nmcmc)
    sigmas[[var.name]] = 1
    }
upper.bound = list("lamW"=Inf, "lamV"=Inf, "lamY"=Inf, "lamEta"=Inf,
    "rhoW1"=1, "rhoW2"=1, "rhoW3"=1, "rhoV1"=1, "rhoV2"=1, "theta"=1)
params = list("lamW"=c(1.1, 1.1), "lamV"=c(1.1, 1.1), "lamY"=c(1.1, 1.1),
    "lamEta"=c(1.1, 1.1), "rhoW1"=c(1.1, 1.1), "rhoW2"=c(1.1, 1.1),
    "rhoW3"=c(1.1, 1.1), "rhoV1"=c(1.1, 1.1), "rhoV2"=c(1.1, 1.1),
    "theta"=c(1.1, 1.1))
# using result 1
#params[["lamEta"]] = params[["lamEta"]] + c(m*(n.eta-p.eta)/2,
#    0.5*t(sim.yStd)%*%(diag(m)-K%*%solve(t(K)%*%K)%*%t(K))%*%sim.yStd)
#params[["lamY"]] = params[["lamY"]] + c((n-rankMatrix(B))/2,
#    0.5*t(field.yStd)%*%(W-W%*%B%*%solve(t(B)%*%W%*%B + diag(ridge, n*2))%*%t(B)%*%W)%*%field.yStd)

window = 100

cov.mat = covmat()
post = compute.posterior(cov.mat)
for (iter in 2:nmcmc){
    for (var.name in label){
        hold = draws[[var.name]][iter-1]
        draws[[var.name]][iter] = hold
        cand = rnorm(1, hold, sigmas[[var.name]])
        if (cand >= 0 && cand <= upper.bound[[var.name]]){
            draws[[var.name]][iter] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                accepts[[var.name]][iter] = 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                draws[[var.name]][iter] = hold
                }
            }
        if (floor(iter/window) == iter/window && iter < burnin){
            sigmas[[var.name]] = sigmas[[var.name]]*autotune(
                mean(accepts[[var.name]][(iter-window+1):iter]),
                k = max(window/50, 1.1))
            }
        if (post == cand.post){
            cat("\n", var.name, rep(" ", 7-nchar(var.name)),
                trailing(mean(accepts[[var.name]][(window*
                floor(iter/window)+1):iter])), " ", trailing(post),
                " ", trailing(sigmas[[var.name]]), sep="")
            }
        }

    par(mfrow=c(10,1), mar=double(4)+0.2, oma=c(3,5,0,0))
    for (var.name in label){
        plot(new.seq(iter, window*2), draws[[var.name]][new.seq(
            iter, window*2)], type='l', xaxt='n')
        abline(v=seq(window, nburn+nmcmc, by=window), lty=2, col='gray')
        }
    mtext(iter, 1, 2)
    mtext(label, 2, 3, at=seq(0.95, 0.05, length=length(label)),
        outer=TRUE, cex=0.8)
    }

# change the code to combine burnin loop with posterior draws
# (set each dampen = Inf or put a flag in the if statements
# for when changing sigma, change draws[[i]] to include burnin
# and posterior draws)

burnin = 500
means = lapply(draws, function(x){mean(x[burnin:niter])})

# initialize paramters, first rows are starting values
niter = 10000
draws = NULL
accepts = NULL
for (var.name in label){
    draws[[var.name]] = NA*double(niter)
    draws[[var.name]][1] = means[[var.name]]
    accepts[[var.name]] = double(niter)
    }

cov.mat = covmat()
post = compute.posterior(cov.mat)
for (iter in 2:niter){
    for (var.name in label){
        hold = draws[[var.name]][iter-1]
        draws[[var.name]][iter] = hold
        cand = rnorm(1, hold, sigmas[[var.name]])
        if (cand >= 0 && cand <= upper.bound[[var.name]]){
            draws[[var.name]][iter] = cand
            cand.covmat = covmat()
            cand.post = compute.posterior(cand.covmat)
            if (log(runif(1)) < cand.post - post){
                accepts[[var.name]][iter] = 1
                cov.mat = cand.covmat
                post = cand.post
            } else {
                draws[[var.name]][iter] = hold
                }
            }
        if (post == cand.post){
            cat("\n", var.name, rep(" ", 7-nchar(var.name)),
                trailing(mean(accepts[[var.name]][1:iter])), " ",
                trailing(post), sep="")
            }
        }

    par(mfrow=c(10,1), mar=double(4)+0.2, oma=c(3,5,0,0))
    for (var.name in label)
        plot(draws[[var.name]][new.seq(iter, 1000)], type='l', xaxt='n')
    mtext(iter, 1, 2)
    mtext(label, 2, 3, at=seq(0.95, 0.05, length=length(label)),
        outer=TRUE, cex=0.8)
    }

lapply(accepts, mean)
lapply(draws, mean)
lapply(draws, quantile, c(0.025, 0.975))



# save the draws and design matrices
write.table(burn.draws, "./out_burn.txt", row.names=FALSE)
write.table(cbind(sim.X, sim.Y), "./out_simul.txt", row.names=FALSE, col.names=FALSE)
write.table(cbind(field.X, field.Y), "./out_field.txt", row.names=FALSE, col.names=FALSE)

burn.draws = read.table("./out_burn.txt", header=TRUE)
sims = read.table("./out_simul.txt")
sim.X = sims[,1:3]
sim.Y = sims[,4]
fields = read.table("./out_field.txt")
field.X = fields[,1:2]
field.Y = fields[,3]
ysim.mean = mean(sim.Y)
ysim.sd = sd(sim.Y)

burn = 500
nmcmc = niter-burn

plot(density(burn.draws[["theta"]]))

plot3d(cbind(field.X, field.Y))



# highest probability density intervals
hpds = lapply(burn.draws, hpd)



# predictions
x.star = matrix(0, 121, 2)
x.star[,1] = rep(seq(0, 1, length=11), each=11)
x.star[,2] = rep(seq(0, 1, length=11), times=11)
x.star = matrix(c(0.2, 4/15, 0.8, 4/15), 2, 2)
npreds = 1
w = rep(list(matrix(0, npreds, nrow(x.star))), nmcmc)
w = matrix(0, nmcmc, nrow(x.star))
for (i in 1:nmcmc){
    pb.linux(i, nmcmc)
    w[i,] = 
mvrnorm(1, double(nrow(x.star)), 1/burn.draws[["lamW"]][i]*
        R.x(cbind(x.star, burn.draws[["theta"]][i]), c(burn.draws[["rhoW1"]][i],
        burn.draws[["rhoW2"]][i], burn.draws[["rhoW3"]][i])))
#   for (j in 1:npreds){
#       pb.linux(i, npreds)
#       w[[i]][j,] = t(chol.Decomp) %*% rnorm(nrow(x.star))
#       }
    }
v = matrix(0, nmcmc, nrow(x.star))
for (i in 1:nmcmc){
    pb.linux(i, nmcmc)
    v[i,] = mvrnorm(1, double(nrow(x.star)), 1/burn.draws[["lamV"]][i]*R.x(x.star,
        c(burn.draws[["rhoV1"]][i], burn.draws[["rhoV2"]][i])))
    }

w = apply(w, 2, mean)
v = apply(v, 2, mean)
z = w*ysim.sd + v + ysim.mean

plot(double(nmcmc)+1,z[,1], xlim=c(1, 2))
points(double(nmcmc)+2,z[,2])
points(1, field.Y[14], col='red')
points(2, field.Y[15], col='red')

plot3d(cbind(x.star, z))
points3d(cbind(field.X, field.Y), col='red')

x.star = t(c(0.5,0.5))
preds = 0
for (i in 1:1){
    sig.11 = solve(covmat(i))
    sig.22 = diag(c(1/burn.draws[["lamV"]][i], 1/burn.draws[["lamW"]][i]))
    sig.12 = cbind(1/burn.draws[["lamV"]][i]*R.cross(cbind(z.vec, z.vec), x.star, c(burn.draws[["rhoV1"]][i],
        burn.draws[["rhoV2"]][i])), 1/burn.draws[["lamW"]][i]*R.cross(cbind(z.vec, z.vec, z.vec), cbind(x.star, burn.draws[["theta"]][i]),
        c(burn.draws[["rhoW1"]][i], burn.draws[["rhoW2"]][i], burn.draws[["rhoW3"]][i])))
    preds = mvrnorm(1, t(sig.12)%*%sig.11%*%z.vec, sig.22-t(sig.12)%*%sig.11%*%sig.12)
    


    }


# [ 230 x 230   230 x 2 ]
# [   2 x 230     2 x 2 ]

w.star = t(apply(w, 1, function(x) x*ysim.sd + ysim.mean))
v.star = t(apply(v, 1, function(x) x*ysim.sd + ysim.mean))
new.star = v.star - 0*w.star

i = 1
sig.11inv = solve(covmat(i))
sig.12 = 
sig.22 = matrix.assemble(c(2,2), diag(burn.draws[["lamV"]][i], p.delta), NA, NA,
    diag(burn.draws[["lamW"]][i], p.eta))


