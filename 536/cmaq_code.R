# gaussian process of several variables (frequentist)
# coverage
# parallelization


library(LatticeKrig)
library(geoR)
library(rgl)
library(foreach)
library(doMC)
registerDoMC(2)
source("~/files/R/pb_linux.R")

ozone = read.csv("./Ozone.csv")
cmaq = read.csv("./CMAQ.csv")
predlocs = read.csv("./PredLocs.csv")

# remove useless variables and reorder some things
ozone = ozone[,-c(2,1,3)]
# longitude, latitude, ozone level
ozone = ozone[,c(2,1,3)]
cmaq = cmaq[,-3]

predlocs = predlocs[,c(3,2)]
names(predlocs) = c("Longitude", "Latitude")

### quiltplot stuff
pdf("./figs/data.pdf", height=6, width=14)
par(mfrow=c(1,2), mar=c(5,5,3,4.5))
zrange = range(ozone[,3]) # range for cmaq is contained in zrange
zrange[1] = min(zrange[1], 0)
zrange[2] = max(zrange[2], 100)
quilt.plot(ozone[,1], ozone[,2], ozone[,3], xlab="Longitude",
    main="Station Measurements", ylab="Latitude", zlim=zrange)
map("state", add=TRUE)
quilt.plot(cmaq[,1], cmaq[,2], cmaq[,3], xlab="Longitude",
    main="CMAQ Measurements", ylab="", zlim=zrange)
map("state", add=TRUE)
dev.off()

### 3d plot stuff
plot3d(cmaq[sample(nrow(cmaq), 10000),])
points3d(ozone, col='red', size=5)


nu = 1
s2.start.val = 35
phi.start.val = 1.5

# get P closest CMAQ values
# 1 degree longitude is about 54 miles
# 1 degree latitude is about 68 miles
# half degree radius should be about 30 miles or so
all.dist = rdist(ozone[,1:2], cmaq[,1:2])




at = 500
tot.X = matrix(0, 800, at)
for (i in 1:800){
    temp_index = which(all.dist[i,] <= 2.5)
    temp_vals = all.dist[i, temp_index]
    tot.X[i,] = cmaq[temp_index[order(temp_vals)[1:at]],3]
#   points3d(cmaq[temp_index[order(temp_vals)[1:at]],], col='green', size=10)
#   points3d(ozone[i,], col='blue', size=10)
    }


### parallelization to compute mean squared error for
### given number of parameters

mse.par = function(P){
    cat(P," ")
    N = 800
    X = cbind(1, tot.X[,1:P])
    field.dat = as.geodata(ozone, data.col=3,
        coords.col=c(1,2))

    gp.field = likfit(field.dat, cov.model="matern", kappa=nu,
        fix.kappa=TRUE, ini.cov.pars=c(s2.start.val, phi.start.val),
        trend=~X-1, messages=FALSE)

    phi = 1/gp.field$phi
    s2 = gp.field$sigmasq
    beta = gp.field$beta
    tau2 = gp.field$tausq

    D = rdist(ozone[,1:2])
    V = s2*Matern(D, alpha=phi, nu=nu)
    inv.V = solve(V + diag(tau2, N))
    predictions = X%*%beta + V%*%inv.V%*%(ozone[,3]- X%*%beta)
    mean((predictions - ozone[,3])^2) # mse output
    }

mse.out=foreach(i=1:100,.combine=rbind)%dopar%mse.par(i)
write.table(mse.out, "./mse.txt", row.names=FALSE, col.names=FALSE)

plot(mse.out)

P = which.min(mse.out)
# about what distance does P correspond to?
# find d where min(counts) ~ P
# what does a value of d=1 translate to miles?
# is that a reasonale stopping point?
counts = double(800)
d = 0.5
for (i in 1:800)
    counts[i] = sum(all.dist[i,] <= d)
range(counts)

P = 299

### parallelization to compute coverage
coverage.par = function(iter){
    cat(iter," ")
    N = 600
    K = 800 - N
    training = sort(sample(800, N))
    testing = (1:800)[-training]
    train.X = cbind(1, tot.X[training,1:P])
    pred.X = cbind(1, tot.X[testing,1:P])

    field.dat = as.geodata(ozone[training,], data.col=3,
        coords.col=c(1,2))

    gp.field = likfit(field.dat, cov.model="matern", kappa=nu,
        fix.kappa=TRUE, ini.cov.pars=c(s2.start.val, phi.start.val),
        trend=~train.X-1, messages=FALSE)

    phi = 1/gp.field$phi
    s2 = gp.field$sigmasq
    beta = gp.field$beta
    tau2 = gp.field$tausq

    D = rdist(rbind(ozone[testing,1:2], ozone[training,1:2]))
    V = s2*Matern(D, alpha=phi, nu=nu) ##V = Sigma_Y
    inv.V = solve(V[K+(1:N),K+(1:N)]+diag(tau2,N))
    EV = pred.X%*%beta + V[1:K,K+(1:N)]%*%inv.V%*%(ozone[training,3]-
        train.X%*%beta)
    cond.Var = diag((V[1:K,1:K]+diag(tau2, K))-V[1:K,K+(1:N)]%*%
        inv.V%*%t(V[1:K,K+(1:N)]))
    upper = qnorm(0.975,mean=EV,sd=sqrt(cond.Var))
    lower = qnorm(0.025,mean=EV,sd=sqrt(cond.Var))
    mean(apply(cbind(lower <= ozone[testing,3], upper >=
        ozone[testing, 3]), 1, all))
    }

coverage.out=foreach(i=1:10,.combine=rbind)%dopar%coverage.par(i)
write.table(coverage.out, "./coverage.txt", row.names=FALSE, col.names=FALSE)

cover = as.vector(read.table("./coverage.txt"))
cover = cover[,1]
pdf("./figs/coverage.pdf")
hist(cover,col="gray",main="Coverage",freq=FALSE,xlab="Coverage")
dev.off()


N = nrow(ozone)
K = nrow(predlocs)
P = 30
X = cbind(1, tot.X[,1:P])

field.dat = as.geodata(ozone, data.col=3,
    coords.col=c(1,2))

gp.field = likfit(field.dat, cov.model="matern", kappa=nu,
    fix.kappa=TRUE, ini.cov.pars=c(s2.start.val, phi.start.val),
    trend=~X-1, messages=FALSE)

phi = 1/gp.field$phi
s2 = gp.field$sigmasq
beta = gp.field$beta
tau2 = gp.field$tausq

### GP conditioned on observations only
D = rdist(ozone[,1:2])
V = s2*Matern(D, alpha=phi, nu=nu) ##V = Sigma_Y
inv.V = solve(V+diag(tau2, N))
EV = X%*%beta + V%*%inv.V%*%(ozone[,3]-X%*%beta)
cond.Var = diag((V+diag(tau2, N))-V%*%inv.V%*%t(V))
upper = qnorm(0.975,mean=EV,sd=sqrt(cond.Var))
lower = qnorm(0.025,mean=EV,sd=sqrt(cond.Var))

std.err = sqrt(diag(solve(t(X)%*%inv.V%*%X)))

# 1: estimates, 2: std.errs, 3: lower 95%, 4: upper 95%
ests = cbind(0, std.err-beta, -qnorm(0.975)*std.err, qnorm(0.975)*std.err)+beta

Rsq = 1 - sum((ozone[,3]-EV)^2)/sum((ozone[,3]-mean(ozone[,3]))^2)
Radj = 1 - (1 - Rsq) * (N-1) / (N-P-1)
# R^2 = 0.922271
# R_adj = 0.9192387

# residuals
pdf("./figs/residuals.pdf")
plot(scale(EV - ozone[,3]), main="Residuals", xlab="",
    ylab="Standardized Residuals", pch=20)
abline(h=c(-3,3), lty=2, col='gray')
dev.off()

write.table(EV, "./ests.txt", row.names=FALSE, col.names=FALSE)


### Predictions
pred.dist = rdist(predlocs, cmaq[,1:2])
tot.pred.X = matrix(0, K, 100)
for (i in 1:K){
    temp_index = which(pred.dist[i,] <= 0.7)
    temp_vals = pred.dist[i, temp_index]
    tot.pred.X[i,] = cmaq[temp_index[order(temp_vals)[1:100]],3]
    }
pred.X = cbind(1, tot.pred.X[,1:P])


D = rdist(rbind(predlocs, ozone[,1:2]))
V = s2*Matern(D, alpha=phi, nu=nu) ##V = Sigma_Y
inv.V = solve(V[K+(1:N),K+(1:N)]+diag(tau2,N))
pred.EV = pred.X%*%beta + V[1:K,K+(1:N)]%*%inv.V%*%(ozone[,3]-X%*%beta)
cond.Var = diag((V[1:K,1:K]+diag(tau2, K))-V[1:K,K+(1:N)]%*%
    inv.V%*%t(V[1:K,K+(1:N)]))
upper = qnorm(0.975,mean=pred.EV,sd=sqrt(cond.Var))
lower = qnorm(0.025,mean=pred.EV,sd=sqrt(cond.Var))

### quiltplot stuff (on predictions)
pdf("./figs/preds.pdf", height=6, width=14)
par(mfrow=c(1,2), mar=c(5,5,3,4.5))
zrange = range(ozone[,3]) # range for pred.EV is contained in zrange
zrange[1] = min(zrange[1], 0)
zrange[2] = max(zrange[2], 100)
quilt.plot(ozone[,1], ozone[,2], ozone[,3], xlab="Easter Egg",
    main="Station Measurements", ylab="Latitude", zlim=zrange)
map("state", add=TRUE)
quilt.plot(predlocs[,1], predlocs[,2], pred.EV, xlab="Longitude",
    main="Predictions", ylab="", zlim=zrange)
map("state", add=TRUE)
dev.off()

pdf("./figs/predsun.pdf", height=6, width=14)
par(mfrow=c(1,2), mar=c(5,5,3,4.5))
zrange = range(ozone[,3]) # range for pred.EV is contained in zrange
zrange[1] = min(zrange[1], 0)
zrange[2] = max(zrange[2], 100)
quilt.plot(predlocs[,1], predlocs[,2], lower, xlab="Longitude",
    main="Lower 95% Predictions", ylab="Latitude", zlim=zrange)
map("state", add=TRUE)
quilt.plot(predlocs[,1], predlocs[,2], upper, xlab="Longitude",
    main="Upper 95% Predictions", ylab="", zlim=zrange)
map("state", add=TRUE)
dev.off()


par(mfrow=c(1,2), mar=c(6,3,6,1))
quilt.plot(pred[,1], pred[,2], EV)
plot3d(pred[,1], pred[,2], EV)
map("state", add=TRUE)
points3d(pred[,1], pred[,2], newEV, col='red')
map("state", add=TRUE)

points3d(ozone[,1], ozone[,2], ozone[,3], col='blue', size=5)
map("state", add=TRUE)



matern.eff.range = function(phi, nu, cor.target = 0.05, guess = 1,
    eps = 1e-6){
    d = guess
    out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
        besselK(2*phi*sqrt(nu)*d, nu)
    adjust = 1
    while (out > cor.target + eps){
        d = d + adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out < cor.target){
            d = d - adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    adjust = 1
    if (d <= 1)
        d = d + 1.1
    while (out < cor.target - eps){
        d = d - adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out > cor.target){
            d = d + adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    return (d)
    }
2*matern.eff.range(phi, nu, guess=6, eps=1e-6)
2*matern.eff.range(phi*2, nu, guess=6, eps=1e-6)



matern.cor.to.phi = function(d, nu, cor.target = 0.05, guess = 1,
    eps = 1e-6){
    phi = guess
    out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
        besselK(2*phi*sqrt(nu)*d, nu)
    adjust = 1
    while (out > cor.target + eps){
        phi = phi + adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out < cor.target){
            phi = phi - adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    adjust = 1
    if (phi <= 1)
        phi = phi + 1.1
    while (out < cor.target - eps){
        phi = phi - adjust
        out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
            besselK(2*phi*sqrt(nu)*d, nu)
        if (out > cor.target){
            phi = phi + adjust
            adjust = adjust/10
            out = 1/(gamma(nu)*2^(nu-1)) * (2*phi*sqrt(nu)*d)^nu *
                besselK(2*phi*sqrt(nu)*d, nu)
            }
        }
    return (phi)
    }
matern.cor.to.phi(5, 1)
