library(MASS)
dat = read.csv("./gdp_data.csv")

set.seed(1)
index = sample(60, 10)
index = 1:60
Y = dat[,3]           # growth (or decline) in GDP from 1960-1996
X = dat[,-c(1,2,3)]   # remove country columns and response
X = cbind(X, ifelse(X$NEWSTATE == 0, 1, 0))
X = cbind(X, ifelse(X$NEWSTATE == 1, 1, 0))
names(X)[c(68, 69)] = c("NEWSTATE0", "NEWSTATE1")
X = X[,-39]
for (i in 0:4)
    X = cbind(X, ifelse(X$ECORG == i, 1, 0))
names(X)[c(69,70,71,72,73)] = c("ECORG0", "ECORG1", "ECORG2",
    "ECORG3", "ECORG4")
X = X[,-15]
label = c("Intercept", names(X))
X = as.matrix(X)
for (i in 1:ncol(X)){
    if (length(unique(X[,i]))!=2)
        X[,i] = scale(X[,i])
    }
X = cbind(1 ,X)

Y.test = Y[index]
Y.train = Y[-index]
X.test = X[index,]
X.train = X[-index,]

n = nrow(X)
m = 10 # size of hold-out
k = ncol(X)
a = 0
b = 0

#lams = as.matrix(read.table("./lam_values.txt", header=F))
#mses = as.matrix(read.table("./mse_values.txt", header=F))
#mses = apply(mses, 2, mean)
#lambda = lams[which.min(mses)]
lambda = 322


niter = 100000
betas = matrix(0, niter, k)
sig2s = double(niter)

sig2s[1] = var(Y)
# betas[1,] = mvrnorm(1, solve(1/sig2s[1]*t(X.train)%*%X.train+diag(lambda, k))%*%t(X.train)%*%Y.train/sig2s[1],
#     solve(1/sig2s[1]*t(X.train)%*%X.train+diag(lambda, k)))
betas[1,] = double(k)

# Gibbs sampler
# for (i in 2:niter){
#     sig2s[i] = 1/rgamma(1, (n-m)/2 + a, 0.5*t(Y.train-X.train%*%betas[i-1,])%*%(Y.train-X.train%*%betas[i-1,])+b)
#     A = solve(1/sig2s[i]*t(X.train)%*%X.train+diag(lambda, k))
#     betas[i,] = mvrnorm(1, A%*%t(X.train)%*%Y.train/sig2s[i], A)
#     }
Y = as.matrix(Y)
for (i in 2:niter){
    sig2s[i] = 1/rgamma(1, n/2 + a, 0.5*t(Y-X%*%betas[i-1,])%*%(Y-X%*%betas[i-1,])+b)
    A = solve(1/sig2s[i]*t(X)%*%X+diag(lambda, k))
    betas[i,] = mvrnorm(1, A%*%t(X)%*%Y/sig2s[i], A)
    }

# calculate MSE
MSE = double(length(Y.test))
for (l in 1:length(Y.test)){
    MSE.temp = double(niter)
    for (i in 1:niter)
        MSE.temp[i] = rnorm(1,X.test[l,] %*% betas[i,],sd=sqrt(sig2s[i]))
    MSE[l] = (Y.test[l] - mean(MSE.temp))^2
    }

MSE = mean(MSE)
apply(betas, 2, quantile, c(0.025,0.975))

orig.dist = apply(apply(betas, 2, function(x) ifelse(x>0, 1, 0)), 2, mean)
distance = abs(apply(apply(betas, 2, function(x) ifelse(x>0, 1, 0)), 2, mean)-0.5)
sort(out)

func.mode = function(x, precision=512)
    density(x, n=precision)$x[which.max(density(x, n=precision)$y)]

ests = apply(betas, 2, func.mode)
ests = cbind(ests, 3)
ests[,2] = ifelse(distance > 0.15, 2, ests[,2])
ests[,2] = ifelse(distance > 0.25, 1, ests[,2])
out.coef = cbind(label, ests)
ests = ests[order(-distance),]
out.coef = out.coef[order(-distance),]
write.table(out.coef, "./coefs.txt", row.names=FALSE, col.names=FALSE)

# cont: one standard deviation increase in X results in a one beta increase in Y
# cate: when X is 1, we expect Y to increase by beta

best = order(-distance)[1:7]

# get the closest values in density() for the intervals
bounds = function(int, dens, x=TRUE){
    if (x) # returns x-value from density
        return(dens$x[which.min(abs(dens$x-int))])
    return(dens$y[which.min(abs(dens$x-int))])
    }

pdf("./figs/betas.pdf")
par(mfrow=c(7,1), mar=double(4)+1.1, oma=c(1,4,1,1))
for (i in 1:7){
    dens = density(betas[,best[i]])
    plot(density(betas[,best[i]]), ylab="", xlab="", main="")
    color = "gold2"
    if (orig.dist[best[i]] < 0.5){
        at = which(dens$x - bounds(0, dens, TRUE) == 0)
        polygon(c(dens$x[1:at], 0), c(dens$y[1:at],0), col=color)
    } else {
        at = which(dens$x - bounds(0, dens, TRUE) == 0)
        polygon(c(dens$x[at:512], 0), c(dens$y[at:512],0), col=color)
        }
    mtext(label[best[i]], side=2, line=2.8)
    }
dev.off()



out.sig2 = func.mode(sig2s, 10000)
prob.int = quantile(sig2s, c(0.025, 0.975))
pdf("./sig_post.pdf")
plot(density(sig2s, n=100000), main=expression(paste("Posterior Distribution of ",
    sigma,""^2)), xlab = expression(paste(sigma,""^2)), xlim=c(-1e-5, 3e-5))
polygon(c(0, density(sig2s, n=100000)$x, 0), c(0, density(sig2s, n=100000)$y, 0), col='cyan3')
dev.off()
plot(density(betas[,71]), type='l')


lambda = 322
out.sige2 = 0.000004231827
