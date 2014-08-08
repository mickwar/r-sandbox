library(MASS)
dat = read.csv("../gdp_data.csv")

set.seed(1)
index = sample(60, 10)
Y = dat[,3]             # growth (or decline) in GDP from 1960-1996
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
X = as.matrix(X)
for (i in 1:ncol(X)){
    if (length(unique(X[,i]))!=2)
        X[,i] = scale(X[,i])
    }
X = cbind(1 ,X)
Y.test.set = Y[index]
Y = Y[-index]
X.test.set = X[index,]
X = X[-index,]

n = nrow(X)
m = 5 # size of hold-out
k = ncol(X)
a = 5
b = 5

library(foreach)
library(doMC)
registerDoMC(32)

doit = function(kfold){
    lambda = seq(10, 100, length=100)
    MSE.total = double(length(lambda))

    train = as.matrix(read.table("./training.txt", header=F))
    test.index = as.numeric(train[kfold,])

    X.train = X[-test.index,]
    Y.train = as.matrix(Y[-test.index])
    X.test = X[test.index,]
    Y.test = as.matrix(Y[test.index])
 
    for (j in 1:length(lambda)){
        niter = 10000
        betas = matrix(0, niter, k)
        sig2s = double(niter)

        sig2s[1] = 0.1
        betas[1,] = mvrnorm(1, solve(1/sig2s[1]*t(X.train)%*%X.train+diag(lambda[j], k))%*%t(X.train)%*%Y.train/sig2s[1],
            solve(1/sig2s[1]*t(X.train)%*%X.train+diag(lambda[j], k)))
        
        # Gibbs sampler
        for (i in 2:niter){
            sig2s[i] = 1/rgamma(1, m/2 + a, 0.5*t(Y.train-X.train%*%betas[i-1,])%*%(Y.train-X.train%*%betas[i-1,])+b)
            A = solve(1/sig2s[i]*t(X.train)%*%X.train+diag(lambda[j], k))
            betas[i,] = mvrnorm(1, A%*%t(X.train)%*%Y.train/sig2s[i], A)
            }

        # calculate MSE
        MSE = double(length(Y.test))
        for (l in 1:length(Y.test)){
            MSE.temp = double(niter)
            for (i in 1:niter)
                MSE.temp[i] = rnorm(1,X.test[l,] %*% betas[i,],sd=sqrt(sig2s[i]))
            MSE[l] = (Y.test[l] - mean(MSE.temp))^2
            }

        MSE.total[j] = mean(MSE)
        }
    return(t(MSE.total))
    }

out = foreach(i=1:10,.combine=rbind) %dopar% doit(i)
write.table(out, "~/files/stat536/LAM_4.txt", col.names=FALSE, row.names=FALSE)
