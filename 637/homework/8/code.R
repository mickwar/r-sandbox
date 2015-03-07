source("~/files/R/mcmc/bayes_functions.R")
library(MASS)
library(truncnorm)

load("~/files/data/casella.Rdata")

new.train = train[,-which(names(train) %in% c("title.review", "review", "name"))]
new.train[is.na(new.train)] = 0

n = NROW(new.train)
Y = new.train$stars.scores
J = max(Y)
# temp = matrix(0, n, 5)
# for (i in 1:n)
#     temp[i, Y[i]] = 1
# Y = temp
# rm(temp)

X = cbind(1, new.train$help.yes, new.train$help.total, new.train$year)
p = ncol(X)
matrix(scale(X), n, p)
X = matrix(scale(X), n, p)
X[,1] = 1

xtx = solve(t(X) %*% X)
xtxx = xtx %*% t(X)


### parameters
nburn = 100
ngibbs = 10000

param.beta = matrix(0, p, nburn + ngibbs)
param.gamma = matrix(0, J-1, nburn + ngibbs)
param.gamma[,1] = sort(runif(J-1, -2, 2))
#matrix(rep(sort(runif(J-1, -2, 2)), 7), J-1, 7)
#param.gamma = matrix(rep(sort(runif(J-1, -2, 2)), nburn + ngibbs), J-1, nburn + ngibbs)

for (j in 1:n){
    if (Y[j] == 1)
        param.z[j,1] = runif(1, param.gamma[Y[j],1]-1, param.gamma[Y[j],1])
    if (Y[j] == J)
        param.z[j,1] = runif(1, param.gamma[Y[j]-1,1], param.gamma[Y[j]-1,1]+1)
    if (!(Y[j] == 1 || Y[j] == J))
        param.z[j,1] = runif(1, param.gamma[Y[j]-1,1], param.gamma[Y[j],1])
    }

z.lower = double(n)
z.upper = double(n)

### gibbs sampling
for (i in 2:(nburn + ngibbs)){
    cat("Iteration:", i, "/", nburn + ngibbs, "\n")
    # beta
    param.beta[,i] = mvrnorm(1, xtxx %*% param.z[,i-1], xtx)

    # z
    for (j in 1:n){
        if (Y[j] == 1){
            z.lower[j] = -Inf
        } else {
            z.lower[j] = param.gamma[Y[j]-1, i-1]
            }
        if (Y[j] == J){
            z.upper[j] = Inf
        } else {
            z.upper[j] = param.gamma[Y[j], i-1]
            }
        }
    param.z[,i] = rtruncnorm(1, mean = X %*% param.beta[,i], sd = 1,
        a = z.lower, b = z.upper)

    # gamma
    z.uniq = sort(unique(c(z.lower, z.upper)))
    for (j in 1:(J-1)){
        param.gamma[j, i] = runif(1,
            min = max(max(param.z[Y == j, i]), z.uniq[j]),
            max = min(min(param.z[Y == j+1, i]), z.uniq[j+1]))
        }
    }

par(mfrow=c(3,1))
plot(param.beta[3,], type='l')
plot(param.z[6,], type='l')
plot(param.gamma[4,], type='l')



