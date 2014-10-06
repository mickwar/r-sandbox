rm(list=ls())
# make the eigen values positive in covariance function
pd.sig = function(x){
  y = eigen(x, symmetric = TRUE)
  if (any(y$values < 0)){
    y$values[y$values < 0] = min(y$values[y$values > 0])
    return (y$vectors %*% diag(y$values) %*% t(y$vectors))
  } else {
    return (x)
  }
}

# expectation-maximization routine
em = function(x, tol = 0.01){
  miss = which(apply(x, 1, function(x) any(is.na(x))))
  if (length(miss) == 0)
    stop("Matrix x contains no missing values.")
  p = rep(list(0), length(miss))
  q = rep(list(0), length(miss))
  for (i in 1:length(miss)){
    q[[i]] = which(is.na(x[miss[i],]))
    p[[i]] = which(!is.na(x[miss[i],]))
  }
  mu = apply(x, 2, mean, na.rm = TRUE)
  sig = pd.sig(var(x, na.rm = TRUE))
  x[which(is.na(x))] = 0
  x.old = 2*x
  count = 0
  while (max(abs(x - x.old)) > tol){
    count = count + 1
    x.old = x
    # expectation
    for (i in 1:length(miss)){
      j = miss[i]
      x[j,q[[i]]] = mu[q[[i]]] + sig[q[[i]],p[[i]]] %*% 
        solve(sig[p[[i]],p[[i]]]) %*% (x[j,p[[i]]] - mu[p[[i]]])
    }
    # maximization
    mu = apply(x, 2, mean)
    sig = var(x)
  }
  cat("Iterations:",count,"\n")
  return (x)
}

#### test on region 2

olive2 <- as.matrix(read.table('~/files/R/666/data/oliver2a.txt',
    stringsAsFactors=FALSE,header=T))
n1 <- nrow(olive2)
p <- ncol(olive2)
max_2 <- em(olive2,tol=1e-6)
mu10 <- c(1300,120,265,7310,820,45,65,28)
S1 <- cov(max_2)
m1 <- apply(max_2,2,mean)
(T2 <- n1 * t(m1 - mu10) %*% solve(S1) %*% (m1 - mu10))
(F <- T2 * ((n1-1) - p + 1) / ((n1-1)*p))

1 - pf(F,p,n1-1 - p +1)

#### test on region 4

olive4 <- as.matrix(read.table('~/files/R/666/data/oliver4a.txt',
    stringsAsFactors=FALSE,header=T))
n2 <- nrow(olive4)
p <- ncol(olive4)

max_4 <- em(olive4,tol=1e-6)
mu20 <- c(1230,105,275,7360,830,41,75,38)
S2 <- cov(max_4)
m2 <- apply(max_4,2,mean)
(T2 <- n2 * t(m2 - mu20) %*% solve(S2) %*% (m2 - mu20))
(F <- T2 * ((n2-1) - p + 1) / ((n2-1)*p))

1 - pf(F,p,n2-1 - p +1)

#### two sample test problem 2

tr <- function(m) sum(diag(m))

Se <- S1/n1 + S2/n2

T2 <- t(m1-m2) %*% solve(Se) %*% (m1-m2)

num <- (tr(Se%*%Se) + (tr(Se))^2)
den <- 1/(n1-1)*( tr(S1/n1)^2 + tr((S1/n1)%*%(S1/n1)) ) + 1/(n2-1)*( tr(S2/n2)^2 + tr((S2/n2)%*%(S2/n2)) )
df.star <- num/den

F <- T2 * ((df.star) - p + 1) / ((df.star)*p)

1 - pf(F,p,df.star - p +1)

a <- solve(Se) %*% (m1 - m2)
a

astar <- diag(sqrt(diag(Se))) %*% a
astar 


### test for missing at random

miss.mat2 <- 1*is.na(olive2)
miss.mat4 <- 1*is.na(olive4)

colnames(miss.mat2) <- colnames(miss.mat4) <- c("a1","a2","a3","a4","a5","a6","a7","a8")

### multinomial test
calc.chi = function(x, n)
    sum((x - n/length(x))^2/(n/length(x)))
miss.mult2 = apply(miss.mat2, 2, sum)
miss.mult4 = apply(miss.mat4, 2, sum)

tt = calc.chi(miss.mult2, sum(miss.mult2))
vec = apply(rmultinom(10000, sum(miss.mult2), prob = rep(1/p, p)),
    2, function(x) calc.chi(x, sum(miss.mult2)))
hist(vec, col='gray'); abline(v=tt)
plot(density(vec), type='l'); abline(v = tt)
mean(vec >= tt)

tt = calc.chi(miss.mult4, sum(miss.mult4))
vec = apply(rmultinom(10000, sum(miss.mult4), prob = rep(1/p, p)),
    2, function(x) calc.chi(x, sum(miss.mult4)))
hist(vec, col='gray'); abline(v=tt)
plot(density(vec), type='l'); abline(v = tt)
mean(vec >= tt)

### test of independence of two subvectors (p.269)
lam.to.f = function(lam, p, nu.h, nu.e){
    w = nu.e + nu.h - 0.5*(p+nu.h+1)
    t = sqrt((p^2 * nu.h^2 - 4)/(p^2 + nu.h^2 - 5))
    df1 = p * nu.h
    df2 = w*t - 0.5*(p*nu.h -2)
    F.stat = (lam^(-1/t) - 1) * df2/df1
#   return (c("F"=F.stat, "crit"=qf(1-0.95, df1, df2, lower.tail = FALSE)))
    return (c("F"=F.stat, "p-val"=pf(F.stat, df1, df2, lower.tail = FALSE)))
    }
des2 <- cbind(max_2,miss.mat2)
des4 <- cbind(max_4,miss.mat4)
des4 = des4[,-13]

S2 = var(des2)
R2 = cor(des2)
(lam2 = det(S2) / (det(S2[1:8,1:8]) * det(S2[9:16,9:16])))
lam.to.f(lam2, 8, 8, n1 - 1 - 8)
S4 = var(des4)
R4 = cor(des4)
(lam4 = det(S4) / (det(S4[1:8,1:8]) * det(S4[9:15,9:15])))
lam.to.f(lam4, 8, 7, n2 - 1 - 7)


m1 <- glm(a1 ~ .,family="binomial",data=as.data.frame(des2[,2:9]))

rbind(predict(m1,newdata=as.data.frame(des2[,2:9]),type="response") > .1, miss.mat2[,1])
