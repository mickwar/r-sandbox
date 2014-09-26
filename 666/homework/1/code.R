### 3.11
dat = read.table("~/files/R/666/data/T3_4_CALCIUM.DAT")
Y = as.matrix(dat[,-1])

n = nrow(Y)
J = matrix(1, n, n)

S = 1/(n-1) * t(Y) %*% (diag(n) - J/n) %*% Y
# cov(Y) is equivalent

# generalized sample variance |S|
det(S)

# total sample variance tr(S)
sum(diag(S))

### 3.15
z = matrix(c(3,-1,2),3,1)
w = matrix(c(-2,3,1),3,1)

Z = Y %*% z
W = Y %*% w

# r_xy (eq 3.13)
sum((Z - mean(Z))*(W - mean(W))) / 
    sqrt( sum( (Z - mean(Z))^2) * sum( (W - mean(W))^2))

# r_xy (eq 3.57)
t(z) %*% S %*% w / sqrt( (t(z) %*% S %*% z) * (t(w) %*% S %*% w))

### 3.17
A = matrix(c(1,2,-1,1,-3,-2,1,2,-3),3,3)
ybar = apply(Y, 2, mean)
zbar = A %*% ybar

Sz = A %*% S %*% t(A)

Dz = diag(sqrt(diag(Sz)))
Rz = solve(Dz) %*% Sz %*% solve(Dz)

Dz %*% Rz %*% Dz

### 3.22
dat = read.table("~/files/R/666/data/T3_9_GLUCOSE.DAT")
names(dat) = c("y1", "y2", "y3", "x1", "x2", "x3")

bar = apply(dat, 2, mean)

bar
cov(dat)

### 4.1
sig1 = matrix(c(14,8,3,8,5,2,3,2,1),3,3)
sig2 = matrix(c(6,6,1,6,8,2,1,2,1),3,3)

14*5*1 + 8*2*3 + 3*8*2 - 14*2*2 - 8*8*1 - 3*5*3
det(sig1)

6*8*1 + 6*2*1 + 1*6*2 - 6*2*2 - 6*6*1 -1*8*1
det(sig2)

d1 = diag(sqrt(diag(sig1)))
d2 = diag(sqrt(diag(sig2)))

solve(d1) %*% sig1 %*% solve(d1)
solve(d2) %*% sig2 %*% solve(d2)

### 4.10
A = matrix(c(1,1,1,-1,1,2), 2,3)
mu = matrix(c(3,1,4),3,1)
sigma = matrix(c(6,1,-2,1,13,4,-2,4,4),3,3)

A %*% mu
A %*% sigma %*% t(A)

D = matrix(c(1,0,0.5, 0, 0, 0.5, 0, 1, 0),3,3)
D %*% mu
D %*% sigma %*% t(D)

### 4.17
mu.y = matrix(c(3,-2),2,1)
mu.x = matrix(c(4,-3,5),3,1)
sig.yy = matrix(c(14,-8,-8,18),2,2)
sig.yx = matrix(c(15,8,0,6,3,-2),2,3)
sig.xx = matrix(c(50,8,5,8,4,0,5,0,1),3,3)

# E(y|x), unspecified x
mu.y - sig.yx %*% solve(sig.xx) %*% mu.x
sig.yx %*% solve(sig.xx)

# cov(y|x)
sig.yy - sig.yx %*% solve(sig.xx) %*% t(sig.yx)




### last problem
set.seed(2)
n = 10
p = 5
# total data matrix
(Y = matrix(rnorm((n+1)*p), n+1, p))

# yesterday's data matrix
X = Y[1:n,]
# new observation
(z = matrix(Y[n+1,], p, 1))

# current covariance matrix
(X.cov = cov(X))
# current inverse
(X.inv = solve(X.cov))

Y.cov = cov(Y)
# target
(Y.inv = solve(Y.cov))

# calcule the c in B + cc'
(c = sqrt(n+1)*(z - apply(Y, 2, mean)) / n)

# note the equality
(n-1)/n*X.cov + c %*% t(c)
Y.cov

B.inv = n/(n-1)*X.inv

# special case of sherman-morrison-woodbury
B.inv - (B.inv %*% c %*% t(c) %*% B.inv) / as.vector(1 + t(c) %*% B.inv %*% c)
# is equal to todays inverse (used a check)
Y.inv




# add new data point
Y = rbind(Y, rnorm(p))
n = n + 1

z = Y[n+1,]
z = sqrt(n+1)*(z - apply(Y, 2, mean)) / n

new.inv = n/(n-1)*Updated
Updated = new.inv - (new.inv %*% z %*% t(z) %*% new.inv) / as.vector(1 + t(z) %*% new.inv %*% z)
solve(cov(Y))
###########
