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
n = 11
p = 5
Y = matrix(runif(n*p), n, p)

D = matrix(sqrt(diag(t(Y) %*% Y - t(Y[1:10,]) %*% Y[1:10,])), p, 1)

t(Y) %*% Y - t(Y[1:10,]) %*% Y[1:10,]
D %*% t(D) 

t(Y[1:10,]) %*% Y[1:10,] + D %*% t(D)
t(Y) %*% Y

Ybar.n = apply(Y, 2, mean)
Ybar.10 = apply(Y[1:10,], 2, mean)

# cov(Y)
(t(Y) %*% Y - n * Ybar.n %*% t(Ybar.n))/(n-1)

# given:
n*Ybar.n%*%t(Ybar.n) - (n-1) * Ybar.10 %*% t(Ybar.10)

H = (n-1)*Ybar.10 %*% t(Ybar.10)
squig = matrix(sqrt(diag(H/n)), p, 1)
(squig + D) %*% (t(squig + D)
 * ((n-1)^2)/(n^2)
n * Ybar.n %*% t(Ybar.n)


S = cov(Y)
S.10 = cov(Y[1:10,])

J = rep(1, n-1)
1/(n-1) * (t(Y[1:10,]) %*% Y[1:10,] + D %*% t(D) - 1/n * (t(Y[1:10,]) %*% J + D) %*%
    (t(J) %*% Y[1:10,] + t(D)))

(S.10 * (n-2) + D %*% t(D))/(n-1)
S

t(Y[1:10,]) %*% J %*% t(D)
t(D %*% t(J) %*% Y[1:10,])
