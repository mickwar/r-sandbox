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
