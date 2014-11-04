##### Part 1
X = read.table("~/files/R/666/data/PSAData.txt", header = TRUE)
xname = c("Intercept", names(X))
X = as.matrix(X)
# use intercept? getting problems in calculating Sxx, Sxy, etc.
X.1 = as.matrix(cbind(1, X))

Y = read.table("~/files/R/666/data/PSAContributions.txt", header = TRUE)
yname = names(Y)
Y = as.matrix(Y)

n = nrow(X)
ybar = apply(Y, 2, mean)

betahat = solve(t(X.1) %*% X.1) %*% t(X.1) %*% Y
E = t(Y - X.1 %*% betahat) %*% (Y - X.1 %*% betahat)
H = t(betahat) %*% t(X.1) %*% Y - n * ybar %*% t(ybar)

evals = eigen(solve(E) %*% H)$values
evals / sum(evals) # first 2 eigen values

### using covariance (because why ever?)
# Sxx = var(X)
# Syy = var(Y)
# Sxy = var(X, Y)
# 
# Sa = solve(Syy) %*% t(Sxy) %*% solve(Sxx) %*% Sxy
# Sb = solve(Sxx) %*% Sxy %*% solve(Syy) %*% t(Sxy)

### using correlation
Rxx = cor(X)
Ryy = cor(Y)
Rxy = cor(X, Y)

# related to the y's
Ra = solve(Ryy) %*% t(Rxy) %*% solve(Rxx) %*% Rxy
# related to the x's
Rb = solve(Rxx) %*% Rxy %*% solve(Ryy) %*% t(Rxy)

eigen(Ra)$vectors[,1:2] # y
eigen(Rb)$vectors[,1:2] # x
# for eigen 1: small x6 is related to large y1 and y5
# for eigen 2: large x4 and small x6 corresponds with large y1 and small y5

### remove one of the variables
# testing whether the column in the beta matrix corresponding
# to the variable is all 0's.
Xr = X.1[,-which(xname == "Pb")]
Br = solve(t(Xr) %*% Xr) %*% t(Xr) %*% Y
Hr = t(betahat) %*% t(X.1) %*% Y - t(Br) %*% t(Xr) %*% Y

# wilks lambda for comparing the full and reduced models
lam = det(E) / det(E + Hr) # distributed Lambda(p,h,n-q-1)
h = 1 # number of remove variables
p = ncol(Y)
q = ncol(X)

nu.h = h
nu.e = n - q - 1

# since nu.h = 1, there is an exact F transformation
F.stat = (1 - lam) / lam * (nu.e - p + 1) / p
pf(F.stat, p, nu.e - p + 1, lower.tail = FALSE)
# crazy small, reject that snitch (i.e., Pb add much)

##### Part 2
dat = matrix(scan("~/files/R/666/data/bodyfat.txt"),
    252, 15, byrow = TRUE)
yname = c("density", "bodyfat")
xname = c("age", "weight", "neck", "chest", "abdomen", "hip",
    "thigh", "knee", "ankle", "biceps", "forearm", "wrist")

Y = dat[,1:2]
X = dat[,3:15]
X.1 = cbind(1, X)

n = nrow(X)
ybar = apply(Y, 2, mean)

betahat = solve(t(X.1) %*% X.1) %*% t(X.1) %*% Y
E = t(Y - X.1 %*% betahat) %*% (Y - X.1 %*% betahat)
H = t(betahat) %*% t(X.1) %*% Y - n * ybar %*% t(ybar)

evals = eigen(solve(E) %*% H)$values
evals / sum(evals) # only the first eigen value

### using correlation
Rxx = cor(X)
Ryy = cor(Y)
Rxy = cor(X, Y)

# related to the y's
Ra = solve(Ryy) %*% t(Rxy) %*% solve(Rxx) %*% Rxy
# related to the x's
Rb = solve(Rxx) %*% Rxy %*% solve(Ryy) %*% t(Rxy)

eigen(Ra)$vectors[,1] # y
as.numeric(eigen(Rb)$vectors[,1]) # x
