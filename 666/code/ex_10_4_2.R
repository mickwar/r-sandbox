dat = read.table("~/files/R/666/data/T10_1_CHEM.DAT")

y = as.matrix(dat[,2:4])
x = as.matrix(cbind(1, dat[,5:7]))

# multivariate
bhat = solve(t(x) %*% x) %*% t(x) %*% y

# univariate on each y
bhat1 = solve(t(x) %*% x) %*% t(x) %*% y[,1]
bhat2 = solve(t(x) %*% x) %*% t(x) %*% y[,2]
bhat3 = solve(t(x) %*% x) %*% t(x) %*% y[,3]

# equivalence (when should they not be?)
cbind(bhat1, bhat2, bhat3)
bhat
