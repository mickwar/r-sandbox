### 5.12
### part a
dat = read.table("~/files/R/666/data/T3_6_PROBE.DAT", row.names = 1)

y = as.matrix(dat)
colnames(y) = paste("X", 1:5, sep="")
rownames(y) = NULL

n = nrow(y)
p = ncol(y)
nu = n - 1
ybar = apply(y, 2, mean)
S = var(y)
Sinv = solve(S)

# hypothesize mean vector
mu0 = c(30, 25, 40, 25, 30)

# Hotellings T^2
(T2 = n * t(ybar - mu0) %*% Sinv %*% (ybar - mu0))

# convert to F-statistic, distributed F(p, nu-p+1)
F.stat = (nu - p + 1)/(nu * p) * T2

# p-value
pf(F.stat, p, nu - p + 1, lower.tail = FALSE)
# low p-value, so we reject

### part b
t.stat = (ybar - mu0) / sqrt(diag(S)/n)
# p-value
data.frame("t.stat"=t.stat, "p-val"=2*pt(abs(t.stat), nu, lower.tail = FALSE),
"Reject"=2*pt(abs(t.stat), nu, lower.tail = FALSE) < 0.05)
# reject mu1 and mu3

### part c
# discriminant function
a = Sinv %*% (ybar - mu0)
D = diag(sqrt(diag(S)))
(a.star = D %*% a)

### 5.16
# two sample comparison
dat = read.table("~/files/R/666/data/T5_5_FBEETLES.DAT")

y = as.matrix(dat[,-1])
colnames(y) = c("sample", paste("X", 1:4, sep=""))
rownames(y) = NULL

y1 = y[y[,1]==1,-1]
y2 = y[y[,1]==2,-1]

n = c(nrow(y1), nrow(y2))
p = ncol(y1)
nu = sum(n)-2
ybar = rbind(apply(y1, 2, mean), apply(y2, 2, mean))
S = list(var(y1), var(y2))
Sinv = list(solve(S[[1]]), solve(S[[2]]))

Spl = ((n[1] - 1) * S[[1]] + (n[2] - 1) * S[[2]]) / (sum(n) - 2)
Spl.inv = solve(Spl)

# Hotellings T^2
(T2 = (n[1]*n[2])/(n[1]+n[2]) * t(ybar[1,] - ybar[2,]) %*% Spl.inv %*% (ybar[1,] - ybar[2,]))

F.stat = (nu - p + 1)/(nu*p) * T2
pf(F.stat, p, nu-p+1, lower.tail = FALSE)

# univariate tests
(t.stat = sqrt(n[1]*n[2]/(n[1]+n[2])) * (ybar[1,] - ybar[2,]) / sqrt(diag(Spl)))

# discriminant function
a = Spl.inv %*% (ybar[1,] - ybar[2,])
a.star = diag(sqrt(diag(Spl))) %*% a

### 5.20
dat = read.table("~/files/R/666/data/T5_8_GOODS.DAT")
y = as.matrix(dat[,-1])
colnames(y) = c("sample", paste("X", 1:4, sep=""))
rownames(y) = NULL

y1 = y[y[,1]==1,-1]
y2 = y[y[,1]==2,-1]

n = c(nrow(y1), nrow(y2))
p = ncol(y1)
nu = sum(n)-2
ybar = rbind(apply(y1, 2, mean), apply(y2, 2, mean))
S = list(var(y1), var(y2))
Sinv = list(solve(S[[1]]), solve(S[[2]]))

Spl = ((n[1] - 1) * S[[1]] + (n[2] - 1) * S[[2]]) / (sum(n) - 2)
Spl.inv = solve(Spl)

# Hotellings T^2
(T2 = (n[1]*n[2])/(n[1]+n[2]) * t(ybar[1,] - ybar[2,]) %*% Spl.inv %*% (ybar[1,] - ybar[2,]))

(F.stat = (nu - p + 1)/(nu*p) * T2)
pf(F.stat, p, nu-p+1, lower.tail = FALSE)

# Nel and Van der Merwe
# compare the standard deviations from each covariance matrix
rbind(sqrt(diag(S[[1]])), sqrt(diag(S[[2]])))

S.e = S[[1]]/n[1] + S[[2]]/n[2]

(T.star = t(ybar[1,] - ybar[2,]) %*% solve(S.e) %*% (ybar[1,] - ybar[2,]))

tr = function(x)
    sum(diag(x))

# what is S.e^2? is it S.e %*% S.e?
den = 0
for (i in 1:length(n))
    den = den + 1/(n[i] - 1) * (tr(S[[i]] %*% S[[i]] / n[i]^2) + tr(S[[i]]/n[i])^2)
(nu.star = (tr(S.e %*% S.e) + tr(S.e)^2) / den)

(F.stat = (nu.star - p + 1)/(nu.star*p) * T.star)
pf(F.stat, p, nu.star-p+1, lower.tail = FALSE)

### 5.21
dat = read.table("~/files/R/666/data/T5_9_ESSAY.DAT")
y = as.matrix(dat[,-1])
colnames(y) = c("y1", "y2", "x1", "x2")
rownames(y) = NULL

d = cbind(y[,1]-y[,3], y[,2]-y[,4])
dbar = apply(d, 2, mean)
S = var(d)
Sinv = solve(S)
n = nrow(d)
p = ncol(d)

(T2 = n * t(dbar) %*% Sinv %*% dbar)
(F.stat = (n - p)/((n-1)*p) * T2)
pf(F.stat, p, n - p, lower.tail=FALSE)

(a = Sinv %*% dbar)
(a.star = diag(sqrt(diag(S))) %*% a)

(t.stat = dbar/sqrt(diag(S)/n))
pt(2*abs(t.stat), n-1, lower.tail = FALSE)

### problem 6
mu = c(3, 11, -4, 1)
sigma = matrix(c(11,3,3,1,3,10,2,-5,3,2,10,-6,1,-5,-6,15),4,4)

a = c(3,-2, 0, 0)
t(a) %*% mu
t(a) %*% sigma %*% a

A = matrix(c(3,0,0,-2,0,0,0,1,0,0,0,1),3,4)
(new.mu = A %*% mu)
(new.sig = A %*% sigma %*% t(A))

new.mu[1] + new.sig[1,2:3] %*% solve(new.sig[2:3,2:3]) %*% (c(1,-2)-new.mu[2:3])
new.sig[1,1] - new.sig[1,2:3] %*% solve(new.sig[2:3,2:3]) %*% new.sig[2:3,1]

det(sigma)
