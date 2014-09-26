### 5.12
### part a
dat = read.table("~/files/R/666/data/T3_6_PROBE.DAT", row.names = 1)

y = as.matrix(dat)
colnames(y) <- rownames(y) <- NULL

n = nrow(y)
p = ncol(y)
nu = n - 1
ybar = apply(y, 2, mean)
S = var(y)
Sinv = solve(S)

# hypothesize mean vector
mu0 = c(30, 25, 40, 25, 30)

# Hotellings T^2
T2 = n * t(ybar - mu0) %*% Sinv %*% (ybar - mu0)

# convert to F-statistic, distributed F(p, nu-p+1)
F.stat = (nu - p + 1)/(nu * p) * T2

# p-value
pf(F.stat, p, nu - p + 1, lower.tail = FALSE)
# low p-value, so we reject

### part b
t.stat = (ybar - mu0) / sqrt(diag(S)/n)
# p-value
2*pt(abs(t.stat), nu, lower.tail = FALSE)
2*pt(abs(t.stat), nu, lower.tail = FALSE) < 0.05
# reject mu1 and mu3
