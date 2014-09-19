# MANOVA chapter 6
# (univariate)
k = 5
n = 5
p = 1
y = matrix(rnorm(n*k*p), k*n, p)
 y[seq(n*(3-1)+1, n*3),] = rnorm(n*p, 1)

ybar = matrix(0, k, p)
ydot = matrix(0, k, p)
for (i in 1:k){
    for (j in 1:p){
        ybar[i,j] = mean(y[seq(n*(i-1)+1, n*i),j])
        ydot[i,j] = sum(y[seq(n*(i-1)+1, n*i),j])
        }
    }
# the ith row in ydot is y_i_dot, the sum of the
# observations from sample i

SSH = sum(ydot^2)/n - (sum(ydot))^2/(k*n)
SSE = sum(y^2) - sum(ydot^2)/n

# using matrix multiplcation
one.k = matrix(1, k, 1)
SSH = (t(ydot) %*% ydot)/n - (t(one.k) %*% ydot %*% t(ydot) %*% one.k) / (k*n)
SSE = t(y) %*% y - (t(ydot) %*% ydot)/n

F.val = (SSH/(k-1)) / (SSE/(k*(n-1)))
# F.val is distributed F(k-1, k*(n-1))
alpha = 0.05
F.star = qf(alpha, k-1, k*(n-1), lower.tail=FALSE)

# reject the null that all means are equivalent when
# F.val is greater than F.star
data.frame(F.val, F.star, "Reject"=F.val > F.star)


# (multivariate) (using data from Example 6.1.7, page 183)
y = read.table("~/files/R/666/data/T6_2_ROOT.DAT")
k = max(y[,1])
n = nrow(y) / k
p = ncol(y) - 1
y = as.matrix(y[,-1])

# remove the V1, V2, ... tags
colnames(y) = NULL

ydot = matrix(0, k, p)
for (i in 1:k)
    for (j in 1:p)
        ydot[i,j] = sum(y[seq(n*(i-1)+1, n*i),j])
# the ith row in ydot is y_i_dot, the sum of the
# observations from sample i

# y_dot_dot
dd = apply(ydot, 2, sum)

# H and E are p x p
# H is the "between" sum of squares on the diagonal
# for each of the p variables, off-diagonals are
# analogous sums of products for each pair of variables
H = (t(ydot) %*% ydot)/n - (dd %*% t(dd)) / (k*n)
# E is "within" sum of squares for each variable
E = t(y) %*% y - (t(ydot) %*% ydot)/n
# the expected value of E = k*(n-1) Sigma

# degrees of freedom
nu_H = k - 1
nu_E = k*(n-1)

# Wilk's big lambda test, 0 <= Lambda <= 1, reject for
# small values of lambda
# eq 6.13
(lambda = det(E) / det(E + H))
system.time(for (i in 1:1000) det(E) / det(E + H))

# eq 6.14 (both lambdas are equivalent)
(lambda = prod(1/(1+eigen(solve(E) %*% H)$values)))
# H %*% solve(E) could also be used instead

# eq 6.15 (approximate F test)
w = nu_E + nu_H - 0.5*(p + nu_H + 1)
t = sqrt((p^2 * nu_H^2 - 4)/(p^2 + nu_H^2-5))
df1 = p * nu_H
df2 = w*t - 0.5*(p * nu_H -2)

F.val = (1-lambda^(1/t))/lambda^(1/t) * df2/df1
alpha = 0.05
F.star = qf(alpha, df1, df2, lower.tail=FALSE)

data.frame(F.val, F.star, "Reject" = F.val > F.star)

# eq 6.16 (approximate Chi-sq test)
chi.val = -(nu_E - 0.5*(p - nu_H + 1)) * log(lambda)
chi.star = qchisq(alpha, p * nu_H)
data.frame(chi.val, chi.star, "Reject" = chi.val > chi.star)
