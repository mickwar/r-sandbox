### simplified expectation maximization (for values missing at random)
# "generate" missing values
gen.miss = function(x, m = 0.01){
    # m is proportion of missing values
    n = nrow(x)
    p = ncol(x)
    x[sample(n*p, floor(n*p*m))] = NA
    return (x)
    }
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


# no correlation among variables
set.seed(1)
n = 100
p = 5
x1 = matrix(rnorm(n*p), n, p)
x1.miss = gen.miss(x1, 0.10)

y1 = em(x1.miss, tol=0.00001)

x1.miss
max(abs(y1 - x1))

# correlation existing
set.seed(1)
n = 20000
p = 10
# uncorrelated
A = diag(p)
# weakly correlated
A = chol(matrix(c(1.00,0.50,0.25,0.12,0.06,
                  0.50,1.00,0.50,0.25,0.12,
                  0.25,0.50,1.00,0.50,0.25,
                  0.12,0.25,0.50,1.00,0.50,
                  0.06,0.12,0.25,0.50,1.00), 5, 5))
# strongly correlated
A = chol(matrix(c(1.00,0.98,0.96,0.94,0.92,
                  0.98,1.00,0.98,0.96,0.94,
                  0.96,0.98,1.00,0.98,0.96,
                  0.94,0.96,0.98,1.00,0.98,
                  0.92,0.94,0.96,0.98,1.00), 5, 5))
# more correlated
A = chol(matrix(c(1.00,0.99,0.99,0.99,0.99,
                  0.99,1.00,0.99,0.99,0.99,
                  0.99,0.99,1.00,0.99,0.99,
                  0.99,0.99,0.99,1.00,0.99,
                  0.99,0.99,0.99,0.99,1.00), 5, 5))

x2 = matrix(rnorm(n*p), n, p) %*% A
x2.miss = gen.miss(x2, 0.20)

y2 = em(x2.miss, tol=0.00001)

missing = which(is.na(x2.miss))
plot(density(y2[missing] - x2[missing]))
curve(dnorm(x, mean(y2[missing]-x2[missing]),
    sd(y2[missing]-x2[missing])), col='red', add=TRUE)
    
range(y2[missing]-x2[missing])
quantile(y2[missing]-x2[missing], c(0.025, 0.975))

###### test for missing at random
### 
dat = read.table("~/files/R/666/data/oliver4a.txt", header=TRUE)
n.miss = sum(is.na(dat))
zero = double(nrow(dat))
for (i in 1:length(zero))
    zero[i] = sum(is.na(dat[i,]))

gen.miss = function(x, k){
    # k is number of missing values
    n = nrow(x)
    p = ncol(x)
    # 1 indicates missing (for zero matrix)
    x[sample(n*p, k)] = 1
    return (x)
    }

M = 1000
# 9 is for 0 through 8, the number of possible missing acids
out = matrix(0, M, ncol(dat)+1)
for (i in 1:M){
    x = gen.miss(matrix(0, nrow(dat), ncol(dat)), n.miss)
    x.miss = apply(x, 1, sum)
    out[i,1:length(table(x.miss))] = table(x.miss)
    }

plot(table(out[,1])/M)
abline(v=(0.1+table(zero)[1]), col='red')
plot(table(out[,2])/M)
abline(v=(0.1+table(zero)[2]), col='red')
plot(table(out[,3])/M)
abline(v=(0.1+table(zero)[3]), col='red')
plot(table(out[,4])/M)
abline(v=(0.1+table(zero)[4]), col='red')
plot(table(out[,5])/M)
abline(v=(0.1), col='red')

pairs(rbind(out[,1:5], c(table(zero), 0, 0)), pch=20,
    col=c(rep("black", M), "red"))
