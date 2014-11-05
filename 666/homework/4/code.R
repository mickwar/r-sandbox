rm(list=ls())
# F approximation for wilks lambda distribution
lam.to.f = function(lam, p, nu.h, nu.e){
    w = nu.e + nu.h - 0.5*(p+nu.h+1)
    t = sqrt((p^2 * nu.h^2 - 4)/(p^2 + nu.h^2 - 5))
    df1 = p * nu.h
    df2 = w*t - 0.5*(p*nu.h -2)
    F.stat = (lam^(-1/t) - 1) * df2/df1
#   return (c("F"=F.stat, "crit"=qf(1-0.95, df1, df2, lower.tail = FALSE)))
    return (c("F"=F.stat, "p-val"=pf(F.stat, df1, df2, lower.tail = FALSE)))
    }

##### Part 1
X = read.table("~/files/R/666/data/PSAData.txt", header = TRUE)
xname = c("Intercept", names(X))
X = as.matrix(X)
# use intercept? getting problems in calculating Sxx, Sxy, etc.
# Yes, use the intercept in the regression, just don't include it
# when doing the canonical correlation analysis
X.1 = as.matrix(cbind(1, X))
colnames(X.1) = xname

Y = read.table("~/files/R/666/data/PSAContributions.txt", header = TRUE)
yname = names(Y)
Y = as.matrix(Y)

p = ncol(Y)
q = ncol(X)
n = nrow(X)
ybar = apply(Y, 2, mean)

betahat = solve(t(X.1) %*% X.1) %*% t(X.1) %*% Y
E = t(Y - X.1 %*% betahat) %*% (Y - X.1 %*% betahat)
H = t(betahat) %*% t(X.1) %*% Y - n * ybar %*% t(ybar)
# the n * ybar %*% t(ybar) term is to subtract out the portion
# due to the intercept

### using covariance (because why ever?)
Sxx = var(X)
Syy = var(Y)
Sxy = var(X, Y)

Sa = solve(Syy) %*% t(Sxy) %*% solve(Sxx) %*% Sxy
Sb = solve(Sxx) %*% Sxy %*% solve(Syy) %*% t(Sxy)

a = eigen(Sa)$vector
c = diag(sqrt(diag(Syy))) %*% a
c[,1:2]

### using correlation
R = cor(cbind(X, Y))
Rxx = cor(X)
Ryy = cor(Y)
Rxy = cor(X, Y)

# related to the y's
Ra = solve(Ryy) %*% t(Rxy) %*% solve(Rxx) %*% Rxy
 solve(Ryy, t(Rxy)) %*% solve(Rxx,Rxy)
# related to the x's
Rb = solve(Rxx) %*% Rxy %*% solve(Ryy) %*% t(Rxy)


# 1: is there a significant relationship between x's and y's
lam.to.f(det(R) / (det(Rxx) * det(Ryy)), p, q, n-1-q)
# reject that there is no linear relationship between x's and y's

# 2: what is the essential dimensionality of the relationship
#    between x's and y's
evals = eigen(solve(E) %*% H)$values
evals / sum(evals) # first 2 eigen values

eigen(Ra)$vectors[,1:2]/c[,1:2] # y
c[,1:2]

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

nu.h = h
nu.e = n - q - 1

# since nu.h = 1, there is an exact F transformation
F.stat = (1 - lam) / lam * (nu.e - p + 1) / p
pf(F.stat, p, nu.e - p + 1, lower.tail = FALSE)
# crazy small, reject that snitch (i.e., Pb add much)
lam.to.f(lam, p, nu.h, nu.e)

##### Part 2
dat = matrix(scan("~/files/R/666/data/bodyfat.txt"),
    252, 15, byrow = TRUE)
yname = c("density", "bodyfat")
xname = c("age", "weight", "neck", "chest", "abdomen", "hip",
    "thigh", "knee", "ankle", "biceps", "forearm", "wrist")

Y = dat[,1:2]
X = dat[,3:15]
X.1 = cbind(1, X)

p = ncol(Y)
q = ncol(X.1)
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


# step function to work for mlm models, beginning with
# only the intercept, ending with all column given in X
mw.step = function(start.index, y, x, method = "both",
    FUN = lm, k = 2, max.step = 1000, output = TRUE){
    if (missing(start.index))
        start.index = 1
    x = as.matrix(x)
    if (!all(x[,1] == 1))
        x = cbind(1, x)
    y = as.matrix(y)
    n = nrow(x)
    q = ncol(x)
    p = ncol(y)

    var.in = start.index
    var.out = (1:q)[-var.in]

    mod = FUN(y ~ 0 + x[,var.in])
    aic = extractAIC(mod, k = k)[2L]
    test.aic = rep(Inf, q)
    if (output)
        cat("\n    Starting AIC:", aic, "\n")

    count = 0
    while (count < max.step){
        flag = 0
        count = count + 1
        if (output)
            cat("\n    Iteration:", count, "\n")
        # for forward selection
        if (method == "forward" || method == "both"){
            for (check in var.out){
                test.aic[check] = extractAIC(FUN(y ~ 0 + x[,c(var.in, check)]), k = k)[2L]
                }
            if (min(test.aic) < aic){
                var.out = var.out[-which(var.out == which.min(test.aic))]
                var.in = c(var.in, which.min(test.aic))
                aic = min(test.aic)
                if (output){
                    cat("        AIC:", aic, "\n")
                    cat("        Add variable:", which.min(test.aic), "\n")
                    }
            } else {
                flag = flag + 1
                if (output){
                    cat("        AIC:", aic, "\n")
                    cat("        No variable added\n")
                    }
                }
            test.aic = rep(Inf, q)
            }
        # backward selection
        if (method == "backward" || method == "both"){
            temp.in = var.in
            # don't check the most recently added variable (flag = 1 if none removed)
            if (method == "both" && flag == 0)
                temp.in = temp.in[-length(temp.in)]
            for (check in temp.in){
                test.aic[check] = extractAIC(FUN(y ~ 0 + x[,var.in[-which(check == temp.in)]]), k = k)[2L]
                }
            if (min(test.aic) < aic){
                var.in = var.in[-which(var.in == which.min(test.aic))]
                var.out = c(var.out, which.min(test.aic))
                aic = min(test.aic)
                if (output){
                    cat("        AIC:", aic, "\n")
                    cat("        Remove variable:", which.min(test.aic), "\n")
                    }
            } else {
                flag = flag + 1
                if (output){
                    cat("        AIC:", aic, "\n")
                    cat("        No variable removed\n")
                    }
                }
            }
        if ((method == "forward" || method == "backward") && flag == 1)
            count = max.step
        if (method == "both" && flag == 2)
            count = max.step
        }
    var.in
#    FUN(y ~ 0 + x[,var.in])
    }

# the index of the variables kept after backward selecting 
index = sort(mw.step(1:q, Y, X.1, k = log(n), method = "backward"))

### using the reduced set of variables
X.1 = X.1[,index]
X = X.1[,-1]

p = ncol(Y)
q = ncol(X.1)
n = nrow(X)
ybar = apply(Y, 2, mean)

betahat = solve(t(X.1) %*% X.1) %*% t(X.1) %*% Y
E = t(Y - X.1 %*% betahat) %*% (Y - X.1 %*% betahat)
H = t(betahat) %*% t(X.1) %*% Y - n * ybar %*% t(ybar)

evals = eigen(solve(E) %*% H)$values
evals / sum(evals) # only the first eigen value

# using correlation
Rxx = cor(X)
Ryy = cor(Y)
Rxy = cor(X, Y)

# related to the y's
Ra = solve(Ryy) %*% t(Rxy) %*% solve(Rxx) %*% Rxy
# related to the x's
Rb = solve(Rxx) %*% Rxy %*% solve(Ryy) %*% t(Rxy)

eigen(Ra)$vectors[,1] # y
as.numeric(eigen(Rb)$vectors[,1]) # x
