### Frequentist LASSO example
library(glmnet)

### Randomize a data set
set.seed(1)
n = 50
p = 100
X = matrix(rnorm(n*p), n, p)
beta = c(rep(3, round(p/4)), rep(1, round(p/4)),
    rep(-1, round(p/4)), rep(0, round(p/4)))
y = X %*% beta + rnorm(n)

### Get lambda via cross-validation (10-fold)
lambda = cv.glmnet(X, y, nfolds = 10)$lambda.min

### Fit model, get beta coefficients, get predictions
mod = glmnet(X, y, family = "gaussian", lambda = lambda, alpha = 1)
beta.hat = coef(mod)
y.pred = predict(mod, newx = X)

### Observed vs. fitted plot
plot(y, y.pred)
abline(0, 1)

### Compare beta.hat to truth
plot(beta.hat, pch = 16)
points(beta, col = 'dodgerblue', pch = 15)

