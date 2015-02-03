logit = function(x)
    log(x / (1-x))
logistic = function(x)
    exp(x) / (1 + exp(x))

li = c(8, 8, 10, 10, 12, 12, 12, 14, 14, 14, 16, 16, 16, 18, 20, 20, 20, 22, 22, 24, 26, 28, 32, 34, 38, 38, 38)
y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0)
cbind(y, li)


mod = glm(y ~ li, family = binomial(link = "logit"))
summary(mod)
exp(coef(mod))

phat = logistic(cbind(1, li) %*% coef(mod))

p0 = seq(1, 0, length = 101)
sens = double(length(p0))
spec = double(length(p0))
for (i in 1:length(p0)){
    yhat = ifelse(phat >= p0[i], 1, 0)
    sens[i] = sum(y * yhat) / sum(y)
    spec[i] = sum((1-y) * (1-yhat)) / sum(1-y)
    }
cbind(1-spec, sens, p0)

calc_auc = function(sens, spec){
    y = sens
    x = 1 - spec
    auc = double(length(y)-1)
    # this only works for certain functions
    for (i in 1:(length(y)-1))
        auc[i] = y[i]*(x[i+1]-x[i]) + 0.5*(y[i+1]-y[i])*(x[i+1]-x[i])
    return (sum(auc))
    }
calc_auc(sens, spec)

plot(1-spec, sens, type='l', main = paste0("ROC Curve -- AUC = ",
    round(calc_auc(sens, spec),3)),
    xlab = "1 - Specificity", ylab = "Sensitivity")
lines(c(0,1), c(0,1), lty=2, col='gray', lwd=2)
