calc_auc = function(sens, spec){
    y = sens
    x = 1 - spec
    auc = double(length(y)-1)
    # this only works for certain functions
    for (i in 1:(length(y)-1))
        auc[i] = y[i]*(x[i+1]-x[i]) + 0.5*(y[i+1]-y[i])*(x[i+1]-x[i])
    return (sum(auc))
    }


dat = read.table("~/files/data/health.txt", header=TRUE)

y = ifelse(dat$sick == "none", 0, 1)
x = dat[,2:ncol(dat)]
x = x[,-2]

int.mod = glm(y ~ 1, data = x, family = binomial)
full.mod = glm(y ~ ., data = x, family = binomial)

mod = step(int.mod, scope = list("lower"=int.mod, "upper"=full.mod))
summary(mod)
# protein is getting included for some reason even
# though it has massive standard error

mod = glm(y ~ 1 + dairy + egg + cheese + meat + vegetable,
    data = x, family = binomial)
summary(mod)

phat = predict(mod, type="response")

p = seq(1, 0, length = 101)
sens = double(length(p))
spec = double(length(p))
for (i in 1:length(p)){
    yhat = ifelse(phat >= p[i], 1, 0)
    sens[i] = sum(y * yhat) / sum(y)
    spec[i] = sum((1-y) * (1-yhat)) / sum(1-y)
    }

### sensitivty at p0 = 0.5 (p0[51] = 0.5)
cfix = 83
yhat = ifelse(phat >= p[cfix], 1, 0)
# sens
c(sum(y * yhat), sum(y))
sum(y * yhat) / sum(y)
# spec
c(sum((1-y) * (1-yhat)), sum(1-y))
sum((1-y) * (1-yhat)) / sum(1-y)

calc_auc(sens, spec)

plot(1-spec, sens, type='l', main = paste0("ROC Curve -- AUC = ",
    round(calc_auc(sens, spec),3)), lwd = 2,
    xlab = "1 - Specificity", ylab = "Sensitivity")
points(1-spec[cfix], sens[cfix], col='red', lwd=3, pch=20)
lines(c(0,1), c(0,1), lty=2, col='darkgray', lwd=2)


n = length(y)
m = 101
lower = min(p)
upper = max(p)
cutoff = seq(lower,upper,length=m)
pos.rate = double(m)
neg.rate = double(m)
error = double(m)

train.n = 0

for (i in 1:m){
    t.pos = sum(y == 1 & phat > cutoff[i])
    t.neg = sum(y == 0 & phat < cutoff[i])
    total.neg = (n-train.n)-sum(y)
    total.pos = sum(y)

    pos.rate[i] = 1-t.pos/total.pos
    neg.rate[i] = 1-t.neg/total.neg
    error[i] = 1-(t.pos+t.neg)/(n-train.n)
    }

plot(cutoff, pos.rate, type='s', col="darkblue", ylab="Error Rate",
    xlab="Threshold", main="Predictive Error Rates")
points(cutoff, neg.rate, type='s', col="red")
points(cutoff, error, type='s', col="green", lwd=2)
points(cutoff[m-cfix], pos.rate[m-cfix], pch=20, col="darkblue", lwd=5)
points(cutoff[m-cfix], neg.rate[m-cfix], pch=20, col="red", lwd=5)
points(cutoff[m-cfix], error[m-cfix], pch=20, col="green", lwd=5)
legend(0.3, 1, legend=c("False Positive", "False Negative", "Overall Error"),
    col=c("darkblue", "red", "green"), lwd=c(1,1,2), lty=1, cex=1.1)

