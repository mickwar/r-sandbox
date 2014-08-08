library(AUC)

n = nrow(dat)
train.n = ceiling(n * 0.5)
m = 100 # number of cross-validation loops to do
        # (don't make this too high since it will take longer to run,
        # but the higher the more accurate your estimates)
k = 101 # number of cutoffs to use (equally spaced from 0 to 1)
l = 50
pos.rate = matrix(0, m, k)
neg.rate = matrix(0, m, k)
error = matrix(0, m, k)
all.roc = matrix(0, m, l)

auc.total = double(m)
# cross-validation loop
for (i in 1:m){
    training.obs = sample(1:n,train.n)
    test.obs = c(1:n)[-training.obs]
    training.data = dat[training.obs, ]
    test.data = dat[test.obs, ]

    # model with significant parameters
    train.mod = glm(ben ~ [your variables and intercept], family = binomial,
        data = dat, subset = training.obs)

    risk.score = predict(train.mod, newdata=test.data)

    # calcuate predicted probabilites
    p = exp(risk.score)/(1+exp(risk.score))

    # for assessing model adequacy
    ROC = roc(p, as.factor(test.data[,1]))
    auc.total[i] = auc(ROC)
    app = approx(ROC$fpr, ROC$tpr, xout = seq(0,1, length=l))
    all.roc[i,] = app$y

    lower = 0
    upper = 1
    cutoff = seq(lower,upper,length=k)

    # cutoff loop
    for (j in 1:k){
        t.pos = sum(test.data[ ,1] == 1 & p > cutoff[j])
        t.neg = sum(test.data[ ,1] == 0 & p < cutoff[j])
        total.neg = n-length(table(training.obs))-sum(test.data[ ,1])
        total.pos = sum(test.data[ ,1])

        pos.rate[i, j] = 1-t.pos/total.pos
        neg.rate[i, j] = 1-t.neg/total.neg
        error[i, j] = 1-(t.pos+t.neg)/(n-length(table(training.obs)))
        }
    }

# remove problematic loops
remmy = which(is.na(pos.rate[,1]))
if (length(remmy) > 0){
    pos.rate = pos.rate[-remmy,]
    neg.rate = neg.rate[-remmy,]
    error = error[-remmy,]
    }

plot(cutoff, apply(pos.rate, 2, mean), type='l', col="darkblue", lwd=2,
    ylab="Error Rate", xlab="Threshold", main="Predictive Error Rates")
points(cutoff, apply(pos.rate, 2, quantile, 0.975), type='l',
    col="blue", lty=2, lwd=0.5)
points(cutoff, apply(pos.rate, 2, quantile, 0.025), type='l',
    col="blue", lty=2, lwd=0.5)

points(cutoff, apply(neg.rate, 2, mean), type='l', col="red2", lwd=2)
points(cutoff, apply(neg.rate, 2, quantile, 0.975), type='l',
    col="red", lty=2, lwd=0.5)
points(cutoff, apply(neg.rate, 2, quantile, 0.025), type='l',
    col="red", lty=2, lwd=0.5)

points(cutoff, apply(error, 2, mean), type='l', col="darkgreen", lwd=2)
points(cutoff, apply(error, 2, quantile, 0.975), type='l',
    col="green", lty=2, lwd=0.5)
points(cutoff, apply(error, 2, quantile, 0.025), type='l',
    col="green", lty=2, lwd=0.5)

legend(0.4, 1, col=c("darkblue", "red2", "darkgreen"), lwd=c(1,1,2),
    legend=c("False Positive","False Negative","Overall Error"),
    cex = 1.1, lty = 1)

plot(seq(0, 1, length=l), apply(all.roc, 2, mean), type='l', lwd=2,
    ylab="True Positive Rate", col='blue', ylim=c(0, 1),
    main="Receiver Operating Characteristic Curve",
    xlab="False Positive Rate")
points(seq(0, 1, length=l), apply(all.roc, 2, quantile, 0.025),
    type='l', lty=2, col='lightblue')
points(seq(0, 1, length=l), apply(all.roc, 2, quantile, 0.975),
    type='l', lty=2, col='lightblue')
abline(0, 1, lty=2, col="gray", lwd=2)

merror = apply(error, 2, mean)
lerror = apply(error, 2, quantile, 0.025)
uerror = apply(error, 2, quantile, 0.975)
at = which.min(merror)
auc = mean(auc.total)

# false positive - predicting a "success" (1) when really a failure (0)
# false positive - predicting a "failure" (0) when really a success (1)

# estimates for your error rates
# 1: lower bound on overall error
# 2: average overall error
# 3: upper bound on overall error
# 4: area under the ROC curve (closer to 1 the better, if around 0.5
#    or below you're in trouble)
# 5: the cutoff to use
c(lerror[at], merror[at], uerror[at], auc, cutoff[at])
