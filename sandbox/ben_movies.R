library(AUC)
library(splines)
library(car)

# data cleaning
dat = read.table("../data/movie.txt", header=T)
pred = read.table("../data/movie_pred.txt", header=T)
pred_movies = pred[,1]
pred = pred[,-1]
movies = dat[,2]
dat = dat[,-2]

dat = cbind(dat, "kids_I"=is.na(dat$kids_S)*1)
for (i in 2:4)
    dat[,i] = ifelse(is.na(dat[,i]), 0, dat[,i])
dat = cbind(dat, "rot_I"=is.na(dat$rot_top)*1)
for (i in 6:8)
    dat[,i] = ifelse(is.na(dat[,i]), 0, dat[,i])

# add indicators to the prediction set as well
pred = cbind(pred, "kids_I"=is.na(pred$kids_S)*1)
for (i in 1:3)
    pred[,i] = ifelse(is.na(pred[,i]), 0, pred[,i])
pred = cbind(pred, "rot_I"=is.na(pred$rot_top)*1)
for (i in 5:7)
    pred[,i] = ifelse(is.na(pred[,i]), 0, pred[,i])

##### deal with collineraity via orthogonalization
# orthogonalize the same variables in the prediction set with the
# data set. (is this how to do prediction when orthogonalizing
# variables?)
n = nrow(dat)
m = nrow(pred)
temp1 = c(dat$rot_all, pred$rot_all)
temp2 = c(dat$rot_top, pred$rot_top)
pred$rot_all = as.matrix((temp1 - temp2 %*% solve(t(temp2) %*% temp2) %*% 
    t(temp2) %*% temp1)[(n+1):(n+m)])
names(pred)[which(names(pred) == "rot_all")] = "rot_allO"

temp1 = c(dat$rot_top, pred$rot_top)
temp2 = c(dat$rot_aud, pred$rot_aud)
pred$rot_top = as.matrix((temp1 - temp2 %*% solve(t(temp2) %*% temp2) %*% 
    t(temp2) %*% temp1)[(n+1):(n+m)])
names(pred)[which(names(pred) == "rot_top")] = "rot_topO"

temp1 = c(dat$rot_aud, pred$rot_aud)
temp2 = c(dat$imdb, pred$imdb)
pred$rot_aud = as.matrix((temp1 - temp2 %*% solve(t(temp2) %*% temp2) %*% 
    t(temp2) %*% temp1)[(n+1):(n+m)])
names(pred)[which(names(pred) == "rot_aud")] = "rot_audO"


ortho = dat$rot_all - dat$rot_top %*% solve(t(dat$rot_top) %*%
    dat$rot_top) %*% t(dat$rot_top) %*% dat$rot_all
dat$rot_all = ortho
names(dat)[which(names(dat) == "rot_all")] = "rot_allO"

ortho = dat$rot_top - dat$rot_aud %*% solve(t(dat$rot_aud) %*%
    dat$rot_aud) %*% t(dat$rot_aud) %*% dat$rot_top
dat$rot_top = ortho
names(dat)[which(names(dat) == "rot_top")] = "rot_topO"

ortho = dat$rot_aud - dat$imdb %*% solve(t(dat$imdb) %*%
    dat$imdb) %*% t(dat$imdb) %*% dat$rot_aud
dat$rot_aud = ortho
names(dat)[which(names(dat) == "rot_aud")] = "rot_audO"

pairs(dat[,c(2,3,4,6,7,8,10)])

# functions
mw.smooth = function(x, y, d=1){
    kern = function(x, y, d)
        exp(-1/(2*d^2)*(x-y)^2)
    m = 8*length(x)
    outx = seq(min(x), max(x), length=m)
    outy = double(m)
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i], d))/
            sum(kern(x, outx[i], d))
    return (list("x"=outx, "y"=outy))
    }

# assumption checking
# check for non-monotonicity
# possible issues in kids_P, rot_all
n = nrow(dat)

CC = dat$kids_S
SS = sort(unique(CC))
means = double(length(SS))
for (i in 1:length(means))
    means[i] = mean(dat[which(CC == SS[i]), 1])
smooth = mw.smooth(SS, means, 0.5)
plot(CC+rnorm(n,0,0.15),dat$ben+rnorm(n,0,0.015),pch=20,cex=0.5)
lines(smooth$x, smooth$y, col='red')

CC = dat$rot_top
SS = sort(unique(CC))
means = double(length(SS))
for (i in 1:length(means))
    means[i] = mean(dat[which(CC == SS[i]), 1])
smooth = mw.smooth(SS, means, 5)
plot(CC+rnorm(n,0,0.15),dat$ben+rnorm(n,0,0.015),pch=20,cex=0.5)
lines(smooth$x, smooth$y, col='red')

CC = dat$rot_all
SS = sort(unique(CC))
means = double(length(SS))
for (i in 1:length(means))
    means[i] = mean(dat[which(CC == SS[i]), 1])
smooth = mw.smooth(SS, means, 5)
plot(CC+rnorm(n,0,0.15),dat$ben+rnorm(n,0,0.015),pch=20,cex=0.5)
lines(smooth$x, smooth$y, col='red')

CC = dat$rot_aud
SS = sort(unique(CC))
means = double(length(SS))
for (i in 1:length(means))
    means[i] = mean(dat[which(CC == SS[i]), 1])
smooth = mw.smooth(SS, means, 5)
plot(CC+rnorm(n,0,0.15),dat$ben+rnorm(n,0,0.015),pch=20,cex=0.5)
lines(smooth$x, smooth$y, col='red')



full.mod = glm(ben ~ kids_P + rot_allO + kids_S +
    kids_V + rot_topO + rot_audO + MPAA + kids_I +
    rot_I + live_action, data=dat, family=binomial)
null.mod = glm(ben ~ 1, data=dat, family=binomial)

mod = step(null.mod, scope=list(lower=null.mod, upper=full.mod),
    k=log(nrow(dat)), data=dat, direction="both", family=binomial)
summary(mod)
vif(mod)
vif(full.mod)

risk.score = predict(mod, newdata=dat)
p = exp(risk.score)/(1+exp(risk.score))

data.frame(movies, p)

risk.score = predict(mod, newdata=pred)
p = exp(risk.score)/(1+exp(risk.score))

n = nrow(dat)
train.n = ceiling(n * 0.5)
m = 100
k = 100
l = 50
pos.rate = matrix(0, m, k)
neg.rate = matrix(0, m, k)
error = matrix(0, m, k)
all.roc = matrix(0, m, l)

auc.total = double(m)
for (i in 1:m){
    # strict subset
    training.obs = sample(1:n,train.n)
    # bootstrap
#   training.obs = sample(n, train.n, replace=TRUE)
    test.obs = c(1:n)[-training.obs]
    training.data = dat[training.obs, ]
    test.data = dat[test.obs, ]

    # model with significant parameters
#   train.mod = glm(ben ~ 1 + kids_S + rot_topO, family = binomial,
#       data = dat, subset = training.obs)
    # model that just fits everything
    train.mod = glm(ben ~ kids_P + rot_allO + kids_S +
        kids_V + rot_topO + rot_audO + MPAA + kids_I +
        rot_I + live_action, data=dat, family=binomial)

    risk.score = predict(train.mod, newdata=test.data)

    p = exp(risk.score)/(1+exp(risk.score))
    ROC = roc(p, as.factor(test.data[,1]))
    auc.total[i] = auc(ROC)
#   plot(ROC$fpr, ROC$tpr, ylab="True Positive Rate", type='s', lwd=2,
#       main="Receiver Operating Characteristic Curve",
#       xlab="False Positive Rate")
#   abline(0, 1, lty=2, col="gray")
    app = approx(ROC$fpr, ROC$tpr, xout = seq(0,1, length=l))
    all.roc[i,] = app$y
#   points(app$x, app$y, ylab="True Positive Rate", type='l', lwd=2,
#       main="Receiver Operating Characteristic Curve", col='blue',
#       xlab="False Positive Rate")
#   abline(0, 1, lty=2, col="gray")

#   lower = min(p)
#   upper = max(p)
    lower = 0
    upper = 1
    cutoff = seq(lower,upper,length=k)

    for (j in 1:k){
        t.pos = sum(test.data[ ,1] == 1 & p >= cutoff[j])
        t.neg = sum(test.data[ ,1] == 0 & p < cutoff[j])
        total.neg = n-length(table(training.obs))-sum(test.data[ ,1])
        total.pos = sum(test.data[ ,1])

        pos.rate[i, j] = 1-t.pos/total.pos
        neg.rate[i, j] = 1-t.neg/total.neg
        error[i, j] = 1-(t.pos+t.neg)/(n-length(table(training.obs)))
        }

#   readline()
    }
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
c(lerror[at], merror[at], uerror[at], auc, cutoff[at])

merror2 = apply(error, 2, mean)
lerror2 = apply(error, 2, quantile, 0.025)
uerror2 = apply(error, 2, quantile, 0.975)
at2 = which.min(merror2)
auc2 = mean(auc.total)

rbind(c(lerror[at], merror[at], uerror[at], auc, cutoff[at]),
    c(lerror2[at2], merror2[at2], uerror2[at2], auc2, cutoff[at2]))

# false positive - predicting a "like" when not really liked
# false negative - predicting a "dislike" when really liked

risk.score = predict(mod)
p1 = exp(risk.score)/(1+exp(risk.score))

risk.score = predict(full.mod)
p2 = exp(risk.score)/(1+exp(risk.score))
p2 = ifelse(p2 < 1e-5, 0, p2)
data.frame(p1, p2, dat[,1], movies)





