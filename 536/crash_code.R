# logistic regression
# roc curve / auc
# predictive error rates (across thresholds)

library(car) # for vif
library(AUC)

dat = read.csv("./crash.csv")

# data cleaning

# remove unnecessary variables
# all years are 2012
dat = dat[,-2]

# get rid of hour 99 observations (only 23, shouldn't
# be a problem)
dat = dat[-which(dat$Hour == 99),]
# remove 9999 Mod_year, only 4 observations
dat = dat[-which(dat$Mod_year == 9999),]
# remove the 2 unknowns from light status
dat = dat[-which(dat$Light == "Unknown"),]

# combine DUI, 0 is no dwi, 1 is at least 1 dwi
# only 7 with 2 dwis and 1 with 3 dwis
dat$DWI = ifelse(dat$DWI > 0, 1, 0)
dat$DWI = factor(dat$DWI)

# combine car_types (some groups have low number of
# observations
levels(dat$Car.Type) = c(levels(dat$Car.Type), "Car", "Midsize")
for (i in 1:nrow(dat)){
    if (dat$Car.Type[i] == "2-Door Sedan" || dat$Car.Type[i] == "4-Door Sedan" ||
        dat$Car.Type[i] == "Convertible")
        dat$Car.Type[i] = "Car"
    if (dat$Car.Type[i] == "Jeep" || dat$Car.Type[i] == "Small Truck" ||
        dat$Car.Type[i] == "Station Wagon" || dat$Car.Type[i] == "Large Van" ||
        dat$Car.Type[i] == "Minivan")
        dat$Car.Type[i] = "Midsize"
    }
dat$Car.Type = factor(dat$Car.Type)

# combine dawn and dusk (similar visibility conditions)
levels(dat$Light) = c(levels(dat$Light), "Dawn/Dusk")
for (i in 1:nrow(dat))
    if (dat$Light[i] == "Dawn" || dat$Light[i] == "Dusk")
        dat$Light[i] = "Dawn/Dusk"
dat$Light = factor(dat$Light)

# combine some seat belt groups ?

# treat speed limit as continuous

# check for non-linear (do not need to check with categorical)

# check for collinearity (how to check for categorical variables?)
# possible collinearity between:
#   height and weight
#   hour and light

# consider combining month and day of month to make a holiday indicator

# plot(jitter(dat$Hour), jitter(dat$Fatal), col=dat$Fatal+1)

# variable selection, hybrid selection, BIC
full.mod = glm(Fatal ~ ., data=dat, family=binomial)
null.mod = glm(Fatal ~ 1, data=dat, family=binomial)

mod = step(null.mod, scope=list(lower=null.mod, upper=full.mod),
    k=log(nrow(dat)), data=dat, direction="both", family=binomial)

bounds = cbind(coef(mod) - qnorm(0.975)*coef(summary(mod))[,2],
    coef(mod) + qnorm(0.975)*coef(summary(mod))[,2])

pdf("./figs/nonlinear.pdf")
speeds = sort(unique(dat$Speed.Limit))
means = double(length(table(dat$Speed.Limit)))
for (i in 1:length(means))
    means[i] = mean(dat[which(dat$Speed.Limit == speeds[i]),1])
plot(jitter(dat$Speed.Limit, 1.0), jitter(dat$Fatal, 0.1),
    pch=20, cex=0.1, xlab="Speed Limit", ylab="Fatality")
lines(speeds, means, col='red')
dev.off()

vif(mod)

pdf("./figs/speedyear.pdf", width=12, height=6)
par(mfrow=c(1,2))
speed.dat = dat[1:100,]
for (i in 1:ncol(speed.dat))
    speed.dat[,i] = speed.dat[2,i]
a = min(dat$Speed.Limit)
b = max(dat$Speed.Limit)
speed.dat$Speed.Limit = seq(a, b, length=100)
speed.preds = 1/(1+exp(-predict(mod, newdata=speed.dat)))
plot(seq(a, b, length=100), speed.preds, type='l',
    xlab="Speed Limit", ylab="Predicted Probability of Death", main="Speed Limit Risk")

year.dat = dat[1:100,]
for (i in 1:ncol(year.dat))
    year.dat[,i] = year.dat[2,i]
a = min(dat$Mod_year)
b = max(dat$Mod_year)
year.dat$Mod_year = seq(a, b, length=100)
year.preds = 1/(1+exp(-predict(mod, newdata=year.dat)))
plot(seq(a, b, length=100), year.preds, type='l',
    xlab="Model Year", ylab="", main="Vehicle Model Year Risk")
dev.off()

n = nrow(dat)
train.n = ceiling(n * 0.75)
m = 100

auc.total = double(m)
for (i in 1:m){
    training.obs = sample(1:n,train.n)
    test.obs = c(1:n)[-training.obs]
    training.data = dat[training.obs, ]
    test.data = dat[test.obs, ]

    train.mod = glm(formula = Fatal ~ Drugs + Belt + Speed.Limit + Speed.Related + 
        Light + Drink + Distracted + Mod_year, family = binomial, data = dat)

    risk.score = predict(train.mod, newdata=test.data)

    p = exp(risk.score)/(1+exp(risk.score))
    ROC = roc(p, as.factor(test.data[,1]))
    auc.total[i] = auc(ROC)
    }
mean(auc.total)

pdf("./figs/roc.pdf")
plot(ROC$fpr, ROC$tpr, main="Receiver Operating Characteristic Curve",
    ylab="True Positive Rate", xlab="False Positive Rate", type='s', lwd=2)
abline(0, 1, lty=2, col="gray")
legend(0.45, 0.2, legend="Average AUC = 0.92", cex=1.5, box.lty=0)
dev.off()


m = 500
lower = min(p)
upper = max(p)
cutoff = seq(lower,upper,length=m)
pos.rate = double(m)
neg.rate = double(m)
error = double(m)

for (i in 1:m){
    t.pos = sum(test.data[ ,1] == 1 & p > cutoff[i])
    t.neg = sum(test.data[ ,1] == 0 & p < cutoff[i])
    total.neg = (n-train.n)-sum(test.data[ ,1])
    total.pos = sum(test.data[ ,1])

    pos.rate[i] = 1-t.pos/total.pos
    neg.rate[i] = 1-t.neg/total.neg
    error[i] = 1-(t.pos+t.neg)/(n-train.n)
    }

pdf("./figs/error.pdf")
plot(cutoff, pos.rate, type='s', col="darkblue", ylab="Error Rate",
    xlab="Threshold", main="Predictive Error Rates")
points(cutoff, neg.rate, type='s', col="red")
points(cutoff, error, type='s', col="green", lwd=2)
legend(0.6, 1, legend=c("False Positive", "False Negative", "Overall Error"),
    col=c("darkblue", "red", "green"), lwd=c(1,1,2), lty=1, cex=1.1)
dev.off()
