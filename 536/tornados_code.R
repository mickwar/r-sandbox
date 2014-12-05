# multiple classification
# random forests
# cross-validation

library(randomForest)
library(tree)

dat = read.csv("~/files/R/data/tornados_data.csv",header=TRUE)

dat = dat[-c(120,122,113,116,158,168,148,149,
    400,187,188,200,202,209,179,180,936),]

dat$hour = as.numeric(substr(dat$Time,1,
    nchar(as.character(dat$Time))-6))
dat$TimeofDay = numeric(nrow(dat))
for (i in 1:nrow(dat)){
	if(dat$hour[i]>=22 || dat$hour[i]<=6)
		dat$TimeofDay[i]='Night'
	if(dat$hour[i]>=7 && dat$hour[i]<=12)
		dat$TimeofDay[i]='Morning'
	if(dat$hour[i]>=13 && dat$hour[i]<=21)
		dat$TimeofDay[i]='Day'
    }
dat$TimeofDay = as.factor(dat$TimeofDay)
dat$Month = as.factor(dat$Month)
dat = dat[,-c(1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 19)]



# cross-validation for choosing m
set.seed(1) # this seed gaurentees train contains 2 observations of F-4
train = sample(nrow(dat), ceiling(nrow(dat)*0.75))
sum(dat[train,2] == 4) # 2 observations
error = double(ncol(dat)-1)
for (i in 1:(ncol(dat)-1)){ # ncol(dat)-1 is max number of predictors
    mod = randomForest(factor(Fscale) ~ ., data=dat,
        importance=TRUE, ntree=500, mtry=i, subset=train)
    error[i] = mod$err[500,1]
    }
which.min(error)
  m = 3
mod = randomForest(factor(Fscale) ~ ., data=dat,
    importance=TRUE, ntree=500, mtry=3, subset=train)
yhat = predict(mod, newdata = dat[-train,])
mean(yhat == dat[-train,2]) # test set classification rate

importance(mod)

pdf("./figs/error.pdf")
plot(mod, main="Error Rates", ylim=c(-0.15,1))
abline(h=1-mean(yhat == dat[-train,2]), col="gray20", lty=2)
legend(0, 0.05, legend=c('Out-of-Bag','F0','F1'),
    lty=1, lwd=2, col=1:3, box.lty=0)
legend(150, 0.05, legend=c('F2','F3','F4'), box.lty=0,
    lty=1, lwd=2, col=4:6)
legend(250, 0.05, legend="Test set", box.lty=0, lty=2, col="gray20")
legend(0,0.05,legend=rep(paste(rep("  ",43), collapse=""),3),lty=0)
dev.off()

pdf("./figs/importance.pdf")
par(mar=c(5.1,5.1,4.1,2.1))
barplot(sort(importance(mod)[,7], decreasing=TRUE), horiz=TRUE,
    cex.names=1, col='dodgerblue', main='Variable Importance',
    xlab='Mean Decrease Gini Index', las=1)
dev.off()

pdf("./figs/onetree.pdf", width=8, height=6)
tree = tree(factor(Fscale)~.,data=dat)
par(oma=c(0,0,0,0), mar=c(0,0,2.5,0))
plot(tree, lwd=2)
text(tree, pretty=0)
title("One Tree Hill")
dev.off()

pdf("./figs/pairs.pdf")
pairs(dat[,5:9])
dev.off()

summary(tree)
