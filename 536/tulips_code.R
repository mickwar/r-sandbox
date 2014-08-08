# logistic regression
# "non-linearity" (i.e. non-monotonicity, for logistic)
# natural splines
# roc / auc
# bootstrapping

library(splines)
library(car)
library(AUC)
library(boot)
set.seed(1)

dat = read.csv("./tulips.csv")

# remove year/day variables
dat = dat[,-c(2,3)]
dat$Germinated = as.numeric(dat$Germinated)-1

# doesn't say much
#pairs(dat)

mw.smooth = function(x, y, method="gaussian"){
    if (method=="gaussian")
        kern = function(x, y, d=16/16)
            exp(-1/(2*d^2)*(x-y)^2)
    m = 8*length(x)
    outx = seq(min(x), max(x), length=m)
    outy = double(m)
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i]))/
            sum(kern(x, outx[i]))
    return (list("x"=outx, "y"=outy))
    }

calc.ideal = function(predictions, drop = 0.05,
    seq.x = seq(0, 12, length=100)){
    n = length(seq.x)
    p = predictions
    mode.ind = which.max(p)
    lower.ind = mode.ind
    upper.ind = mode.ind
    while (lower.ind > 1 && p[lower.ind] > p[mode.ind] - drop)
        lower.ind = lower.ind - 1
    while (upper.ind < n && p[upper.ind] > p[mode.ind] - drop)
        upper.ind = upper.ind + 1
    return (c(seq.x[lower.ind], seq.x[mode.ind], seq.x[upper.ind]))
    }


calc.full = function(i, j){
    keep.auc = double(6)
    ind = c(index[[i]], index[[j]])
    for (knots in 0:(length(keep.auc)-1)){
        temp.auc = double(100)
        for (k in 1:length(temp.auc)){
            pop1.train = sample(1:210, ceiling(210*0.85))
            pop1.test = (1:210)[-pop1.train]
            pop2.train = sample(1:210, ceiling(210*0.85))
            pop2.test = (1:210)[-pop.train]
            pop2.train = pop2.train + 210
            pop2.test = pop2.test + 210
            ind.train = ind[c(pop1.train, pop2.train)]
            ind.test = ind[-c(pop1.train, pop2.train)]
            mod = glm(Germinated ~ factor(Population)+ns(ChillingTime,
                df=1+knots)+factor(Population)*ns(ChillingTime,df=1+knots),
                data = dat[ind.train,], family=binomial)
            p = 1/(1+exp(-predict(mod, newdata=dat[ind.test,])))
            temp.auc[k] = auc(roc(p, as.factor(dat[ind.test,3])))
            }
        keep.auc[knots+1] = mean(temp.auc)
        }
    knot = which.max(keep.auc)-1
    mod = glm(Germinated ~ factor(Population)+ns(ChillingTime,
        df=1+knot)+factor(Population)*ns(ChillingTime,df=1+knot),
        data = dat[ind,], family=binomial)
    return (list("mod"=mod, "knot"=knot))
    }

calc.mode = function(x, precision=512, method="density"){
    d = density(x, n = precision)
    dx = d$x
    dy = d$y
    dx[which.max(dy)]
    }

########
# non-linearity in chilling time (across all populations)
# non-linearity as in non-monotonic
chilling = sort(unique(dat$ChillingTime))
means = double(length(table(dat$ChillingTime)))
for (i in 1:length(means))
    means[i] = mean(dat[which(dat$ChillingTime == chilling[i]),3])
# plot(jitter(dat$ChillingTime, 1.5), jitter(dat$Germinated, 0.20),
#     pch=20, cex=0.1)
allsmooth = mw.smooth(chilling, means)
# lines(allsmooth$x, allsmooth$y, col='red')

# chilling time across all populaions
chilling = sort(unique(dat$ChillingTime))
means = matrix(0, 12, length(table(dat$ChillingTime)))
for (i in 1:nrow(means)){
    for (j in 1:ncol(means)){
        means[i, j] = mean(dat[which(dat$ChillingTime ==
            chilling[j]  & dat$Population == i),3])
        }
    }
pdf("./figs/monotone.pdf", width=12, height=6)
seq.x = seq(0, 12, length=100)
par(mfrow=c(1,2))
plot(seq.x, 1/(1+exp(-(seq.x-6))), type='l', xlab="", ylab="",
    main="Logistic Function")
plot(sample(dat$ChillingTime, 840)+rnorm(840,0,0.15), sample(
    dat$Germinated, 840)+rnorm(840,0,0.015), pch=20, cex=0.1,
    xlab="Chilling Time", ylab="Germinated", main=paste("Smoothed",
    "Germination Rates"))
smoothed = apply(means, 1, function(x) mw.smooth(chilling, x))
for (i in 1:12)
    lines(smoothed[[i]]$x, smoothed[[i]]$y,
        col=rainbow(nrow(means))[i], lwd=1)
lines(allsmooth$x, allsmooth$y, col='black', lwd=3)
dev.off()
##########
   
##########
# overall model (across all populations)
keep.auc = double(6)
for (knots in 0:(length(keep.auc)-1)){
    temp.auc = double(100)
    for (k in 1:length(temp.auc)){
        train = sample(nrow(dat), ceiling(nrow(dat)*0.75))
        test = (1:nrow(dat))[-train]
        all.mod = glm(Germinated ~ Population+ns(ChillingTime,
            df=1+knots)+Population*ns(ChillingTime,df=1+knots),
            data = dat[train,], family=binomial)
        p = 1/(1+exp(-predict(all.mod, newdata=dat[test,])))
        temp.auc[k] = auc(roc(p, as.factor(dat[test,3])))
        }
    keep.auc[knots+1] = mean(temp.auc)
    }
knot = which.max(keep.auc)-1
all.mod = glm(Germinated ~ Population+ns(ChillingTime,
    df=1+knot)+Population*ns(ChillingTime,df=1+knot),
    data = dat, family=binomial)
red.mod = glm(Germinated ~ ns(ChillingTime, df=1+knot),
    data=dat, family=binomial)
all.pval = anova(red.mod, all.mod, test="Chisq")[2,5]

##########

index = NULL
for (i in 1:12)
    index[[i]] = which(dat$Population == i)
# pop.train = sample(210, ceiling(210*0.75))
# pop.test = (1:210)[-train]
# dat[index[[1]][pop.train],]

# full model, each get own curve (more parameters)
# reduced model, get same curve (less parameters)

# same/different chilling time effect in populations
# ignore population 12, why? obvious, it's a bad population
# bonferroni: choose(11, 2) = 55 tests
pvals = matrix(Inf, 11, 11)
for (i in 1:10){
    for (j in (i+1):11){
        ind = c(index[[i]], index[[j]])
        full.mod = calc.full(i, j)
        red.mod = glm(Germinated ~ ns(ChillingTime, df=1+
            full.mod$knot), data=dat[ind,], family=binomial)
        pvals[i,j] = anova(red.mod, full.mod$mod, test="Chisq")[2,5]
        pvals[j,i] = pvals[i,j]
        }
    }
tests = (pvals < 0.05/choose(11,2))*1+3*diag(11)
pdf("./figs/chill_test.pdf")
plot(1:11, 1:11, type='l', xlim=c(1,11), ylim=c(1,11),
    xlab="Population", ylab="Population", axes=FALSE,
    main="Likelihood Ratio", lwd=2)
for (i in 1:11)
    for (j in 1:11)
        points(i, j, pch=15, col=4-tests[i,j], cex=3)
axis(1, at=1:11)
axis(2, at=1:11)
mtext("Green: difference -- Blue: similarity", line=0.5)
dev.off()

A = c(2,3,4)
B = c(6,7,10,11)
new.index = NULL
for (i in 1:4)
    new.index[[i]] = index[[c(1,5,8,9)[i]]]
new.index[[5]] = c(index[[A[1]]], index[[A[2]]], index[[A[3]]])
new.index[[6]] = c(index[[B[1]]], index[[B[2]]], index[[B[3]]],
    index[[B[4]]])
# models for each population
mod = NULL
# knot.keep = c(3,4,5,3,5,4,3,0,2,5,5,0)
# pop.train = sample(210, ceiling(210*0.85))
# pop.test = (1:210)[-train]
# for (i in 1:11){
#     mod[[i]] = glm(Germinated ~ ns(ChillingTime,df=1+knot.keep[i]),
#         data=dat[index[[i]][pop.train],], family=binomial)
#     }
# mod[[12]] = glm(Germinated ~ ChillingTime, data=dat[index[[12]][
#     pop.train],], family=binomial)
knot.keep = double(6)
auc.keep = double(6)
for (i in 1:6){ # mod[[12]] done separately
    keep.auc = double(6)
    for (knots in 0:5){
        temp.auc = double(100)
        set.seed(1)
        for (j in 1:100){
            pop.train = sample(length(new.index[[i]]), ceiling(
                length(new.index[[i]])*0.85))
            pop.test = (1:length(new.index[[i]]))[-pop.train]
            mod[[i]] = glm(Germinated ~ ns(ChillingTime,df=1+knots),
                data=dat[new.index[[i]][pop.train],],
                family=binomial)
            p = 1/(1+exp(-predict(mod[[i]], newdata=dat[new.index[[
                i]][pop.test],])))
            temp.auc[j] = auc(roc(p, as.factor(dat[new.index[[i]][
                pop.test],3])))
            }
        keep.auc[1+knots] = mean(temp.auc)
        }
    knot.keep[i] = which.max(keep.auc)-1
    auc.keep[i] = max(keep.auc)
    mod[[i]] = glm(Germinated ~ ns(ChillingTime, df=which.max(
        keep.auc)), data=dat[new.index[[i]],],
        family=binomial)
    }

pdf("./figs/model_fit.pdf")
seq.x = seq(0, 12, length=100)
chill.x = dat[1:length(seq.x),]
plot(sample(dat$ChillingTime, 840)+rnorm(840,0,0.15), sample(
    dat$Germinated, 840)+rnorm(840,0,0.015), pch=20, cex=0.1,
    xlab="Chilling Time", ylab="Probability of Germination",
    main="Model Fit Across Populations")
for (i in 1:6){
    chill.x$ChillingTime = seq.x
    p = 1/(1+exp(-predict(mod[[i]], newdata=chill.x)))
    points(seq.x, p, type='l', lwd=1, col=rainbow(nrow(means))[i])
    points(seq.x[seq(1,100,by=6)], p[seq(1,100,by=6)],
        col=rainbow(nrow(means))[i], pch=c(1,5,8,9,"A","B")[i])
    }
dev.off()

# get uncertainties on ideal chilling time for each population
seq.x = seq(0, 12, length=1000)
bounds = matrix(0, 6, 4)
chill.x = dat[1:length(seq.x),]
for (i in 1:6){
#   chill.x[,1] = i
    chill.x$ChillingTime = seq.x
    p = 1/(1+exp(-predict(mod[[i]], newdata=chill.x)))
    bounds[i,] = c(calc.ideal(p, drop=0.05, seq.x), max(p))
    }
bounds = round(bounds, 3)
# write.table(bounds, "./bounds.txt", row.names=FALSE,
#     col.names=FALSE)


### bootstrapping
chill.y = dat[1:2,]
chill.y$ChillingTime = c(10, 8)
B = 1000
maxes = matrix(0, 6, B)
probs = matrix(0, 6, B)
diffs = matrix(0, 6, B)
for (i in 1:6){
    for (b in 1:B){
        bsamp = sample(new.index[[i]], ceiling(length(new.index[[
            i]])*0.5), replace=TRUE)
        mod[[i]] = glm(Germinated ~ ns(ChillingTime,
            df=1+knot.keep[i]), data=dat, subset=bsamp,
            family=binomial)
        p = 1/(1+exp(-predict(mod[[i]], newdata=chill.x)))
        maxes[i,b] = seq.x[which.max(p)]
        probs[i,b] = max(p)
        diffs[i,b] = diff(1/(1+exp(-predict(mod[[i]], newdata=chill.y))))
        }
    }

bounds = cbind(t(apply(maxes, 1, quantile, c(0.025, 0.975))),
    apply(maxes, 1, calc.mode))

prob.bounds = cbind(t(apply(probs, 1, quantile, c(0.025, 0.975))),
    apply(probs, 1, calc.mode))

diff.bounds = cbind(t(apply(diffs, 1, quantile, c(0.025, 0.975))),
    apply(diffs, 1, calc.mode))




