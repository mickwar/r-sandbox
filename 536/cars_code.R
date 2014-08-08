# non-linearity
# natural splines
# variable dichotomization
# mean-squared error (for knot selection)
# coverage (repeated random sampling)


library(splines)

count = function(x){ # basically equivalent to table()
    y = unique(x)
    out = double(length(y))
    for (i in 1:length(y))
        out[i] = sum(x == y[i])
    names(out) = y
    return (out)
    }

X = read.csv("./cars.csv")

Price = X$Price
Miles = X$Miles
n = nrow(X)

# fix the problems
X = X[,-c(1,2,3,4,13)] # ID, model, price, age (months), cylinders
    # too many different models, age is redundant, and
    # cylinders only takes on one value
X$cc[81] = 1600 # fix cc mistake

# combine the low numbered colors into a single group
X$Color = as.vector(X$Color)
X$Color = ifelse(X$Color == "Yellow", "Other", X$Color)
X$Color = ifelse(X$Color == "Violet", "Other", X$Color)
X$Color = ifelse(X$Color == "Beige", "Other", X$Color)
# combine weights into fewer groups
plot(X$Weight, Price)
X$Weight = ifelse(X$Weight >= 1200, 3, X$Weight)
X$Weight = ifelse(X$Weight >= 1100, 2, X$Weight)
X$Weight = ifelse(X$Weight >= 1000, 1, X$Weight)
# combine cc into fewer groups
plot(X$cc, Price)
X$cc = ifelse(X$cc >= 1800, 4, X$cc)
X$cc = ifelse(X$cc >= 1500, 3, X$cc)
X$cc = ifelse(X$cc >= 1350, 2, X$cc)
X$cc = ifelse(X$cc >= 1300, 1, X$cc)
# combine HP into fewer groups
plot(X$HP, Price)
X$HP = ifelse(X$HP >= 100, 4, X$HP)
X$HP = ifelse(X$HP >= 90, 3, X$HP)
X$HP = ifelse(X$HP >= 80, 2, X$HP)
X$HP = ifelse(X$HP >= 60, 1, X$HP)
# combine doors into fewer groups
plot(X$Doors, Price)
X$Doors = ifelse(X$Doors == 2, 3, X$Doors)

# set categorical variables as factors
X[,1] = as.factor(X[,1])
for (i in 3:20)
    X[,i] = as.factor(X[,i])

# look at the data, most are categorical, but you can get
# an idea of what's happening. The last 11 variables just
# 0 or 1, so not necessary to really look at.
par(mfrow=c(3,3))
for (i in 2:10)
    plot(X[,i], Price, xlab=names(X)[i])

# look at the non-linearity between miles and price
pdf("./figs/nonlinear.pdf")
plot(X$Miles, Price, xlab="Miles", ylab="Price ($)")
dev.off()

# determine the number of knots to use
MSE = double(6)
for (knots in 0:(length(MSE)-1)){
    spl = ns(X$Miles, df=1+knots)
    mod = lm(Price ~ spl)
    pred = predict(mod)
    MSE[1+knots] = mean((Price - pred)^2)
    }
knots = 2 # "smallest" MSE

pdf("./figs/nonlinear.pdf")
plot(Miles, Price, ylab="Price ($)")
ord = order(Miles)
pred = predict(lm(Price ~ ns(Miles, df=1+knots)))
points(Miles[ord], pred[ord], col='red', type='l', lwd=2)
dev.off()

# add natural spline functions to data matrix
X = cbind(X, ns(X$Miles, df=1+knots))
# remove miles
X = X[,-2]
# rename columns
names(X)[20:22] = c("nsMiles1", "nsMiles2", "nsMiles3")

# plot number of knots vs MSE
pdf("./figs/knots.pdf")
plot(0:(length(MSE)-1), MSE, type='b', xlab="Number of Knots",
    ylab="MSE", main="Choosing Number of Knots for Natural Spline",
    lwd=2, col="darkturquoise")
rect(par("usr")[1], par("usr")[3], par("usr")[2],
    par("usr")[4], col = "gray90")
points(0:(length(MSE)-1), MSE, type='b', lwd=2, col="darkcyan")
dev.off()

# variable selection (Price ~ ns(X$Miles, df=1+knots) ensures the spline
# is used in the model)
lowerMod = lm(Price ~ 1 + nsMiles1 + nsMiles2 + nsMiles3, data=X)
upperMod = lm(Price ~ ., data=X)
mod = step(lowerMod, scope=list(lower=lowerMod, upper=upperMod), k=log(n),
    direction="both", data=X, silent=TRUE)


# compute the coverage with repeated random test sampling
coverage = double(10000)
for (i in 1:length(coverage)){
    index = sample(n, ceiling(n/10))
    pred.X = X[-index,]
    pred.Y = Price[-index]
    test.X = X[index,]
    test.Y = Price[index]
    mod = lm(pred.Y ~ nsMiles1 + nsMiles2 + nsMiles3 + Mfg_Year +
        Automatic_airco + cc + HP + Airco + Mfr_Guarantee + Weight +
        Doors + Automatic + Powered_Windows, data=pred.X)
    pred = predict(mod, newdata=test.X, interval="prediction")
    coverage[i] = mean(apply(cbind(test.Y, pred[,2:3]), 1,
        function(x) x[1] >= x[2] && x[1] <= x[3]))
    }
pdf("./figs/coverage.pdf")
hist(coverage, col='gray', xlab="Coverage", ylab="Density",
    main="Distribution of Coverage", freq=F)
dev.off()

quantile(coverage, c(0.025, 0.975))

pdf("./figs/residuals.pdf")
plot(rstandard(mod), pch=20, ylab="Standardized Residual")
abline(h=c(-3,3), col='gray', lty=2)
dev.off()
#at = identify(rstandard(mod))

bounds = cbind(coef(mod) - qt(0.975, mod$df)*coef(summary(mod))[,2],
    coef(mod) + qt(0.975, mod$df)*coef(summary(mod))[,2])
