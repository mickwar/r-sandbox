library(ISLR)
library(tree)
library(randomForest)

### Beginning of chapter 8 in the 536 book
dat = Hitters # Hitters is from ISLR package
dat = dat[!is.na(dat$Salary),]
names(dat)

pdf("./figs/ex_tree.pdf", width = 18, height = 9)
par(mfrow=c(1,2), mar=c(5.1,4.6,4.1,1.6))
# figure 8.1
tree.hitters = tree(log(Salary) ~ Hits + Years, data = dat)
tree.hitters = prune.tree(tree.hitters, best = 3)
plot(tree.hitters, lwd=2)
text(tree.hitters, cex=2.5)

# figure 8.2
plot(dat$Years, dat$Hits, pch=20, xlab = "Years", ylab = "Hits",
    cex.lab = 2.0, col = "darkorange")
text(2, 117.5, expression(R[1]), cex = 3)
text(11.8, 52.5, expression(R[2]), cex = 3)
text(11.8, 170.5, expression(R[3]), cex = 3)
abline(v=4.5, lwd=3, col="darkgreen")
lines(c(4.5, 50), c(117.5, 117.5), lwd=3, col="darkgreen")
dev.off()

### comparison of linear models and random forests
set.seed(1)
x = runif(100, min = 0, max = 20)
y = double(length(x))
y[x > 0 & x <= 5] = 7 + rnorm(sum(x > 0 & x <= 5))
y[x > 5 & x <= 10] = 3 + rnorm(sum(x > 5 & x <= 10))
y[x > 10 & x <= 15] = 12 + rnorm(sum(x > 10 & x <= 15))
y[x > 15 & x <= 20] = 11 + rnorm(sum(x > 15 & x <= 20))

(gen.tree = tree(y ~ x))
#plot(gen.tree, lwd=2)
#text(gen.tree, cex = 3)

gen.lm = lm(y ~ x)

pdf("figs/compare_tree_lm.pdf", width = 18, height = 9)
par(mfrow=c(1,2))
plot(x, y, pch=20, cex = 1.5)
abline(v = c(5.03505, 10.1633, 14.894), col = 'darkgreen', lwd=3, lty=2)
lines(c(-5, 5.03505), c(7.297, 7.297), lwd=2, col='blue')
lines(c(5.03505, 10.1633), c(2.979, 2.979), lwd=2, col='blue')
lines(c(10.1633, 14.894), c(12.170, 12.170), lwd=2, col='blue')
lines(c(14.894, 25), c(10.610, 10.610), lwd=2, col='blue')

plot(x, y, pch=20, cex = 1.5)
abline(coef(gen.lm), col='red', lwd=3)
dev.off()

### multivariate tree stuff
library(mvpart) # multivariate regression trees
                # multivariate partitioning

data(spider) # what does this do?
mvpart(data.matrix(spider[,1:12]) ~ twigs + water,
    data = spider, method = "mrt", dissim = "euc")

mvpart(data.matrix(spider[,1:12]) ~ water + sand + moss +
    reft + twigs + herbs, data = spider)#, method = "mrt", dissim = "euc")


### read in ben's movie data
source("./read_data.R")


n = nrow(X)
p = ncol(X)

# divide into training and test sets (the random forest
# mechanism will do a bootstrap sample on the training set)
set.seed(1)
train.ind = sample(n, floor(n*0.5))

mod = randomForest(factor(Y) ~ ., data = X, importance = TRUE,
    ntree=500, mtry=floor(sqrt(p)), subset=train.ind)

plot(mod$err[,1], type='l')
yhat = predict(mod, newdata = data.frame(X[-train.ind,]))
1-mean(yhat == Y[-train.ind])

barplot(sort(importance(mod)[,4]), horiz = TRUE, las = 1,
    xlab = "Mean Decrease Gini Index", cex.names = 1,
    col = "dodgerblue", main = "Variable Importance")

# error rate for a single tree
tre = tree(factor(Y) ~ ., data = data.frame(X), subset = train.ind)
plot(tre)
text(tre)

yhat.tre = predict(tre, newdata = data.frame(X[-train.ind,]))
1-mean(yhat.tre == Y[-train.ind])
