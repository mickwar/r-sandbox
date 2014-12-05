### load and clean the data
source("./clean.R")
library(randomForest)

n = nrow(X)
p = ncol(X)

# divide into training and test sets (the random forest
# mechanism will do a bootstrap sample on the training set)
set.seed(1)
train.ind = sort(sample(n, floor(n*0.57)))

set.seed(1)
mod = randomForest(factor(Y) ~ ., data = X, importance = TRUE,
    ntree=500, mtry=floor(sqrt(p)), subset=train.ind, proximity = TRUE)
yhat = predict(mod, newdata = data.frame(X[-train.ind,]))
# test set error rate
(test.error = 1-mean(yhat == Y[-train.ind]))

# confusion matrix for the test set
# assumed predicted[i] corresponds with true[i]
confuse = function(true, predicted){
    names = names(table(true))
    k = length(names)
    out = matrix(0, k, k)
    for (i in 1:length(true)){
        atx = which(true[i] == names)
        aty = which(predicted[i] == names)
        out[atx, aty] = out[atx, aty] + 1
        }
    errors = double(k)
    for (j in 1:k)
        errors[j] = 1 - out[j,j] / sum(out[j,])
    out = cbind(out, errors)
    rownames(out) = names
    colnames(out)[1:k] = names
    return (out)
    }
confuse(Y[-train.ind], yhat)

# error rate for a single tree
library(tree)
tre = tree(factor(Y) ~ ., data = data.frame(X), subset = train.ind)
yhat.tre = predict(tre, newdata = data.frame(X[-train.ind,]))
one.error = 1-mean(yhat.tre == Y[-train.ind])

pdf("./figs/forest_error.pdf", height = 7.5, width = 7.5)
par(mar = c(5.1,5.1,4.1,4.1))
plot(mod, type='l', cex.lab = 1.5, cex.main = 2, lwd = 1.5, lty=1,
    main = "Average classification error rate")
abline(h = test.error, lwd=2, lty=2)
abline(h = one.error, lwd=2, lty=2, col="red")
mtext("OOB", 1, at = 548, line = -9.0)
mtext("Test set", 1, at = 555, line = -10.5)
mtext("One tree", 1, at = 555, line = -25.0)
mtext("Class 1", 1, at = 555, line = -2.5)
mtext("Class 3", 1, at = 555, line = -28.0)
mtext("Class 4", 1, at = 555, line = -21.5)
dev.off()

# confusion matrix (on the training set)
mod$conf

# variable importance
pdf("./figs/importance2.pdf", height = 7.5, width = 15)
par(mar=c(5.1,5.1,4.1,2.1), mfrow=c(1, 2))
set.seed(1)
mod = randomForest(factor(Y) ~ ., data = X, importance = TRUE,
    ntree=200000, mtry=floor(sqrt(p)), subset=train.ind, proximity = TRUE)
barplot(sort(importance(mod)[,4]), horiz = TRUE, las = 1,
    xlab = "Mean Decrease in Gini Index", cex.names = 1, cex.main = 2,
    col = "dodgerblue", main = "Variable Importance", cex.lab = 1.5)
# dev.off()
### a single tree
# pdf("./figs/mov_tree.pdf")
tre = tree(factor(Y) ~ ., data = data.frame(X))
plot(tre, lwd = 1)
text(tre, cex = 1.5)
dev.off()

# proximity (something is going on here, the proximity matrix
# is supposed to be positive definite. only include training set?)
set.seed(1)
full.mod = randomForest(factor(Y) ~ ., data = X, importance = TRUE,
    ntree=500, mtry=floor(sqrt(p)), proximity = TRUE)
eig = eigen(full.mod$proximity)$values
det(full.mod$proximity)
