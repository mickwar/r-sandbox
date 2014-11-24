library(ISLR)
library(tree)
library(randomForest)

dat = Hitters # Hitters is from ISLR package

tree(dat[,19] ~ dat[,c(2,7)])
