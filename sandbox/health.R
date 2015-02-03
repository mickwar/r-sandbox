dat = read.table("~/files/data/health.txt", header=TRUE)

y = ifelse(dat$sick == "none", 0, 1)
x = dat[,2:ncol(dat)]

mod = glm(y ~ 0 + time, data = x, family = binomial)
summary(mod)

apply(x, 2, sum)
