dat = read.table("~/files/data/health.txt", header=TRUE)

y = ifelse(dat$sick == "none", 0, 1)
x = dat[,2:ncol(dat)]

int.mod = glm(y ~ 1, data = x, family = binomial)
full.mod = glm(y ~ ., data = x, family = binomial)
summary(mod)

apply(x, 2, sum)

mod = step(int.mod, scope = list("lower"=int.mod, "upper"=full.mod), k = log(length(y)))

