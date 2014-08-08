x = seq(-3, 3, length=25)
y = 2 * x + 1 + rnorm(length(x))

plot(x, y, pch=20)

mod = lm(y ~ x)
abline(coef(mod)[1], coef(mod)[2], lwd=2)

est = summary(mod)[[4]]

# bootstrap samples
B = 100000
x0 = double(B) # intercepts
x1 = double(B) # slopes
se0 = double(B) # s.e. of intercept
se1 = double(B)
for (b in 1:B){
    bootstrap = sample(length(x), replace=TRUE)
    bootmod = lm(y[bootstrap] ~ x[bootstrap])
    ss = summary(bootmod)[[4]]
    x0[b] = ss[1,1]
    x1[b] = ss[2,1]
    se0[b] = ss[1,2]
    se1[b] = ss[2,2]
    }

rbind(t(c(mean(x0), mean(se0))), t(c(mean(x1), mean(se1))))
est

hist(x0, breaks=100)
hist(x1, breaks=100)
hist(se0, breaks=100)
hist(se1, breaks=100)
