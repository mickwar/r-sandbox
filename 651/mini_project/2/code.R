library(truncnorm)

dat = read.table("./faculty.dat")
colnames(dat)=NULL
dat = as.matrix(dat)

f = function(x)
    0.05*dnorm(x, 4.5, 0.25) + 0.40*dnorm(x, 5.4, 0.30) + 0.55*dnorm(x, 6.15, 0.30)

x = seq(0, 7, length=1000)

plot(density(dat))
curve(dnorm(x, mean(dat), sd(dat)), col='blue', add=TRUE)
curve(dtruncnorm(x, 1, 7, mean(dat), sd(dat)), col='green', add=TRUE)
points(x, f(x), type='l', col='red')
