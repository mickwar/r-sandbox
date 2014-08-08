source("./textplot.R")
textplot(seq(-4.5, 4.5, length=100), rep(0, 100))
textplot(1:8)#, box=FALSE)


### using multiple characters in textplot()
x = density(rnorm(100000))
y = x$y
x = x$x
hlx = seq(min(x), max(x), length=1000)
hly = rep(0.2, length=length(hlx))
# whichever set of points comes first, that pch is on top
xaxt = as.character(signif(seq(min(x), max(x), by=1),2))
textplot(c(x, hlx), c(y, hly), pch=c(rep("*", length(x)), rep("-", length(hlx))),
    xaxt = xaxt, yaxt = seq(0, 0.4, by = 0.1))
textplot(c(hlx, x), c(hly, y), pch=c(rep("-", length(hlx)), rep("*", length(x))),
    xlim=c(-10,10))

textplot(x <- density(runif(190, 0.999999999, 1)))
textplot(x, ylim=range(x$y))
textplot(x)
textplot(seq(-4.4, 4.4, length=100), rep(0, 100))

# plot with vertical lines, notice that length(pch) = length(x) = length(y)
x = rgamma(10000, 2, 0.8)-2
ints = quantile(x, c(0.025, 0.975))
dens = density(x)
zeroy = seq(0, dens$y[which.min(abs(dens$x))], length=100)
zerox = rep(0, length(zeroy))
lby = seq(0, dens$y[which.min(abs(dens$x-ints[1]))], length=100)
lbx = rep(ints[1], length(lby))
rby = seq(0, dens$y[which.min(abs(dens$x-ints[2]))], length=100)
rbx = rep(ints[2], length(rby))
x = c(zerox, lbx, rbx, dens$x)
y = c(zeroy, lby, rby, dens$y)
pch = c(rep("o", length(zerox)), rep("|", length(lbx)+length(rbx)),
    rep("*", length(dens$x)))
textplot(x, y, pch=pch, ylim=c(0, max(dens$y)), xlim=range(dens$x))
