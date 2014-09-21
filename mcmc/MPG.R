# the inverse gamma density
igpdf = function(x, shape, scale)
    (1/(gamma(shape)*scale^shape))*x^(-shape-1)*exp(-1/(x*scale))
    
# returns the posterior parameters for mu, depending
# on the current value of sigma2, mumu and sigma2mu are known
updatemu = function(data, sigma2, mumu, sigma2mu){
    xbar = mean(data)
    n = length(data)
    precision = 1/sigma2
    precisionmu = 1/sigma2mu
    postmu = ((xbar*n*precision+mumu*precisionmu)/(n*precision+precisionmu))
    postvar = (1/(n*precision+precisionmu))
    return (c(postmu, postvar))
    }

# returns the posterior parameters for sigma2, depending
# on the current value of mu, priorSh and priorSc are known
updatesigma2 = function(data, mu, priorSh, priorSc){
    n = length(data)
    SSx = sum((data-mu)^2)
    postSh = n/2+priorSh
    postSc = 2*priorSc/(priorSc*SSx+2)
    return (c(postSh, postSc))
    }

# function that contains the looping, at each iteration the
# parameter is update based on the most recent value of the
# other parameters (only two parameters in this case)
gibbs = function(data, mumu, sigma2mu, prsh, prsc, loops=10000){
    # initialize the matrix
    out = matrix(0, nrow=loops, ncol=2)
    # choose staring location for mu to be the mean
    out[1,1] = mean(data)
    # sigma2 is not initialize since we'll update first thing
    # in the loop
    for (i in 1:loops){
        post = updatesigma2(data,out[i-1,1],prsh,prsc)
        # 1/rgamma yields randoms draws from inverse gamma
        # get random draw for sigma2
        out[i,2] = 1/rgamma(1,shape=post[1],scale=post[2])
        post = updatemu(data,out[i,2],mumu,sigma2mu)
        # get random draw for mu
        out[i,1] = rnorm(1,post[1],sqrt(post[2]))
        }
    return (out)
    }

hpd = function(x, prob = 0.95, precision = 1000){
    range = seq(0, 1-prob, length=precision)
    range = cbind(range, range+prob)
    best = range[which.min(apply(range, 1, function(y)
        diff(quantile(x, y)))),]
    return (quantile(x, best))
    }

# get the closest values in density() for the intervals
bounds = function(int, dens, x=TRUE){
    if (x) # returns x-value from density
        return(dens$x[which.min(abs(dens$x-int))])
    return(dens$y[which.min(abs(dens$x-int))])
    }

### Miles per gallon for my Saturn.  Assumed to be normally distributed.
### The prior mean is distribution as normally and the prior variance
### follows an inverse gamma distribution.  Thus, the mean will indicate
### the miles per gallon.

### Miles and Gallons Data from my Saturn.  Includes both highway and
### surface street miles.
mpg.dat=c(209.5/6.171, 296.0/9.664, 240.5/9.867, 221.8/10.231,
        319.3/10.404, 287.5/9.543, 307.3/9.911, 227.9/9.405,
        318.1/9.812, 309.4/9.634, 269.2/9.271, 297.9/10.334,
        163.3/9.913, 300.0/10.12, 355.1/10.384, 286.6/9.611,
        330.6/10.390, 301.0/9.956, 316.2/10.020, 366.9/9.942,
        262.1/10.055, 215.8/10.048, 283.1/8.690, 357.6/9.879,
        242.5/10.121, 188.9/8.485, 311.3/9.870, 225.3/10.236,
        259.8/10.304, 264.8/9.904, 277.3/10.465, 272.5/11.277)
source("~/files/R/dates.R")
dates = dates[41135:41835,]
dates$Day = as.numeric(dates$Day)
days = data.frame("Day"=c("15", "20", "20", " 9", "10", "21", " 7",
    "22", "22", "30", "30", " 9", "26", "24", "29", "10", "12", "24",
    "26", "26", " 1", "12", "16", "26", " 5", "23", " 6", "18"),
    "Month"=c("August", "August", "September", "November", "November",
    "November", "December", "December", "December", "December",
    "December", "January", "March", "April", "April", "May", "May",
    "May", "May", "May", "July", "August", "August", "August",
    "October", "November", "January", "March"),
    "Year"=c(2012, 2012, 2012, 2012, 2012, 2012, 2012, 2012,
    2012, 2012, 2012, 2013, 2013, 2013, 2013, 2013, 2013, 2013,
    2013, 2013, 2013, 2013, 2013, 2013, 2013, 2013, 2014, 2014))
n = length(mpg.dat)
indexed = numeric(n)
for (i in 1:n)
    indexed[i] = which(apply(dates[,2:4], 1, function(x) all(days[i,] == x)))

# time plot
plot(indexed, mpg.dat, pch=20, type='b')

### My prior beliefs about the data, which would be normal
mpg.mumu = 30
mpg.s2mu = 10
mpg.prsh = 20
mpg.prsc = 0.002

mpgpost = gibbs(mpg.dat, mpg.mumu, mpg.s2mu, mpg.prsh, mpg.prsc, 100000)

### Prior and Poster Plots of Mean and Variance (who cares?)
# plot(density(mpgpost[,1]),lwd=2,col='red',xlim=c(20,40))
# curve(dnorm(x,mpg.mumu,sqrt(mpg.s2mu)),lwd=2,add=T)
# plot(density(mpgpost[,2]),lwd=2,col='green',xlim=c(10,50))
# curve(igpdf(x,mpg.prsh,mpg.prsc),lwd=2,add=T)

### posterior predictive
preds = rnorm(nrow(mpgpost),  mpgpost[,1], sqrt(mpgpost[,2]))
png("~/files/R/figs/mpghist.png")
par(mar=c(4.1,3.1,3.1,1.1))
hist(mpg.dat, breaks=10, col='gray', main="", xlab="MPG", ylab="",
    xlim=c(min(mpg.dat)-3,max(mpg.dat)+3),freq=FALSE) 
points(density(preds), col='red', type='l', lwd=2)
legend(15, 0.14, c("Posterior Predictive", "Data"), lty=c(1, 1),
    col=c('red', "gray"), lwd=c(2,10), cex=1.3)
dev.off()

### posterior distribution on mu
png("~/files/blog/mpg/mupost.png", width=8, height=8, units="in",
    res=300)
densmu = density(mpgpost[,1])
int = quantile(mpgpost[,1],c(0.025,0.975))
par(mar=c(4.1,3.1,3.1,1.1))
plot(densmu,lwd=2,col='blue', xlab=expression(mu),
    main=expression(paste("Posterior of ", mu)), cex.main=1.3,
    cex.lab=1.3, ylab="")
polygon(density(mpgpost[,1]), col='lightblue', border='blue')
lines(rep(bounds(int[1], densmu), 2), c(0, bounds(int[1],
    densmu, FALSE)), col='green', lwd=2)
lines(rep(bounds(int[2], densmu), 2), c(0, bounds(int[2],
    densmu, FALSE)), col='green', lwd=2)
legend(24.2, 0.42, "95% Probability Interval", lty=1,
    col='green', lwd=2)
dev.off()

### Point Estimates for the Mean Parameter
# These should all be about the same since it's normally distributed
# Mean
mean(mpgpost[,1])
# Median
quantile(mpgpost[,1],0.5)
# Mode
density(mpgpost[,1])$x[which.max(density(mpgpost[,1])$y)]
abline(v=quantile(mpgpost[,1],c(0.025,0.975)), col='green',
    lwd=2)

### Point Estimates for the Variance Parameter
# Mean
mean(mpgpost[,2])
# Median
quantile(mpgpost[,2],0.5)
# Mode
density(mpgpost[,2])$x[which.max(density(mpgpost[,2])$y)]


### Posterior Probability Intervals on the Mean Parameter
# 90%
quantile(mpgpost[,1],c(0.05,0.95))
# 95%
quantile(mpgpost[,1],c(0.025,0.975))
# 99%
quantile(mpgpost[,1],c(0.005,0.995))


### Posterior Probability Intervals on the Variance Parameter
# 90%
quantile(mpgpost[,2],c(0.05,0.95))
# 95%
quantile(mpgpost[,2],c(0.025,0.975))
# 99%
quantile(mpgpost[,2],c(0.005,0.995))

### Mean
### As of 31 Dec 2012, the 99% posterior probability interval is (25.8, 32.8)
### As of 01 Mar 2013, the 99% posterior probability interval is (25.6, 31.9)
### As of 14 May 2013, the 99% posterior probability interval is (26.0, 32.0)



mean(770/mpgpost[,1])*3.5
#    ^-miles to drive ^--price of a gallon of gas
# is about how much the trip will cost

### hypothesis test on the skewness of the data
skew = function(x) mean((x-mean(x))^3) / var(x)^(3/2)

# using the posterior predictive
pred.mean = mean(preds)
pred.sd = sd(preds)
observed.skew = skew(mpg.dat)
null.skew = double(100000)
for (i in 1:length(null.skew))
    null.skew[i] = skew(rnorm(length(mpg.dat), pred.mean, pred.sd))
png("~/files/blog/mpg/simu.png")
par(mar=c(4.1,3.1,3.1,1.1))
hist(null.skew, col='gray', freq=FALSE, main="Simulation-based",
    xlab="Sample Skewness", ylab="", cex.main=1.3, cex.lab=1.3)
abline(v=c(-1,1)*observed.skew, col='green', lwd=2)
dev.off()
mean(null.skew < observed.skew) + mean(null.skew > -observed.skew)

test.skew = double(100000)
for (i in 1:length(null.skew))
    test.skew[i] = skew(rnorm(28, 5, 100))
hist(test.skew, col='gray', freq=FALSE, main="Simulation-based",
    xlab="Sample Skewness", ylab="", cex.main=1.3, cex.lab=1.3)
mean(test.skew < observed.skew) + mean(test.skew > -observed.skew)

# using bootstrapping
boot.skew = double(100000)
for (i in 1:length(boot.skew))
    boot.skew[i] = skew(sample(mpg.dat, replace=TRUE))
png("~/files/blog/mpg/boot.png")
par(mar=c(4.1,3.1,3.1,1.1))
hist(boot.skew, col='gray', freq=FALSE, xlab="Sample Skewness",
    main="Bootstrapping", cex.main=1.3, cex.lab=1.3)
abline(v=quantile(boot.skew, c(0.025, 0.975)), col='green', lwd=2)
dev.off()
quantile(boot.skew, c(0.025, 0.975))

