### example from the appendix in space and space-time 
### modeling using process convolutions
library(nlme) # load lme()

# a fake dataset to make the bumps with
n = 30
m = 7

# create sites s
sb = c(1, 10)
s = seq(sb[1], sb[2], length=n)

# create the data y
e1 = rnorm(n, mean=0, sd=0.1)
e2 = cos(s/10*2*pi*4)*0.2
y = sin(s/10*2*pi)+e2+e1
plot(s, y, col="gray70", pch=20)

# locations of support points
w = seq(sb[1]-2, sb[2]+2, length=m)

# width of kernel
sdkern = 2

# create the matrix K
K = matrix(NA, nrow=n, ncol=m)
for (i in 1:m)
    K[,i] = dnorm(s, mean=w[i], sd=sdkern)

# create a data frame to hold th edata
df1 = data.frame(y=y, K=K, sub=1)
df1$sub = as.factor(df1$sub)

# now fit a mixed model using lme
a1 = lme(fixed=y~1, random=list(sub=pdIdent(~K-1)),
    data=df1, na.action=na.omit)

# obtain and plot the fitted values
a1p = as.vector(predict(a1, df1))
lines(s, a1p, lty=1, col='red', lwd=2)

# now a multiscale version
m1 = c(7, 14, 28)   # number of support points at each scale
w = list(NULL)      # list to hold the support locations

# make support locations w
for (i in 1:3)
    w[[i]] = seq(sb[1]-2, sb[2]+2, length=m1[i])
sdkern = c(2, 1, 0.5)   # kernel width by scale

# generate K matricies for each sacle
K1 = matrix(NA, nrow=n, ncol=m1[1])
for (i in 1:m1[1])
    K1[,i] = dnorm(s, mean=w[[1]][i], sd=sdkern[1])
K2 = matrix(NA, nrow=n, ncol=m1[2])
for (i in 1:m1[2])
    K2[,i] = dnorm(s, mean=w[[2]][i], sd=sdkern[2])
K3 = matrix(NA, nrow=n, ncol=m1[3])
for (i in 1:m1[3])
    K3[,i] = dnorm(s, mean=w[[3]][i], sd=sdkern[3])

# create dataframe df2
df2 = data.frame(y=y, K1=K1, K2=K2, K3=K3, sub=1)

# fit mixed effects model
a2 = lme(fixed=y~1, random=list(sub=pdIdent(~K1-1),
    sub=pdIdent(~K2-1), sub=pdIdent(~K3-1)), data=df2,
    na.action=na.omit)

# get predictions
a2p = as.vector(predict(a2, df2))

# plot it
plot(s, y, pch=20, col="gray70")
lines(s, a1p, lty=1, col='red', lwd=2)
lines(s, a2p, lty=2, col='darkblue', lwd=2)


### attempt at the lattice of iid gamma(2,1) r.v. mentioned
### on page 4
n = 30
xx = seq(0, 5, length=n)
d = max(xx) - min(xx)
beyond = 0.1
# for the kernels:
# x is the locations for data points
# y is the locations at which to predict
kern.gaussian = function(x, y, delta=d/16) # gaussian
    exp(-1/(2*delta^2)*(x-y)^2)
kern.window = function(x, y, delta=d/16) # window
    1*(abs(x-y) <= delta)
kern.epanech = function(x, y, delta=d/16) # epanechnikov
    0.75*(1-(abs(x-y)/delta)^2)*(abs(x-y) < delta)
kern.bisq = function(x, y, delta=d/16) # bisquare
    (1-(abs(x-y)/delta)^2)^2*(abs(x-y) < delta)
kern.tric = function(x, y, delta=d/16) # tricube
    (1-(abs(x-y)/delta)^3)^3*(abs(x-y) < delta)
kern.wendl = function(x, y, delta=d/4){ # wendland
    d=abs(x-y)/delta
    (1-d)^6*(35*d^2+18*d+3)/3*(d < delta)
    }

kern1 = function(x, y) kern.gaussian(x, y, d/16)
kern2 = function(x, y) kern.bisq(x, y, d/8)

yy = rgamma(n, 2, 1)
yy = rexp(n, 1)
plot(xx, yy, type='h', xlim=c(min(xx)-beyond*d, max(xx)+beyond*d))

m = 4*n
smox1 = seq(min(xx)-beyond*d, max(xx)+beyond*d, length=m)
smoothed1 = double(m)
for (i in 1:m)
    smoothed1[i] = sum(yy*kern1(xx, smox1[i]))/
        sum(kern1(xx, smox1[i]))
points(smox1, smoothed1, type='l', col='red')

smox2 = seq(min(xx)-beyond*d, max(xx)+beyond*d, length=m)
smoothed2 = double(m)
for (i in 1:m)
    smoothed2[i] = sum(yy*kern2(xx, smox2[i]))/
        sum(kern2(xx, smox2[i]))
points(smox2, smoothed2, type='l', col='blue')

mean(smoothed1[smox1 >= min(xx) & smox1 <= max(xx)])
mean(smoothed2[smox2 >= min(xx) & smox2 <= max(xx)])


