mw.smooth = function(x, y){
    # gaussian kernel, delta is the bandwidth parameter
    kern = function(x, y, delta=1/24)
        exp(-1/(2*delta^2)*(x-y)^2)
    # m = number of points to do the smoothing,
    # should be greater than length(x) to look nice
    m = 8*length(x)
    outx = seq(min(x), max(x), length=m)
    outy = double(m)
    # loop to compute the weighted value at each new x
    for (i in 1:length(outy))
        outy[i] = sum(y*kern(x, outx[i]))/
            sum(kern(x, outx[i]))
    return (list("x"=outx, "y"=outy))
    }

png("~/files/R/figs/smooth.png")
xx = seq(0, 1, length=15)
yy = rgamma(length(xx), 2, 1)
plot(xx, yy, type='h')
ss = mw.smooth(xx, yy)
lines(ss$x, ss$y, type='l', col='red')
dev.off()

### smoothing a lattice of iid gamma(2,1) r.v. mentioned
### on page 4 of higdon (2002)
n = 30
xx = seq(0, 5, length=n)

# setting d to be the range of the data, related to delta
d = max(xx) - min(xx)

# how to to go beyond the range of the data
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

# for comparing two different kernels
kern1 = function(x, y) kern.gaussian(x, y, d/16)
kern2 = function(x, y) kern.bisq(x, y, d/8)

# yy = rexp(n, 1)
yy = rgamma(n, 2, 1)
plot(xx, yy, type='h', xlim=c(min(xx)-beyond*d, max(xx)+beyond*d))

m = 4*n
# smox = the values of x at which to make predictions
smox1 = seq(min(xx)-beyond*d, max(xx)+beyond*d, length=m)
# smoothed = the predicted values
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

# mean(smoothed1[smox1 >= min(xx) & smox1 <= max(xx)])
# mean(smoothed2[smox2 >= min(xx) & smox2 <= max(xx)])
