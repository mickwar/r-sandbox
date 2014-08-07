calc.mode.d = function(x, precision=512, method="density"){
    d = density(x, n = precision)
    dx = d$x
    dy = d$y
    dx[which.max(dy)]
    }
calc.mode.l = function(x, precision=512, method="loess"){
    d = density(x, n = precision)
    dx = d$x
    dy = d$y
    l = loess(dy ~ dx)
    l$y[which.max(l$x)]
    }
calc.mode.old = function(x, precision=512)
    density(x, n=precision)$x[which.max(density(x, n=precision)$y)]

# calculating the mode using loess seems to give the most
# accurate results, but takes considerably longer when
# using higher precision on the density (i.e., using high
# n in density()), recommend using the loess version
# with density(x, n = 512), the default

x = rnorm(10000)
ns = seq(100, 10000, by=100)
mode.out = matrix(0, nrow=3, ncol=length(ns))
time.out = matrix(0, nrow=3, ncol=length(ns))
for (i in 1:length(ns)){
    time.out[1, i] = system.time(mode.out[1, i] <- calc.mode.d(x, ns[i]))[3]
    time.out[2, i] = system.time(mode.out[2, i] <- calc.mode.l(x, ns[i]))[3]
    time.out[3, i] = system.time(mode.out[3, i] <- calc.mode.old(x, ns[i]))[3]
    }

par(mfrow=c(2,1), mar=double(4)+2, oma=double(4)+0.5)
plot(mode.out[1,], type='l', ylim=c(min(mode.out), max(mode.out)))
points(mode.out[2,], type='l', col="red")
points(mode.out[3,], type='l', col="blue")

plot(time.out[1,], type='l', ylim=c(min(time.out),max(time.out)))
points(time.out[2,], type='l', col="red")
points(time.out[3,], type='l', col="blue")
