# code to generate a mass spectrometry-like output

tangent = function(x, FUN, del=0.5, eps=0.0000001){
    slope = 0
    newslope = Inf
    while (abs(slope - newslope) > eps){
        newslope = slope
        del = del/2
        slope = (FUN(x+del)-FUN(x-del))/(2*del)
        }
    intercept = FUN(x) - slope*x
    return (c(slope, intercept))
    }

xx = seq(0, 2.5, length=100)
yy = seq(0, 2, length=9)
plot(xx, dexp(xx), type='l',ylim=c(0,1))
points(yy, dexp(yy))
abline(tangent(yy[5],dexp)[2], tangent(yy[5],dexp)[1], lty=2)

dat = c(rnorm(800, 80, 0.5),
        rnorm(400, 60, 0.5),
        rnorm(40, 83, 0.7),
        rnorm(100, 98, 0.5),
        rnorm(150, 110, 0.3),
        rnorm(190, 135, 0.2),
        rnorm(30, 15, 0.2),
        rnorm(50, 195, 0.3),
        rnorm(180, 140, 0.3),
        rnorm(15, 80, 10),  # add noise
        rnorm(15, 60, 10),
        rnorm(15, 100, 10),
        rnorm(15, 110, 10),
        rnorm(15, 135, 10),
        rnorm(15, 15, 10),
        rnorm(15, 195, 10),
        rnorm(15, 140, 10))

dat2 = c(rnorm(280, 80, 0.5),
        rnorm(300, 60, 0.5),
        rnorm(60, 83, 0.7),
        rnorm(50, 98, 0.5),
        rnorm(620, 110, 0.3),
        rnorm(390, 135, 0.2),
        rnorm(140, 15, 0.2),
        rnorm(280, 195, 0.3),
        rnorm(180, 140, 0.3),
        rnorm(15, 80, 10),  # add noise
        rnorm(15, 60, 10),
        rnorm(15, 100, 10),
        rnorm(15, 110, 10),
        rnorm(15, 135, 10),
        rnorm(15, 15, 10),
        rnorm(15, 195, 10),
        rnorm(15, 140, 10))

pdf("./mspec1.pdf", width=8, height=5)
par(mfrow=c(1,2), mar=c(5,4,1.5,1), oma=c(0,0,2.5,0), cex.main=1.5,
    cex.sub=1)
hist(dat, xlim=c(0, 200), col='black', breaks=1000,
    xlab = "Mass / Charge", main=expression(paste("At [e]"[0]," = ",x[1])),
    cex.main = 1, ylab = "Relative Abundance", freq=FALSE, ylim=c(0,0.35))
grid(col='darkgray', lty=2)
hist(dat2, xlim=c(0, 200), col='black', breaks=1000,
    xlab = "Mass / Charge", main=expression(paste("At [e]"[0]," = ",x[2])),
    cex.main = 1, ylab = "", freq=FALSE, ylim=c(0,0.35))
grid(col='darkgray', lty=2)
title("Mass Spectrometry Data", outer=T)
dev.off()

# to "match" with tripredAll.pdf
# at [e]=1.5e10
dat = c(rnorm(1140, 45, 0.5),    # blue   (~ 0.75)
        rnorm(165, 60, 0.5),    # red    (~ 0.15)
        rnorm(085, 90, 0.5),     # green  (~ 0.07)
        rnorm(035, 105, 0.5),    # black  (~ 0.03)
        rnorm(020, 115, 0.5),   # orange (~ 0.00)
        rnorm(15, 45, 20),  # add noise
        rnorm(15, 60, 20),
        rnorm(15, 90, 20),
        rnorm(15, 105, 20),
        rnorm(15, 115, 20))

# at [e]=5.0e10
dat2= c(rnorm(060, 45, 0.5),    # blue   (~ 0.06)
        rnorm(270, 60, 0.5),    # red    (~ 0.27)
        rnorm(530, 90, 0.7),     # green  (~ 0.53)
        rnorm(130, 105, 0.5),    # black  (~ 0.14)
        rnorm(005, 115, 0.3),   # orange (~ 0.00)
        rnorm(15, 45, 20),  # add noise
        rnorm(15, 60, 20),
        rnorm(15, 90, 20),
        rnorm(15, 105, 20),
        rnorm(15, 115, 20))


pdf("./mspec2.pdf", width=8, height=5)
par(mfrow=c(1,2), mar=c(5,4,1.5,1), oma=c(0,0,2.5,0), cex.main=1.5,
    cex.sub=1)
hist(dat, xlim=c(0, 200), col='black', breaks=1000,
    xlab = "Mass / Charge", main=expression(paste("At [e]"[0]," = ",x[1])),
    cex.main = 1, ylab = "Relative Abundance", freq=FALSE, ylim=c(0,1.00))
grid(col='darkgray', lty=2)
hist(dat2, xlim=c(0, 200), col='black', breaks=1000,
    xlab = "Mass / Charge", main=expression(paste("At [e]"[0]," = ",x[2])),
    cex.main = 1, ylab = "", freq=FALSE, ylim=c(0,1.00))
grid(col='darkgray', lty=2)
title("Mass Spectrometry Data", outer=T)
dev.off()

