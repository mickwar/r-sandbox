library(rgl)

n = 20
s = 0
xmax = 55.125
ymax = 17.0
zmax = 43.375
amax = max(xmax, ymax, zmax)
oldcol = 'blue'
newcol = 'red'
plot3d(0, 0, 0, xlim = c(0-s, amax+s), ylim = c(0-s, amax+s), zlim = c(0-s, amax+s),
    type='n', ax = FALSE, xlab = "x", ylab = "y", zlab = "z")

#a = rbind(cbind(seq(0, 1.375, length = n), 0, 0),
#cbind(seq(0, 1.375, length = n), 17, 0),
#cbind(0, seq(0, 17, length = n), 0),
#cbind(1.375, seq(0, 17, length = n), 0))
#
#lines3d(a, type='shade')
#shade3d(a)

lines3d(col=oldcol,seq(0, 1.375, length = n), 0, 0)
lines3d(col=oldcol,seq(0, 1.375, length = n), 17, 0)
lines3d(col=oldcol,0, seq(0, 17, length = n), 0)
lines3d(col=oldcol,1.375, seq(0, 17, length = n), 0)
lines3d(col=oldcol,0, 0, seq(0, 30, length = n))
lines3d(col=oldcol,1.375, 0, seq(0, 28.625, length = n))
lines3d(col=oldcol,0, 17,  seq(0, 30, length = n))
lines3d(col=oldcol,1.375, 17, seq(0, 28.625, length = n))

lines3d(col=oldcol,seq(xmax-1.375, xmax, length = n), 0, 0)
lines3d(col=oldcol,seq(xmax, xmax-1.375, length = n), 17, 0)
lines3d(col=oldcol,xmax, seq(0, 17, length = n), 0)
lines3d(col=oldcol,xmax-1.375, seq(0, 17, length = n), 0)
lines3d(col=oldcol,xmax, 0, seq(0, 30, length = n))
lines3d(col=oldcol,xmax-1.375, 0, seq(0, 28.625, length = n))
lines3d(col=oldcol,xmax, 17,  seq(0, 30, length = n))
lines3d(col=oldcol,xmax-1.375, 17, seq(0, 28.625, length = n))

lines3d(col=oldcol,seq(0, xmax, length = n), 0, 30)
lines3d(col=oldcol,seq(0, xmax, length = n), 0, 28.625)
lines3d(col=oldcol,seq(0, xmax, length = n), ymax, 30)
lines3d(col=oldcol,seq(0, xmax, length = n), ymax, 28.625)
lines3d(col=oldcol,0, seq(0, ymax, length = n), 30)
lines3d(col=oldcol,0, seq(0, ymax, length = n), 28.625)
lines3d(col=oldcol,xmax, seq(0, ymax, length = n), 30)
lines3d(col=oldcol,xmax, seq(0, ymax, length = n), 28.625)
lines3d(col=oldcol,1.325, seq(0, ymax, length = n), 28.625)
lines3d(col=oldcol,xmax-1.325, seq(0, ymax, length = n), 28.625)

lines3d(col=newcol,0, ymax, seq(30, zmax, length=n))
lines3d(col=newcol,1.325, ymax, seq(30, 42, length=n))
lines3d(col=newcol,0, ymax-10.5, seq(30, zmax, length=n))
lines3d(col=newcol,1.325, ymax-10.5, seq(30, 42, length=n))

lines3d(col=newcol,1.325, ymax-1.325, seq(30, 42, length=n))
lines3d(col=newcol,0, ymax-1.325, seq(30, 42, length=n))
lines3d(col=newcol,0, ymax-10.5+1.325, c(30, 42))
lines3d(col=newcol,1.325, ymax-10.5+1.325, c(30, 42))
lines3d(col=newcol,c(0, 1.325, 1.325, 0), c(ymax-10.5, ymax-10.5, ymax-10.5+1.325, ymax-10.5+1.325), 30)
lines3d(col=newcol,c(0, 1.325, 1.325, 0), c(ymax, ymax, ymax-1.325, ymax-1.325), 30)

lines3d(col=newcol,34, ymax, seq(30, zmax, length=n))
lines3d(col=newcol,34-1.325, ymax, seq(30, 42, length=n))
lines3d(col=newcol,34, ymax-10.5, seq(30, zmax, length=n))
lines3d(col=newcol,34-1.325, ymax-10.5, seq(30, 42, length=n))

lines3d(col=newcol,34-1.325, ymax-1.325, seq(30, 42, length=n))
lines3d(col=newcol,34, ymax-1.325, seq(30, 42, length=n))
lines3d(col=newcol,34, ymax-10.5+1.325, c(30, 42))
lines3d(col=newcol,34-1.325, ymax-10.5+1.325, c(30, 42))
lines3d(col=newcol,c(34, 34-1.325, 34-1.325, 34,34), c(ymax-10.5, ymax-10.5, ymax-10.5+1.325, ymax-10.5+1.325,ymax-10.5), 30)
lines3d(col=newcol,c(34, 34-1.325, 34-1.325, 34,34), c(ymax, ymax, ymax-1.325, ymax-1.325,ymax), 30)

lines3d(col=newcol,seq(0, 34, length = n), ymax -10.5, 42)
lines3d(col=newcol,seq(0, 34, length = n), ymax -10.5, zmax)
lines3d(col=newcol,seq(0, 34, length = n), ymax, 42)
lines3d(col=newcol,seq(0, 34, length = n), ymax, zmax)

lines3d(col=newcol,0, seq(ymax - 10.5, ymax, length = n), zmax)
lines3d(col=newcol,0, seq(ymax - 10.5, ymax, length = n), 42)
lines3d(col=newcol,34, seq(ymax - 10.5, ymax, length = n), zmax)
lines3d(col=newcol,34, seq(ymax - 10.5, ymax, length = n), 42)
lines3d(col=newcol,34-1.325, seq(ymax - 10.5, ymax, length = n), 42)
lines3d(col=newcol,1.325, seq(ymax - 10.5, ymax, length = n), 42)
