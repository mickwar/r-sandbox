# original
autotune = function(accept, target = 0.25, k = 2.5){
    k = c((1-1/k)/(cosh(target)-1), (k-1)/(cosh(target-1)-1))
    1+sign(accept-target)*(cosh(accept-target)-1)*
        k[(sign(accept-target)+1)/2+1]
    }

# symmetric
autotune2 = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)

system.time(for (i in 1:500) autotune(xx))
system.time(for (i in 1:500) autotune2(xx))


### I want blue and green to be symmetric, but the upward
### cosh is different that the lower. Even the red (tar=0.5),
### isn't symmetric with it's inverse.
xx = seq(0, 1, length=10000)
tar = 0.5
k = 15
plot(xx,autotune(xx,t=tar,k=k),type='l',col='red',lwd=2)
points(xx,1/autotune(xx,t=tar,k=k),type='l',col='pink2',lwd=2)
points(xx,autotune2(xx,t=tar,k=k),type='l',col='red',lwd=2,lty=2)
points(xx,1/autotune2(xx,t=tar,k=k),type='l',col='pink2',lwd=2,lty=2)
tar = 0.5
k=10
points(xx,autotune(xx,t=tar,k=k),type='l',col='blue',lwd=2)
points(xx,1/autotune(xx,t=tar,k=k),type='l',col='lightblue3',lwd=2)
points(xx,autotune2(xx,t=tar,k=k),type='l',col='blue',lwd=2,lty=2)
points(xx,1/autotune2(xx,t=tar,k=k),type='l',col='lightblue3',lwd=2,lty=2)
tar = 0.5
k=5
points(xx,autotune(xx,t=tar,k=k),type='l',col='green',lwd=2)
points(xx,1/autotune(xx,t=tar,k=k),type='l',col='lightgreen',lwd=2)
points(xx,autotune2(xx,t=tar,k=k),type='l',col='green',lwd=2,lty=2)
points(xx,1/autotune2(xx,t=tar,k=k),type='l',col='lightgreen',lwd=2,lty=2)

# with lower k, there appears to be symmetry, not so with higher

# possible solution: instead of using a separate cosh for x < target,
# use the inverse of when x > target, scaled to fit in [0, targer).
