acc=1000
tint=c(0,20*pi)
dots=matrix(0,acc,2)
funx=function(t)
    sin(t)+0.5*sin(5*t)+0.25*cos(2.3*t)
funy=function(t)
    cos(t)+0.5*cos(5*t)+0.25*sin(2.3*t)
dots[,1]=funx(seq(tint[1],tint[2],length=acc))
dots[,2]=funy(seq(tint[1],tint[2],length=acc))
plot(dots[,1],dots[,2],type='l')
plot(dots[,1],type='l')
plot(dots[,2],type='l')

### over time
go=function(dots,n){
    plot(dots[,1],dots[,2],xlim=c(-2,2),ylim=c(-2,2), type='l',
        col='gray')
    for (time in 1:n)
        points(dots[1:time,1],dots[1:time,2],type='l',
            xlim=c(-2,2),ylim=c(-2,2), col='blue')
    }
go2=function(dots,n){
    plot(dots[,1],dots[,2],xlim=c(-2,2),ylim=c(-2,2), type='l',
        col='gray')
    for (time in 2:n)
        points(dots[(time-1):time,1],dots[(time-1):time,2],type='l',
            col='red', lwd=2)

    }

system.time(go(dots,acc))  # ~ 8.2 seconds
system.time(go2(dots,acc)) # ~ 2.4 second
