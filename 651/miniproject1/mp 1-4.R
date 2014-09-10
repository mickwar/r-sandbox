### Problem 4

### Normal (unknown mean) with Normal prior

path="C:/Users/Mickey/Desktop/Downloads/School/2013 A Winter/651/MP 4/"

exammu=read.table(paste(path,"normalmean.dat",sep=""))
exammu=t(exammu)

priorA=85
priorB=4
sigma2=81

n=length(exammu)
ybar=mean(exammu)

postA=(n*ybar*priorB+priorA*sigma2)/(n*priorB+sigma2)
postB=(sigma2*priorB)/(n*priorB+sigma2)

png(paste(path,"4priorplot.png",sep=""),width=720,height=400)
xx=seq(0,100,.1)
curve(dnorm(x,priorA,sqrt(priorB)),from=0,to=100,n=1000,main="",
	xlab="",ylab="",ylim=c(0,0.4))
polygon(x=xx,
	y=dnorm(xx,priorA,sqrt(priorB)),
	col='lightgray',border='black')
dev.off()

png(paste(path,"4postplot.png",sep=""),width=720,height=400)
xx=seq(0,100,.1)
curve(dnorm(x,postA,sqrt(postB)),from=0,to=100,n=1000,main="",
	xlab="",ylab="",ylim=c(0,0.4))
polygon(x=xx,
	y=dnorm(xx,postA,sqrt(postB)),
	col='lightgray',border='black')
dev.off()

#mean/median/mode
postA

#variance
postB

#standard deviation
sqrt(postB)

#highest posterior density interval
qnorm(c(0.025,0.975),postA,postB)

# Other priors

png(paste(path,"4otherpriors.png",sep=""),width=720,height=400)
curve(dnorm(x,85,sqrt(4)),xlim=c(0,100),n=1000,main="",ylab="",
	ylim=c(0,0.4),type='l',lwd=2)
curve(dnorm(x,90,sqrt(1)),type='l',col='red',add=T,lwd=2,n=1000)
curve(dnorm(x,75,sqrt(25)),type='l',col='blue',add=T,lwd=2,n=1000)
curve(dnorm(x,80,sqrt(16)),type='l',col='green',add=T,lwd=2,n=1000)
dev.off()