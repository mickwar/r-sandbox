### Problem 4

### Normal (unknown mean) with Normal prior

path="~/files/R/651/mini_project/1/"

exammu=read.table(paste(path,"data/normalmean.dat",sep=""))
exammu=t(exammu)

priorA=85
priorB=4
sigma2=81

n=length(exammu)
ybar=mean(exammu)

calc.postA = function(priorA, priorB)
    (n*ybar*priorB+priorA*sigma2)/(n*priorB+sigma2)
calc.postB = function(priorB)
    (sigma2*priorB)/(n*priorB+sigma2)

postA = calc.postA(priorA, priorB)
postB = calc.postB(priorB)

#pdf(paste(path,"figs/4priorplot.pdf",sep=""),9,5)
#xx=seq(0,100,.1)
#curve(dnorm(x,priorA,sqrt(priorB)),from=70,to=100,n=1000,main="",
#    xlab="",ylab="",ylim=c(0,0.4))
#polygon(x=xx,
#    y=dnorm(xx,priorA,sqrt(priorB)),
#    col='lightgray',border='black')
#dev.off()

pdf(paste(path,"figs/4postplot.pdf",sep=""),9,5)
xx=seq(0,100,.1)
curve(dnorm(x,postA,sqrt(postB)),from=75,to=95,n=1000,main="",
    xlab="",ylab="",ylim=c(0,0.32))
polygon(x=xx,
    y=dnorm(xx,postA,sqrt(postB)),
    col='lightgray',border='black')
curve(dnorm(x,priorA,sqrt(priorB)),from=70,to=100,n=1000,main="",
    xlab="",ylab="",ylim=c(0,0.4), add=TRUE, lty=2)
abline(v=ybar, col='red')
legend(90, 0.32, c("Prior", "Posterior", "MLE"), lty=c(2,1,1), col=c("black","black","red"), cex=1.5)
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
pdf(paste(path,"figs/4otherpriors.pdf",sep=""),9,5)
pa = calc.postA(85, 4)
pb = calc.postB(4)
curve(dnorm(x,pa,sqrt(pb)),xlim=c(70,100),n=1000,main="",ylab="",
    ylim=c(0,0.5),type='l',lwd=2)
pa = calc.postA(90, 1)
pb = calc.postB(1)
curve(dnorm(x,pa,sqrt(pb)),type='l',col='red',add=T,lwd=2,n=1000)
pa = calc.postA(75, 25)
pb = calc.postB(25)
curve(dnorm(x,pa,sqrt(pb)),type='l',col='blue',add=T,lwd=2,n=1000)
pa = calc.postA(80, 16)
pb = calc.postB(16)
curve(dnorm(x,pa,sqrt(pb)),type='l',col='green',add=T,lwd=2,n=1000)
dev.off()
