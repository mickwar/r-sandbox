### Problem 5

### Normal (unknown variance) with an Inverse Gamma prior

path="~/files/R/651/mini_project/1/"

examvar=read.table(paste(path,"data/normalvariance.dat",sep=""))
examvar=t(examvar)

igpdf=function(x,alpha,beta)
    (beta^alpha)/(gamma(alpha))*(1/x)^(alpha+1)*exp(-beta/x)


ig2pdf=function(x,alpha,beta){
    out=beta
    toA=1
    toB=0
    for (a in 1:(alpha-1)){
        out=out/(alpha-a)
        if (toB<ceiling(1*beta/3)){
            out=out*exp(-1/x)
            toB=toB+1
            }
        if (toA<ceiling(1*alpha/3)){
            out=out*beta
            toA=toA+1
            }
        }
    for (a in 1:(alpha+1)){
        out=out*(1/x)
        if (toB<ceiling(2*beta/3)){
            out=out*exp(-1/x)
            toB=toB+1
            }
        if (toA<ceiling(2*alpha/3)){
            out=out*beta
            toA=toA+1
            }
        }
    for (ab in 1:max(alpha,beta)){
        if (toB<ceiling(beta)){
            out=out*exp(-1/x)
            toB=toB+1
            }
        if (toA<ceiling(alpha)){
            out=out*beta
            toA=toA+1
            }
        }
    return(out)
    }


### the gamma function gets too big so it freaks out,
### try a metropolis sampler

priorA=2
priorB=50
mu=87

postA=priorA+length(examvar)/2
postB=priorB+0.5*sum((examvar-mu)^2)

#png(paste(path,"5priorplot.png",sep=""),width=720,height=400)
#xx=seq(0,300,.1)
#plot(xx,igpdf(xx,priorA,priorB),xlim=c(0,300),
#    main="",ylab="",xlab="",ylim=c(0,0.04),type='l')
#polygon(x=c(xx,300),
#    y=c(igpdf(xx,priorA,priorB),0),
#    col='lightgray',border='black')
#dev.off()

mle = sum((examvar-mu)^2) / length(examvar)

pdf(paste(path,"figs/5postplot.pdf",sep=""),9,5)
xx=seq(0,300,.1)
plot(xx,igpdf(xx,postA,postB),xlim=c(0,300),
    main="",ylab="",xlab="",ylim=c(0,0.028),type='l')
polygon(x=c(xx,300),
    y=c(igpdf(xx,postA,postB),0),
    col='lightgray',border='black')
points(xx,igpdf(xx,priorA,priorB),xlim=c(0,300),
    main="",ylab="",xlab="",ylim=c(0,0.04),type='l',lty=2)
abline(v=mle, col='red')
legend(220, 0.028, c("Prior", "Posterior", "MLE"), lty=c(2,1,1), col=c("black","black","red"), cex=1.5)
dev.off()

### Estimates
samp=1/rgamma(10000,postA,postB)

#mean
postB/(postA-1)

#median
quantile(samp,0.5)

#mode
postB/(postA+1)

#variance
(postB^2)/(((postA-1)^2)*(postA-2))

#standard deviation
sqrt((postB^2)/(((postA-1)^2)*(postA-2)))

# central 95%
quantile(samp, c(0.025, 0.975))

#highest posterior density interval
n=5000
p=0.95
ints=matrix(NA,n,2)
for (i in 0:(n-1)){
    l=(1-p)/(n-1)*i
    ints[i+1,1]=quantile(samp,l)
    ints[i+1,2]=quantile(samp,l+p)
    }
len=ints[,2]-ints[,1]
ints[which.min(len),]
min(len)

# Other priors


add.A = length(examvar)/2
add.B = 0.5*sum((examvar-mu)^2)

pdf(paste(path,"figs/5otherpriors.pdf",sep=""),9,5)
xx=seq(0,500,.1)
plot(xx,igpdf(xx,2+add.A,50+add.B),xlim=c(0,200),
    main="",ylab="",xlab="",ylim=c(0,0.03),type='l',lwd=2)
lines(xx,igpdf(xx,2+add.A,150+add.B),lwd=2,col='red')
lines(xx,igpdf(xx,4+add.A,60+add.B),lwd=2,col='blue')
lines(xx,igpdf(xx,5+add.A,120+add.B),lwd=2,col='green')
dev.off()
