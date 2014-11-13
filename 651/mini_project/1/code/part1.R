### Problem 1

### Binomial Likelihood with a Beta Prior

path="~/files/R/651/mini_project/1/"

bonds=read.table(paste(path,"data/binomial2.dat",sep=""))
#V1 is date
#V2 is number of at-bats in the game (n)
#V3 is number of homeruns in the game (y)

# Choosing the priors.  From the previous year (2000), Bonds had 49 homeruns
# with 480 at-bats.  His percentage of homeruns that year is 11%.  However, this
# percentage has been increasing since he began his professional career in 1986.
# So we will want to choose priors that produce a mean of about 10%.  Since
# Bonds appears to have been improving during his career, we need the variance
# to be high enough to take into account this change of homerun rates.

priorA=5
priorB=45

n=sum(bonds$V2)
y=sum(bonds$V3)
postA=priorA+y
postB=priorB+n-y

pdf(paste(path,"figs/1postplot.pdf",sep=""), 9, 5)
curve(dbeta(x,postA,postB),n=1001,main="",
    ylab="",xlab="", from=0, to=0.4)
polygon(x=seq(0,0.4,0.001),
    y=dbeta(seq(0,0.4,0.001),postA,postB),
    col='darkgray',border='black')
curve(dbeta(x,priorA,priorB),n=1001,ylim=c(0,25),
    main="", add=TRUE,
    ylab="",xlab="", lty=2)
abline(v=y/n, col='red')
legend(0.3, 25, c("Prior", "Posterior", "MLE"), lty=c(2,1,1), col=c("black","black","red"), cex=1.5)
dev.off()

### Estimates

#mean
postA/(postA+postB)

#median
qbeta(0.5,postA,postB)

#mode
(postA-1)/(postA+postB-2)

#variance
(postA*postB)/((postA+postB)^2*(postA+postB+1))

#standard deviation
sqrt((postA*postB)/((postA+postB)^2*(postA+postB+1)))

# central 95% interval
qbeta(c(0.025, 0.975), postA, postB)

# minimum interval approximation (HPD)
n=10000
p=0.95
ints=matrix(NA,n,2)
for (i in 0:(n-1)){
    l=(1-p)/(n-1)*i
    ints[i+1,1]=qbeta(l,postA,postB)
    ints[i+1,2]=qbeta(l+p,postA,postB)
    }
len=ints[,2]-ints[,1]
pdf(paste(path,"figs/1intervals.pdf",sep=""), 9, 5)
plot(seq(1,n,1)/n*(1-p),len[seq(1,n,1)],ylim=c(0,0.1),type='l',
    main="",
    ylab="Interval Length",
    xlab="Starting Interval Value")
dev.off()
ints[which.min(len),]    #this is the approximation of the narrowest interval
min(len)

### Other priors
# mean closer to that of all previous seasons, lighter weights
a1=1
b1=14

# only takes into account when bonds was with the giants 1993 - 2000, medium weights
a2=4
b2=45

# only data from 2000 season used, stronger weights
a3=10
b3=90

n=sum(bonds$V2)
y=sum(bonds$V3)
pdf(paste(path,"figs/1otherpriors.pdf",sep=""), 9, 5)
curve(dbeta(x,priorA+y,priorB+n-y),lwd=2,xlim=c(0,.30),ylim=c(0,30),
    main="",xlab="",ylab="")
curve(dbeta(x,a1+y,b1+n-y),lwd=2,col='red',add=T)
curve(dbeta(x,a2+y,b2+n-y),lwd=2,col='blue',add=T)
curve(dbeta(x,a3+y,b3+n-y),lwd=3,col='green',add=T)
dev.off()
