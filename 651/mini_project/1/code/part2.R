### Problem 2

### Poisson with a Gamma prior

path="~/files/R/651/mini_project/1/"

comps=read.table(paste(path,"data/poisson.dat",sep=""))
comps=t(comps)

# The expert could only say that they thought there "ought" to be
# no more than 10 failures in one month.  Poisson distributions have
# a mean of lambda, so we will want to choose lamba to be smaller than
# 10.  Lambda of 4 seems to be an appropriate choice.  When doing a
# random plot, there were very few if any instances when the randomized
# poisson produce a value greater than 10.  Lambda of 5 did produce
# values greater than 10 more oftern, so we want to choose parameters
# for our prior (Gamma) so the mean is about 4.

# After thinking about it more, I'm going to choose 18,3.  This is
# because of the experts lack of confidence by saying "ought", so
# with these priors, we are saying there is more variability, that there
# is a higher probability of more than 10 computers failing in a month.

priorA=18
priorB=3

postA=priorA+sum(comps)
postB=priorB+length(comps)

#pdf(paste(path,"figs/2priorplot.pdf",sep=""),9,5)
#polygon(x=seq(0,12,0.001),
#    y=dgamma(seq(0,12,0.001),priorA,priorB),
#    col='lightgray',border='black')
#dev.off()

pdf(paste(path,"figs/2postplot.pdf",sep=""),9,5)
curve(dgamma(x,postA,postB),n=1001,xlim=c(0,12),
    main="",ylab="",xlab="",ylim=c(0,1.2))
polygon(x=seq(0,12,0.001),
    y=dgamma(seq(0,12,0.001),postA,postB),
    col='lightgray',border='black')
curve(dgamma(x,priorA,priorB),n=1001, lty=2,
    main="",ylab="",xlab="", add=TRUE)
abline(v=mean(comps), col='red')
legend(8, 1.1, c("Prior", "Posterior", "MLE"), lty=c(2,1,1), col=c("black","black","red"), cex=1.5)
dev.off()

### Estimates

#mean
postA/postB

#median
qgamma(0.5,postA,postB)

#mode
(postA-1)/postB

#variance
postA/postB^2

#standard deviation
sqrt(postA/postB^2)

#highest posterior density interval
qgamma(c(0.025, 0.975), postA, postB)

n=10000
p=0.95
ints=matrix(NA,n,2)
for (i in 0:(n-1)){
    l=(1-p)/(n-1)*i
    ints[i+1,1]=qgamma(l,postA,postB)
    ints[i+1,2]=qgamma(l+p,postA,postB)
    }
len=ints[,2]-ints[,1]
pdf(paste(path,"figs/2intervals.pdf",sep=""))
plot(seq(1,n,1)/n*(1-p),len[seq(1,n,1)],ylim=c(0.5,3),type='l',
    main="",
    ylab="Interval Length",
    xlab="Starting Interval Value")
dev.off()
ints[which.min(len),]    #this is the approximation of the narrowest interval
min(len)


# Other priors

pdf(paste(path,"figs/2otherpriors.pdf",sep=""),9, 5)
curve(dgamma(x,16+sum(comps),4+length(comps)),add=FALSE,col='green',xlim=c(0,15),lwd=2,    # more confidence in experts (real likely less than 10)
    main="",xlab="",ylab="")
curve(dgamma(x,8+sum(comps),2+length(comps)),n=1001,add=TRUE,col='blue',lwd=2)            # should have less than 10 failures per month
curve(dgamma(x,8+sum(comps),1+length(comps)),n=1001,add=TRUE,col='red',lwd=2)            # feels experts were too shaky in their certainty of number of failures
curve(dgamma(x,18+sum(comps),3+length(comps)),n=1001,add=TRUE,col='black',lwd=2)        # loosely accepts experts opinion
dev.off()
