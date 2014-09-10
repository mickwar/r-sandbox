### Problem 3

### Exponential with an Inverse Gamma prior

path="C:/Users/Mickey/Desktop/Downloads/School/2013 A Winter/651/MP 1/"

river=read.table(paste(path,"exponential.dat",sep=""))
river=t(river)

igpdf=function(x,alpha,beta){
	((beta^alpha)/(gamma(alpha)))*(1/x)^(alpha+1)*exp(-beta/x)
}

### since the data make the gamma function and the beta^alpha 
### very large, we quickly reach to either 0 or infinity at any
### value of x.  So "ig2pdf" takes all the multiplications in the pdf
### separately, making the big numbers smaller and preventing
### errors in reaching to 0 or infinity to quickly.
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

priorA=5
priorB=1600

postA=priorA+length(river)
postB=priorB+sum(river)

png(paste(path,"3priorplot.png",sep=""),width=720,height=400)
xx=seq(0,1800,.1)
plot(xx,igpdf(xx,priorA,priorB),xlim=c(0,1500),
	main="",ylab="",xlab="",ylim=c(0,0.005),type='l')
polygon(x=xx,
	y=igpdf(xx,priorA,priorB),
	col='lightgray',border='black')
dev.off()

png(paste(path,"3postplot.png",sep=""),width=720,height=400)
xx=seq(0,1000,1)
plot(xx,ig2pdf(xx,postA,postB),xlim=c(0,1000),
	main="",ylab="",xlab="",ylim=c(0,0.01),type='l')
polygon(x=xx,
	y=ig2pdf(xx,postA,postB),
	col='lightgray',border='black')
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
png(paste(path,"3intervals.png",sep=""))
plot(seq(1,n,1)/n*(1-p),len[seq(1,n,1)],ylim=c(189,200),type='l',
	main="",
	ylab="Interval Length",
	xlab="Starting Interval Value")
dev.off()
ints[which.min(len),]	#this is the approximation of the narrowest interval
min(len)

# Other priors

png(paste(path,"3otherpriors.png",sep=""),width=720,height=400)
xx=seq(0,1500,.1)
plot(xx,igpdf(xx,5,1600),xlim=c(0,1500),		#relatively confident that average length is about 200 miles
	main="",ylab="",xlab="",ylim=c(0,0.005),type='l',lwd=2)
lines(xx,igpdf(xx,2,1000),lwd=2,col='red')	#unsure of what lengths could be, thinks could be higher
lines(xx,igpdf(xx,9,3000),lwd=2,col='blue')
lines(xx,igpdf(xx,6,4000),lwd=2,col='green')
dev.off()