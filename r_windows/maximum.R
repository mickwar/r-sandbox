weibpdf=function(t,alpha,beta){
	(beta/t)*(t/alpha)^beta*exp(-(t/alpha)^beta)
	}


### maximize parameters with derivatives?

var1=function(x){
	dnorm(x,10,2)
	}

var2=function(y){
	dnorm(y,10,5)
	}

# We want to find the values of x and y that will make var1*var2 the highest
curve(var1(x),from=0,to=10,ylim=c(0,25))
curve(var2(x),from=0,to=10,add=T)
curve(var1(x)*var2(x),from=0,to=10,add=T,lwd=2)


limx=c(4,16)
limy=c(-5,25)
acc=0.1	#accuracy
ranx=seq(limx[1],limx[2],acc)
rany=seq(limy[1],limy[2],acc)
maxfunc=function(x,y){
	var1(x)*var2(y)
	}
out=matrix(NA,length(ranx),length(rany))

x1=rnorm(1,5,3)

maximize=function(func){
	

	}






