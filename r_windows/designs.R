### Stratified Sampling

require(akima)
require(rgl)
library(gregmisc)
library(randtoolbox)

strat.field=function(size,dimensions,domain,mode){
	if (missing("domain")){
		domain=list()
		for (d in 1:dimensions)
			domain[[d]]=c(0,1)
		}
	if (!is.list(domain))
		domain=list(domain)
	if (len
	if (length(domain) != dimensions)
		stop("Invalid 'domain' size.")
	for (d in 1:dimensions){
		if (domain[[d]][1] >= domain[[d]][2])
			stop("Invalid 'domain' range.")
		}
	if (missing("mode")){
		mode="random"
	} else {
		set=c("random","center","sobol")
		if (!any(mode==set))
			mode="random"
		}
	if (mode != "sobol"){
		out=matrix(NA,size^dimensions,dimensions)
		for (d in 1:dimensions){
			for (i in 1:(size^dimensions)){
				step=(domain[[d]][2]-domain[[d]][1])/size
				adjusted=ceiling(i/size^(dimensions-d)-1)%%size
				if (mode=="random"){
					out[i,d]=runif(1,
						domain[[d]][1]+step*adjusted,
						domain[[d]][1]+step*(adjusted+1))
					}
				if (mode=="center"){
					out[i,d]=domain[[d]][1] + step*adjusted + step/2
					}
				}
			}
	} else {
		require(randtoolbox)
		print("Note: sobol mode requires domain [0,1]^D")
		out=sobol(size,dimensions)
		}
	return(out)
	}

p=strat.field(10,2,c(0,1),"center")
plot(p,xlim=c(min(p[,1],p[,2])-0.05,max(p[,1],p[,2])+0.05),
	ylim=c(min(p[,1],p[,2])-0.05,max(p[,1],p[,2])+0.05))

gridpoints.3d=function(size,domain=c(0,1),acc=100){
	nlines=3*(size-1)^2+6*(2*(size-1))  # n middle lines + n face lines
	step=(domain[2]-domain[1])/size
	out=matrix(NA,nlines*acc,3)
	# Face lines
	x=rescale(1:acc,domain[1],domain[2])
	r=0
	for (face in 1:6){
		for (side in 1:2){
			for (line in 1:(size-1)){
				for (row in 1:acc){
					r=r+1
					out[r,permutations(3,3)[face,1]]=x[row]
					out[r,permutations(3,3)[face,2]]=domain[side]-line*step*(-1)^side
					out[r,permutations(3,3)[face,3]]=domain[side]
					}
				}
			}
		}
	# Middle Section
	for (face in 1:3){
		for (line in 1:(size-1)^2){
			for (row in 1:acc){
				r=r+1
				out[r,permutations(3,3)[face,1]]=x[row]
				out[r,permutations(3,3)[face,2]]=
					domain[1]+step*permutations((size-1),2,repeats.allowed=T)[line,1]
				out[r,permutations(3,3)[face,3]]=
					domain[1]+step*permutations((size-1),2,repeats.allowed=T)[line,2]
				}
			}
		}
	return(out)
	}

rescale=function(x,newmin,newmax){
	(x-min(x))*(newmax-newmin)/(max(x)-min(x))+newmin
	}

size=3
dim=4
min=3
max=7
stratified=strat.field(size,dim,c(min,max))
grid=NULL
grid=gridpoints.3d(dim,c(min,max),50)

#plot(stratified,xlim=c(min-0.05,max+0.05),ylim=c(min-0.05,max+0.05),pch=20)
#other=matrix(runif(dim*size^dim,min,max),size^dim,dim)

plot3d(stratified,xlim=c(min-0.05,max+0.05),
	ylim=c(min-0.05,max+0.05),
	zlim=c(min-0.05,max+0.05))
plot3d(grid,col='gray',add=T,lwd=0.1)



points(other)
abline(v=seq(min,max,(max-min)/size),col='gray')
abline(h=seq(min,max,(max-min)/size),col='gray')



### Latin hypercude designs

latin=function(size){
	out=matrix(1,size,size)
	out[,1]=sample(size,size)
	for (i in 2:size){
		while (any(out[,i]==out[,1:(i-1)])){
			out[,i]=sample(size,size)
			}
		}
	return(out)
	}

lat.2.coor=function(latin){
	nx=dim(latin)[1]
	ny=dim(latin)[2]
	out=matrix(0,nx*ny,3)
	at=0
	for (row in 1:nx){
		for (col in 1:ny){
			at=at+1
			out[at,1]=row
			out[at,2]=col
			out[at,3]=latin[row,col]
			}
		}
	return(out)
	}

fin = lat.2.coor(latin(4))
plot3d(fin,xlim=c(0.5,4.5),ylim=c(0.5,4.5),zlim=c(0.5,4.5),lwd=3)

perm=function(size,dimensions=2){
	out_mat=matrix(0,size,size)
	for (i in 1:size){
		n=sample(size+1-i,1)
		j=1
		while (j <= n){
			if (any(out_mat[j,1:i]==1)){
				n=n+1
				}
			j=j+1
			}
		out_mat[n,i]=1
		}
	out_vec=which(out_mat==1,arr.ind=TRUE,useNames=FALSE)
	return(list(out_mat,out_vec))
	}

ortho=function(vec){
	if (!is.null(dim(vec))){
		vec=vec[,1]
		}
	out=matrix(0,length(vec),length(vec))
	for (i in 1:length(vec)){
		out[vec[i],i]=1
		}
	return(out)
	}


### Distance: Maximin (maximizing the minimum distance)
maximin=function(default=0,n,dimension=2,domain=c(0,1),
	maxdif=sqrt(n)/100, showplot=TRUE, sobol=TRUE, track.weight=3){
	d=dimension
	out=default
	if (!is.matrix(default)){
		# initial scatter
		if (sobol){
			require(randtoolbox)
			out=sobol(n,d)
			}
		if (!sobol)
			out=matrix(runif(n*d),n,d)
		if (d==1)
			out=matrix(out,n,1)
		}
	all.mu = out
	all.sig = double(n)-1
	all.sig = matrix(rep(-1,n*d),n,d)
	if (showplot){
		if (d==1)
			plot(out,integer(n),xlim=c(-0.05,1.05),ylim=c(-0.5,0.5))
		if (d==2)
			plot(out,xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
		if (d > 2)
			print ("Plots not available for dimensions higher than 2.")
		}
	#while ( !all(abs(out - all.mu)<=all.sig) ){
	niter = 0
	while (niter < 10000){
		niter = niter+1
		pair.index=coor.out(out)
		pair.dist=min(dist(out))
		pair.index=sample(pair.index)
		direction=out[pair.index[1],]-out[pair.index[2],]
		scalar=runif(1,0,maxdif)
		for (i in 1:d){
			out[pair.index[1],i]=out[pair.index[1],i]+scalar*direction[i]
			out[pair.index[1],i]=stay.in.bounds(out[pair.index[1],i])
			all.mu[pair.index[1],i]=
				(track.weight*all.mu[pair.index[1],i]+
					out[pair.index[1],i])/(track.weight+1)
			all.sig[pair.index[1],]=sd(c(all.mu[pair.index[1],i],
				all.mu[pair.index[1],i],out[pair.index[1],i]))
			}
		if (showplot){
			if (d==1){
				plot(out,integer(n),xlim=c(-0.05,1.05),ylim=c(-0.5,0.5))
				points(out[pair.index,],c(0,0),col='red',pch=20)
				}
			if (d==2){
				plot(out,xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
				points(out[pair.index,],col='red',pch=20)
				}
			}
		#print(abs(out - all.mu)<all.sig)
		#print(all.sig)
		}
	return(out)
	}
### Gets the two index of the set of points that 
### have the smallest distance
coor.out=function(mat){
	n=dim(mat)[1]
	point=c(2,1)
	x=which.min(dist(mat))
	while (x>n-point[2]){
		x=x-(n-point[2])
		point=point+1
		}
	point[1]=(point[1]:n)[x]
	return(point)
	}
stay.in.bounds=function(x,bounds=c(0,1)){
	if (x < bounds[1]){
		x = bounds[1]
		}
	if (x > bounds[2]){
		x = bounds[2]
		}
	return(x)
	}
difference=function(x){
	n=length(x)
	out=double(n)
	for (i in 1:n){
		out[i]=x[i]-min(x)
		}
	return(out)
	}

design = matrix(runif(100), 50, 2)
plot(design)
out = maximin(0,9,2,maxdif=0.01,sobol=TRUE,track.weight=3,showplot=FALSE)
points(out, col='red',pch=20)

##################################
### choosing a near-maximin subset
circle = function(h,k,r){	# for 2-d
	out = matrix(0, 200, 2)
	out[1:100,1] = seq(-r, r, length=100)
	out[101:200,1] = seq(r, -r, length=100)
	out[,2] = sqrt(r^2 - out[,1]^2)
	out[1:100,2] = -out[1:100,2]
	out[,1] = out[,1]+h
	out[,2] = out[,2]+k
	return(out)
	}

n = 500
k = 2
m = 17

design.set = matrix(runif(n*k), n, k)
test.set = maximin(0,m,k,showplot=FALSE)
radius = min(dist(test.set))/2

plot(design.set, xlim=c(0,1), ylim=c(0,1))
points(test.set, col='red',pch=20)

dist2 = function(x, y){
	A = t(y)-x
	out = double(nrow(y))
	for (i in 1:nrow(y))
		out[i] = sqrt(t(A[,i])%*%A[,i])
	return (out)
	}

TEMP.D = design.set
new = matrix(NA, m, k)
for (i in 1:m){
	readline()
	dists = dist2(test.set[i,], TEMP.D)
	index = which(dists<radius)
	if (length(index) > 0){
		new[i,] = TEMP.D[which.min(dists),]
		white = TEMP.D[index,]
		TEMP.D = matrix(TEMP.D[-index,],ncol=k)
		}
	points(new, col='green',pch=20)
	points(white, col='white')
	circ = circle(test.set[i,1], test.set[i,2], radius)
	points(circ, col='green',type='l')
	}

points(new, col='green', pch=20, xlim=c(0,1),ylim=c(0,1))


# For minimax designs I think you need to choose points to
# calculate the distances used instead of doing it point by point

(out=maximin(0,6,2,iter=3000,maxdif=0.5, acc=-1))
plot(out,type='l')
out2=means(out)
points(out2,type='l',col='red')





out=maximin(0,16,2,iter=100000,maxdif=2, acc=-1)

means=function(values, range=50){
	end=floor(length(values)/range)
	if (end != length(values)/range){
		out=matrix(0,end+1,2)
		out[end+1,1]=length(values)
		out[end+1,2]=mean(values[(1+end*range):length(values)])
	} else {
		out=matrix(0,end,2)
		}
	for (i in 1:end){
		out[i,1]=range*i
		out[i,2]=mean(values[(1+range*(i-1)):(range*i)])
		}
	return(out)
	}

win=matrix(c(
0.5, 1,
0, 0,
1, 0),3,2,byrow=T)

okay=matrix(c(0.792544, 0.0000000,
0.000000, 0.6415444,
1.000000, 1.0000000),3,2,byrow=T)

alright=matrix(c(
0.5965028, 0.43312888,
0.0000000, 0.00000000,
0.6456175, 1.00000000,
1.0000000, 0.35244739,
0.9703870, 0.75771483,
0.0000000, 0.99282493,
0.1932630, 0.36523902,
0.8048835, 0.00000000,
0.4050609, 0.01293363,
0.3252324, 0.75266533),10,2,byrow=T)

while( readline() ){
	1+1}


fine=maximin(0,10)

n=2
out=maximin(strat.field(n,2),n^2,iter=10000)
plot(out[1:10000],type='b')
abline(h=(4+2*sqrt(2)),lty=2)
