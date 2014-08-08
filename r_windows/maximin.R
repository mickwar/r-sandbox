
### Distance: Maximin (maximizing the minimum distance)
maximin=function(default=0,n,dimension=2,domain=c(0,1),
	maxdif=sqrt(n)/100, showplot=TRUE, track.weight=3){
	d=dimension
	out=default
	if (!is.matrix(default)){
		# initial scatter
		out=matrix(runif(n*d),n,d)
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
	count=0
	while ( !row.check(abs(out - all.mu)<=all.sig) || count < 0 ){
		count=count+1
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
			all.sig[pair.index[1],]=
				(track.weight*all.sig[pair.index[1],1] + sd(c(all.mu[pair.index[1],i],
				all.mu[pair.index[1],i],out[pair.index[1],i]))) / (track.weight+1)
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
		print(abs(out - all.mu)<all.sig)
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
row.check=function(x){
	n=dim(x)[1]
	test=TRUE
	for (i in 1:n){
		test = all(test,any(x[i,]))
		if (!test)
			break
		}
	return(test)
	}

maximin(0,7,2,maxdif=0.01,track.weight=3)