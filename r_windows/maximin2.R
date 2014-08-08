
### Distance: Maximin (maximizing the minimum distance)
maximin=function(default=0,n,dimension=2,domain=c(0,1),
	maxdif=sqrt(n)/100, showplot=TRUE, track=10){
	d=dimension
	out=default
	if (!is.matrix(default)){
		# initial scatter
		out=matrix(runif(n*d),n,d)
		}
	if (showplot){
		if (d==1)
			plot(out,integer(n),xlim=c(-0.05,1.05),ylim=c(-0.5,0.5))
		if (d==2)
			plot(out,xlim=c(-0.05,1.05),ylim=c(-0.05,1.05))
		if (d > 2)
			print ("Plots not available for dimensions higher than 2.")
		}
	all.means=matrix(0,n,d)
	all.sigma=matrix(0,n,d)
	track.mat=matrix(0,n*track,d)
	count=integer(n)
	while ( !all(abs(out - all.means)<=all.sigma) ){
		pair.index=sample(coor.out(out))
		count[pair.index[1]]=count[pair.index[1]]+1
		if (count[pair.index[1]] > track)
			count[pair.index[1]] = 1
		direction=out[pair.index[1],]-out[pair.index[2],]
		scalar=runif(1,0,maxdif)
		for (i in 1:d){
			out[pair.index[1],i]=out[pair.index[1],i]+scalar*direction[i]
			out[pair.index[1],i]=stay.in.bounds(out[pair.index[1],i])
			track.mat[track*(pair.index[1]-1)+count[pair.index[1]],i]=out[pair.index[1],i]
			all.means[pair.index[1],i]=mean(track.mat[seq(track*(pair.index[1]-1)+1,track*pair.index[1],length=track),i],na.rm=T)
			all.sigma[pair.index[1],i]=sd(track.mat[seq(track*(pair.index[1]-1)+1,track*pair.index[1],length=track),i],na.rm=T)
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
		#print(out-all.means)
		#print(all.sigma)
		#print(track.mat)
		#print(abs(out - all.means)<=all.sigma)
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

maximin(0,16,2)