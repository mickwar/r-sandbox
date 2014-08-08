setup.colors=function(detail=50){
	out=matrix(NA,detail^3,3)
	row=0
	for (red in seq(0,1,1/(detail-1))){
		for (blue in seq(0,1,1/(detail-1))){
			for (green in seq(0,1,1/(detail-1))){
				row=row+1
				out[row,1]=red
				out[row,2]=blue
				out[row,3]=green
				}
			}
		}
	return(out)
	}

color.plot=function(palette,fixed.limits=TRUE){
	require(rgl)
	n=dim(palette)[1]
	if (fixed.limits){
		plot3d(palette[1:n,],col=rgb(palette[1:n,1],
						palette[1:n,2],
						palette[1:n,3]),
			xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),zlim=c(-0.01,1.01))
		}
	else {
		plot3d(palette[1:n,],col=rgb(palette[1:n,1],
						palette[1:n,2],
						palette[1:n,3]))
		}
	}
	

colors=setup.colors(100)
color.plot(colors)

### changing Red
for (i in seq(0,1,1/(100-1))){
	plot(colors[colors[,1]==i,2:3],col=rgb(i,
		colors[colors[,1]==i,2],
		colors[colors[,1]==i,3]),
		xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),lwd=2,pch=20,
		xlab="Green",ylab="Blue",main=paste0("Red: ",i))
	}

### changing Green
for (i in seq(0,1,1/(100-1))){
	plot(colors[colors[,2]==i,c(1,3)],col=rgb(
		colors[colors[,2]==i,1],
		colors[colors[,2]==i,2],
		colors[colors[,2]==i,3]),
		xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),lwd=2,pch=20,
		xlab="Red",ylab="Blue",main=paste0("Green: ",i))
	}

### changing Blue
for (i in seq(0,1,1/(100-1))){
	plot(colors[colors[,3]==i,c(1,2)],col=rgb(
		colors[colors[,3]==i,1],
		colors[colors[,3]==i,2],
		colors[colors[,3]==i,3]),
		xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),lwd=2,pch=20,
		xlab="Red",ylab="Green",main=paste0("Blue: ",i))
	}


system.time(timed <- setup.colors(10))
system.time(color.plot(timed))

system.time(timed <- setup.colors(20))
system.time(color.plot(timed))

system.time(timed <- setup.colors(30))
system.time(color.plot(timed))

system.time(timed <- setup.colors(40))
system.time(color.plot(timed))

system.time(timed <- setup.colors(50))
system.time(color.plot(timed))