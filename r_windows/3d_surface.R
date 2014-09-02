library(rgl)

cross = function(...){
	vec = list(...)
	d = length(vec)
	N = double(d+2) + 1
	for (i in 1:d)
		N[i+1] = length(vec[[i]])
	out = matrix(0, prod(N), d)
	for (i in 1:d){
		out[,i] = rep(vec[[i]], times=prod(N[1:i]),
			each=prod(N[(i+2):(d+2)]))
		}
	return(out)
	}

n=100

f=function(x,y){x^2+y^2+x^2*y+4}

f=function(x,y){sqrt(x^2+y^2)}

x=seq(-3,3,length=n)
y=x
#y=seq(-1,1,length=n)

out = cbind(cross(x, y), 0)
for (i in 1:n^2){
	out[i,3]=f(out[i,1],out[i,2])
	}

plot3d(out)
