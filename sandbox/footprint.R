### Given an equally-spaced grid of locations and values at each point,
### where is the a x b rectangle that maximizes the sum of the values
### for all possible a x b rectangles?

n = 40
m = 48
xy = expand.grid(1:n, 1:m)

set.seed(3)
z = runif(n*m)
plot(xy, col = rgb(z, 0, 1-z, 2*abs(z-0.5)), pch = 15)

plot(xy, col = rgb((1-z)^(1/1), (1-z)^(1/2), (z)^(3/2)), pch = 15, cex=1.9)
plot(xy, col = rgb((1-z), (1-z), (z)), pch = 15, cex=1.9)




# axb rectangles
a = 6
b = 7
off = as.matrix(expand.grid(0:(a-1), 0:(b-1)))
all.rect = rep(list(matrix(0, a*b, 2)), (n-a+1)*(m-b+1))
ind = which(apply(xy, 1, function(x) max(x[1] + a - 1) %in% 1:n) &
    apply(xy, 1, function(x) max(x[2] + b - 1) %in% 1:m))
grid.ind = rep(list(double(a*b)), length(all.rect))
for (i in 1:length(all.rect)){
    all.rect[[i]] = t(matrix(as.numeric(xy[ind[i],]), 2, a*b)) + off
    grid.ind[[i]] = all.rect[[i]][,1] + (all.rect[[i]][,2]-1)*n
    }

#plot(xy, col = rgb(z, 0, 0, z), pch = 15)
agg = double(length(all.rect))
for (i in 1:length(all.rect)){
#   lines(all.rect[[i]], col = i , lwd = 2)
#   lines(xy[sample(grid.ind[[i]]),], col = rainbow(m)[sample(m, 1)] , lwd = 2)
    agg[i] = sum(z[grid.ind[[i]]])      # maximize over the sum
#   agg[i] = sum(log(z[grid.ind[[i]]])) # maximize over the product
    }
foot.print = which.max(agg)
#all.rect[[foot.print]]
#grid.ind[[foot.print]]

border.points = c(1, a, a*b, 1+a*(b-1), 1)
    
plot(xy, col = rgb(z, 0, 1-z, 2*abs(z-0.5)), pch = 15)
lines(all.rect[[foot.print]][border.points,], lwd = 2)


### Colors
library(MASS)
set.seed(1)
u = runif(2000)
x = rbind(mvrnorm(length(u), c(0, 0), matrix(c(1, 0.0, 0.0, 1),2,2))[which(u < 0.4),],
    mvrnorm(length(u), c(4.5, -3.5), matrix(c(1, 0.6, 0.6, 1),2,2))[which(u >= 0.4 & u < 0.8),],
    mvrnorm(length(u), c(3.0, 0.5), matrix(c(1, -0.7, -0.7, 1),2,2))[which(u >= 0.8),])
plot(x, pch=20)

kd = kde2d(x[,1], x[,2], n = 50)
xy = expand.grid(kd$x, kd$y)
z = c(kd$z)
z = (z - min(z)) / diff(range(z))

plot(xy, col = rgb(
    ifelse(z < 1/3, 1, ifelse(z > 2/3, 0, 3*(2/3 - z))),
    ifelse(z < 1/3, 0.5+z*3/2, 3/2*(1-z)),
    ifelse(z < 1/2, 0, 2*z-1)
    ), pch = 15, cex=1.9)

# How the color components  are changing
qq = seq(0, 1, length = 100)
plot(qq, ifelse(qq < 1/3, 1, ifelse(qq > 2/3, 0, 3*(2/3 - qq))), type='l',
    col = 'red')
lines(qq, ifelse(qq < 1/3, 0.5+qq*3/2, 3/2*(1-qq)), col = 'green')
lines(qq, ifelse(qq < 1/2, 0, 2*qq-1), col = 'blue')

