dat = read.table("~/files/data/simplex_1b.txt")

himmelblau = function(x, y)
    (x^2 + y - 11)^2 + (x + y^2 - 7)^2

x = seq(-5, 5, length=100)
y = seq(-5, 5, length=100)
z = matrix(0, length(x), length(y))
for (i in 1:length(x))
    for (j in 1:length(y))
        z[i,j] = himmelblau(x[i], y[j])

m = nrow(dat) / 3
contour(x, y, z, 50)
for (i in 1:m){
    points(c(3, -2.805118, -3.779310, 3.584428),
        c(2, 3.131312, -3.283186, -1.84126), pch=3, col='red', lwd=3)
    if (i > 1){
        j = 3*(i-2)+1
        lines(dat[c(j:(3*i), j), 1], dat[c(j:(3*i), j), 2],col='white')
        }
    j = 3*(i-1)+1
    lines(dat[c(j:(3*i), j), 1], dat[c(j:(3*i), j), 2])
#   readline()
    Sys.sleep(0.1)
    }

#####

dat = read.table("~/files/data/simplex_2a.txt")

banana = function(x, y)
    (1 - x)^2 + 100*(y - x^2)^2

x = seq(-2, 2, length=100)
y = seq(-1, 2, length=100)
z = matrix(0, length(x), length(y))
for (i in 1:length(x))
    for (j in 1:length(y))
        z[i,j] = banana(x[i], y[j])

m = nrow(dat) / 3
contour(x, y, z, 50)
for (i in 1:m){
    points(1, 1, pch=3, col='red', lwd=3)
    if (i > 1){
        j = 3*(i-2)+1
        lines(dat[c(j:(3*i), j), 1], dat[c(j:(3*i), j), 2],col='white')
        }
    j = 3*(i-1)+1
    lines(dat[c(j:(3*i), j), 1], dat[c(j:(3*i), j), 2])
#   readline()
    Sys.sleep(0.1)
    }
