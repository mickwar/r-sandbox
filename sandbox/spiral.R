spiral = function(nsides, theta, thresh = 1e-3, xlim = c(-0.5, 0.5), ylim=c(-0.5, 0.5)){
    rad = pi/180

    total_angle = 180 * (nsides - 2)
    tau = total_angle / nsides

    poly = matrix(0, nsides, 2)
    for (i in 1:nsides)
        poly[i,] = c(cos((180-tau)*(i-1)*rad), sin((180-tau)*(i-1)*rad))

    plot(0, type='n', xlim=xlim, ylim=ylim)
    polygon(poly)

    # create rotational matrix
    rotate = matrix(c(cos(theta*rad), sin(theta*rad), -sin(theta*rad), cos(theta*rad)), 2, 2)

    # compute scale
    scale = sin(tau/2*rad)/ sin((180-theta-tau/2)*rad)
    
    len = 1

    # iterate downward
    while (len > thresh){
        len = len * scale
        poly = scale * poly %*% rotate
        polygon(poly)
        }
    }

spiral(3, 60)
spiral(4, 45)
spiral(5, 40)
spiral(6, 30)

x = 4
for (i in 1:(180/x*2-1))
    spiral(x, i, 1e-3, c(-0.3,0.3), c(-0.3,0.3))


setEPS()
postscript("../figs/spiral_4b.eps")
pdf("../figs/spiral_4b.pdf")
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
spiral(4, 12.0, 1e-3, c(-0.35,0.35), c(-0.35,0.35))
dev.off()
