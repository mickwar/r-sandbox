spiral = function(nsides, theta, reps = 50){
    rad = pi/180

    total_angle = 180 * (nsides - 2)
    tau = total_angle / nsides

    poly = matrix(0, nsides, 2)
    for (i in 1:nsides)
        poly[i,] = c(cos((180-tau)*(i-1)*rad), sin((180-tau)*(i-1)*rad))

    plot(0, type='n', xlim=c(-1.0, 1.0), ylim=c(-1, 1))
    polygon(poly)

    # create rotational matrix
    rotate = matrix(c(cos(theta*rad), sin(theta*rad), -sin(theta*rad), cos(theta*rad)), 2, 2)

    # compute scale
    scale = sin(tau/2*rad)/ sin((180-theta-tau/2)*rad)

    # iterate downward
    for (i in 1:reps){
        poly = scale * poly %*% rotate
        polygon(poly)
        }
    }

spiral(4, 45, 1000)

x = 8
for (i in 1:(180/x*2-1))
    spiral(x, i, 1000)
