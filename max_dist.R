rmax = function(n, p, q, a, b, c, d){
    u = runif(n)
    w = runif(sum(u <= p))
    y = ifelse(u <= p, 1, rbeta(sum(u > p), a, b))
    x = ifelse(y != 1, 1, ifelse(w <= q, 1, rbeta(sum(w > q), c, d)))
    return (cbind(x, y))
    }

# Marginals and dot plot
plotmax = function(x){
    par(mfrow = c(2,2), mar = c(0,0,0,0), oma = c(5.1, 4.1, 4.1, 2.1))
    dens = density(x[x [,1] != 1, 1])
    ylim = c(0, max(1, dens$y * mean(x[,1] != 1)))
    plot(1, mean(x[,1] == 1), type ='h', xlim = c(0,1), ylim = ylim,
        axes = FALSE)
    lines(dens$x, dens$y * mean(x[,1] != 1))
    box()
    axis(2)

    plot(0, type='n', axes = FALSE, xlab = "", ylab = "")
    plot(x, pch = 20)

    dens = density(x[x[,2] != 1, 2])
    xlim = c(0, max(1, dens$y * mean(x[,2] != 1)))
    plot(0, type ='n', xlim = xlim, ylim = c(0,1), axes = FALSE)
    lines(c(0, mean(x[,2] == 1)), c(1,1))
    lines(dens$y * mean(x[,2] != 1), dens$x)
    box()
    axis(1)
    par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
    }

# Rough sketch of the joint. The 3d plot just shows the marginals (the red may need
# to be reversed), so for it to be the joint, we'd need to scale it down more.
plotmax.3d = function(x){
    require(rgl)
    dens1 = density(x[x[,1] != 1, 1])
    dens2 = density(x[x[,2] != 1, 2])
    zlim = c(0, max(1, dens1$y * mean(x[,1] != 1), dens2$y * mean(x[,2] != 1)))

    plot3d(cbind(1,1,1), type='n', xlim = c(0,1), ylim = c(0,1), zlim = zlim,
        xlab = "x", ylab = "y", zlab = "z")

    lines3d(c(1.05,1.05), c(1,1), c(0, mean(x[,1] == 1)), col = 'red', lwd = 2)
    lines3d(c(1,1), c(1.05,1.05), c(0, mean(x[,2] == 1)), col = 'blue', lwd = 2)

    lines3d(rep(1.05, length(dens1$x)), dens1$x,
        dens1$y * mean(x[,1] != 1), col = 'red', lwd = 2)

    lines3d(dens2$x, rep(1.05, length(dens2$x)),
        dens2$y * mean(x[,2] != 1), col = 'blue', lwd = 2)
    }


n = 100000
p = 0.7
q = 0.3
a = 2
b = 6
c = 7
d = 4

x = rmax(n, p, q, a, b, c, d)

plot(x, pch = 20)



mean(x[,1] == 1)
mean(x[,2] == 1)
mean(apply(x, 1, min) == 1)

# Probability that first column equals 1:   (1 - p) + q*p
# Probability that second column equals 1:  p
# Probability they both equal 1:            p*q


mean(x[,1] != 1)
mean(x[,2] != 1)

# Prob that first is beta:  p*(1-q)
# Prob that second is beta: 1-p
