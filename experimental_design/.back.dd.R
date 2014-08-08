# Find a near-maximin, near-orthogonal subset of design points.
# Motivated by a set of computer experiments containing a very large
# number of runs. We wish to find a subset of the runs which have
# properties similar to Latin hypercubes having maximin and
# orthogonality properties since use of the full dataset will be very
# costly in obtaining estimates for a Gaussian process model where
# matrix inversion is required.

# todo: run diagnostics to see how well the chosen subset
#       meets desired properties (maximin, orthogonal factors, uniform
#       spread over marginals, latin hypercube, no repitition of
#       factor levels, and others?)
msub = function(x, m){
    # x: the design matrix, number of points by number of dimensions
    # m: <= nrow(x), the size of the subset to find
    # note: there seems to be a tendency to choose around the
    #       edges before filling in the middle, but a high enough
    #       m this shouldn't be too much of an issue
    # perhaps other distance metrics may be used. i tried minkowski
    # distance to compute the 'norms' variable below and all choices
    # of p gave the same result, but maybe not so in the j for loop
    euclid = function(x)
        sqrt(sum(x^2))
    n = nrow(x)
    d = ncol(x)
    xmin = apply(x, 2, min)
    xrange = apply(x, 2, function(x) diff(range(x)))
    # scale x to [-0.5, 0.5]^d (really just needs to be centered
    # around 0)
    for (i in 1:d)
        x[,i] = ((x[,i] - xmin[i]) / xrange[i]) - 0.5

    # a default value for m
    if (missing(m))
        m = floor(n^(1/d))
    out = double(m)

    # compute distance from the origin
    norms = apply(x, 1, euclid)
    # take the first point in the subset to be that which is furthest
    # from the origin. the point will be on the edge and more likely
    # to be in a "corner" than, say, on the wider edge of an ellipse
    out[1] = which.max(norms)

    # initialize the dists matrix. these values are the distances
    # from each point in the design and the ones currently selected
    # to be in the subset. the column of Inf's is included so the
    # for loop will always consider dists to be a matrix on not
    # automatically change it to a vector (so apply always works)
    dists = matrix(0, n, m-1)
    dists = cbind(Inf, dists)

    corrs = matrix(0, n, m-1)

    for (j in 2:m){
        dists[,j] = apply(x - matrix(x[out[j-1],], n, d, byrow=TRUE),
            1, euclid)

        # downweight the distances if in a particular dimension the
        # locations are close to previously selected points in the
        # subset. want to encourage latin hypercube similarities.

        # also downweight higher correlated additions

        # in the downweight, i may want to add a small term so
        # the distances don't go exactly to zero
        for (k in 1:d){
            # computing compmax may be unnecessary since all points
            # are divided by that same scaling
            compmax = max(abs(x[out[j-1], k] - x[,k]))
            dists[,j] = dists[,j] * (0.5 + abs(x[out[j-1], k]-x[,k]) /
                compmax)
            }
        if (j > 2){
            for (l in 1:n)
                dists[l,j] = prod(dists[l,j], sqrt(prod(1.5 - abs(
                    cor(x[c(l, out),])) + diag(d))))
            }
        out[j] = which.max(apply(dists[,1:j], 1, min))
        }
    # returns the index of points in the subset
    if (length(unique(out)) < length(out))
        print("Duplicates selected in the subset. Removing.")
    return (unique(out))
    }







# example
n = 50
d = 2
m = 10
x = matrix(runif(n*d), n, d)
# x = matrix(rbeta(n*d, 1, 0.25), n, d)
# x = matrix(mvrnorm(n,c(0.5,0.5),matrix(c(1,0.8,0.8,1),2,2)),n,d)

# takes a bit over a minute (n = 10000, d = 10, m = 500)
system.time(out <- msub(x, m))

plot3d(x)
points3d(x[out,], col='red', size=5)

pairs(x, pch=20)
for (i in 2:m){
    readline()
    pairs(x[out[1:i],], col='red')
    }

# for d = 2
plot(x, pch=20)
points(matrix(x[out,], ncol=2),col='red',cex=1.5, lwd=2)

# see the choice of subset points one at a time
for (i in 1:m){
    points(matrix(x[out[i],], ncol=2),col='red',cex=1.5, lwd=2)
    readline()
    }

# diagnostics should be run on x[out,], it being the subset
# out (from out <- msub(x, m)) is the index of points

###
x = cbind(rep(seq(0, 1, length=11), each=11),
    rep(seq(0, 1, length=11), times=11))
n = nrow(x)
d = ncol(x)

out = msub(x, 10)

plot(x, pch=20)
for (i in 1:length(out)){
    readline()
    points(x[out[i],1],x[out[i],2], col='red', cex=1.5)
    }

cor(x[out,])
# notice the correlation is at 0.47878

dist(x[out,])

y = matrix(c(0,0.2,0,0.5), 2,2)
plot(y, pch=20)
cor(y)
