##### Bayes functions
##### This script contains various functions that have been
##### useful in Bayesian analysis, particularly for posterior
##### distributions from a Metropolis-Hastings algorithm

##### M-H FUNCTIONS #####
# autotune()
# Used to adjust candidate sigmas for normal candidate densities
# Requires the calculation of acceptance rates within some given
# window of mcmc iterations. For example, every 500 draws compute
# the acceptance rate (0 <= A <= 1) for each parameter using the
# last 500 draws. Multiply each candidate sigma by autotune().
# 
# arguments:
# accept - the acceptance rate within a given window
# target - the acceptance rate goal
# k - the maximum the candidate sigma can change by, 1/k is 
#     minimum it can change by. For instance, if the current
#     sigma is 2, and in some window of mcmc iterations every
#     draw was accepted, the new candidate sigma would now
#     be 2*k, which should serve to reduce the acceptance rate.
#     On the other hand, if no draws are accepted, the new candidate
#     sigma would then be 2/k. I recommend k = window/50
#     Larger values of k will change sigma by a larger amount,
#     and vice versa for smaller values (k -> 1, autotune() -> 1,
#     everywhere)
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
##### END M-H FUNCTIONS #####

##### DENSITY FUNCTIONS #####
# a collection of density()-related functions, designed to
# require a vector of random draws (e.g. rnorm(10000)) by
# taking in a density object as the argument

# note: argument "x" refers to vector of random draws
#       argument "dens" refers to density(x)

# calc.mode()
# estimate the mode, method == "loess" gives more accurate
# results, but method == "density" takes less computation time
calc.mode = function(dens, method="loess"){
    dx = dens$x
    dy = dens$y
    if (method == "loess"){
        l = loess(dy ~ dx)
        return (l$x[which.max(l$y)])
        }
    if (method == "density")
        return (dx[which.max(dy)])
    }

# bound()
# get the closest values in dens
# this function is more useful in just getting
# a vertical line from y=0 to the density (for
# say an interval without filling with polygon)
bound = function(x, dens, return.x=TRUE){
    # returns the x-value in dens that is closest
    # to the given x
    if (return.x)
        return(dens$x[which.min(abs(dens$x-x))])

    # returns the y-value in dens at the closest x
    return(dens$y[which.min(abs(dens$x-x))])
    }

# color.den()
# (by Arthur Lui, editted by Mickey)
# Colors area under a density within an interval
# dens has to be a density object
# (note: there seems to be an issue with values to close
#        together, not sure exactly, but had errors)
color.den = function(dens, from, to, col1 = 1, col2 = NULL){
    if (is.null(col2))
        col2 = col1
    index = which(dens$x > from & dens$x < to)
    polygon(c(dens$x[index][1], dens$x[index],
        dens$x[index][length(index)]), c(0, dens$y[index], 0),
        col=col1, border=col2)
    }

# col.mult()
# Multiplies to colors together. If A and B are vector
# in [0,1]^3 (i.e. {R, G, B}), then to multiply the
# colors A and B, do elementwise multiplication.
# col1 and col2 are expected to be given R color names
# or integers of the form 0xRRGGBB (they can really be
# any integer, but I think having the colors in their
# hexadecimal form is easiest to see)
# requires int2rgb()
col.mult = function(col1 = 0x000000, col2 = "black"){
    if (is.character(col1))
        val1 = t(col2rgb(col1) / 255)
    if (is.numeric(col1))
        val1 = t(int2rgb(col1) / 255)
    if (is.character(col2))
        val2 = t(col2rgb(col2) / 255)
    if (is.numeric(col2))
        val2 = t(int2rgb(col2) / 255)
    rgb(val1 * val2)
    }

# int2rgb()
# convert an integer between 0 and 16777215 = 256^3 - 1,
# or between 0 and 0xFFFFFF
# this function is depended upon by col.mult
int2rgb = function(x){
    hex = as.character(as.hexmode(x))
    hex = paste0("#", paste0(rep("0", 6-nchar(hex)), collapse=""), hex)
    col2rgb(hex)
    }

# colgray()
# convert a given RGB color to gray-scale
# "yprime" is the method to account for human vision,
# and is supposed to be best
# rgb is a vector in [0, 1]^3
col.gray = function(rgb, method="yprime"){
    if (method == "yprime"){
        weight = c(0.2126, 0.7152, 0.0722) * rgb
        return (rep(sum(weight), 3))
        }
    if (method == "average"){
        return (rep(mean(rgb), 3))
        }
    if (method == "lightness"){
        return (rep(mean(range(rgb)), 3))
        }
    }
# hpd.plot()
# plot a univariate density with its hpd shaded
# dens     - density object to be plotted
# hpd      - vector containing end points of the hpd region, must
#            satisfy length(hpd) %% 2 == 0
# col1     - the color of the non-shaded portion of the plot
# col2     - the color of the shaded portion, defaults to gray50
# multiply - logical, if true multiply col1 with col2 and set new
#            color to col2, otherwise don't
# border   - color of the density border
# ...      - arguments to pass to plot()
# requires color.den() and col.mult(), which requires int2rgb()
# note: the use of the border argument could be improved (i don't think
#       null would disable the border color if desired)
hpd.plot = function(dens, hpd, col1 = "dodgerblue", col2 = NULL,
    multiply = TRUE, border = "black", ...){
    if (is.null(col2))
        col2 = "gray50"
    if (multiply)
        col2 = col.mult(col1, col2)
    plot(dens, type='n', ...)
    polygon(dens, col=col1, border = NA)
    for (i in 1:(length(hpd)/2))
        color.den(dens, hpd[2*i-1], hpd[2*i], col2)
    color.den(dens, -Inf, min(hpd), col1)
    color.den(dens, max(hpd), Inf, col1)
#   lines(dens, col = border)
    }

# plot.post() original by Arthur Lui (github.com/luiarthur)
# adapted by mickey to work with hpd.plot()
# Plots a shaded posterior density with hpd.plot() and adds a trace
# plot in the corner.
# x          - 
# dens       - density object: density(x)
# hpd        - vector containing end points of the hpd region, must
#              satisfy length(hpd) %% 2 == 0
# right      - logical. should the trace plot be in the top right corner?
#              if false, place in top left corner.
# subfig.mar - margins for the subfigure
# ...        - arguments to pass to hpd.plot() (and plot())
#        hpd.plot() takes arguments hpd, col1, col2, multiply, and border
#        (see hpd.plot() for more information), these must be named
plot.post = function(x, dens, hpd, right = TRUE,
    subfig.mar = c(0.1, 0.1, 1.0, 0.1), ...) {
    if (missing(dens))
        dens = density(x)

    rng = range(dens$y)
    diff = diff(rng)

    opts2 <<- par(no.readonly = TRUE)

    # first plot the hpd
    hpd.plot(dens, hpd, ylim = c(rng[1], rng[2] + diff * 0.5), ...)

    rng.x = range(dens$x)
    x.diff = diff(rng.x)
  
    opts = par(no.readonly = TRUE)
    if (right){
        left = rng.x[1] + x.diff*2/3
        right = rng.x[2]
        right = opts$usr[2]
    } else {
        left = rng.x[1]
        left = opts$usr[1]
        right = rng.x[2] - x.diff*2/3
        }
    par(fig = c(grconvertX(c(left, right), from="user", to="ndc"),
        grconvertY(c(rng[2], rng[2] + diff * 0.5), from="user", to="ndc")),
        mar = subfig.mar, new = TRUE)
    #plot(density(x),col="blue",cex.main=.5,lwd=3)
    plot(x, type="l", col="gray20", cex.main=.5, axes=F, main="Trace Plot")
    axis(1, cex.axis=.5)
    axis(2, cex.axis=.5)
    par(opts)
    }

# need to remove eventually, testing out the plot.post() function
# y = rnorm(10000)
# dens = density(y)
# hpd = hpd.uni(y)
# plot.post(y, dens, hpd, right = F, subfig.mar = c(0,0,0,.1), main = "Posterior for y", xlab = "")

# hpd.uni()
# functions to compute highest posterior density regions
# for unimodal distributions
hpd.uni = function(x, prob = 0.95, precision = 1000){
    range = seq(0, 1-prob, length=precision)
    range = cbind(range, range+prob)
    best = range[which.min(apply(range, 1, function(y)
        diff(quantile(x, y)))),]
    return (quantile(x, best))
    }

# hpd.mult()
# for multi-modal distributions (hpd could be a set)
hpd.mult = function(x, dens, prob = 0.95, klength = 5000){
    c.prob = 1
    temp.prob = 0
    k = seq(min(dens$y), max(dens$y), length=klength)
    # i = 2 to prevent certain problems, test.x4 had an issue
    # with the probability 0 region in the middle, (doesn't always
    # occur) perhaps fix by doing f(x) > k, instead of f(x) >= k?
    i = 2
    zeros = function(y, k, return.max.min = FALSE){
        # y is expected to be density(x)$y
        out = NULL
        int = NULL
        for (i in 2:(length(y)-1)){
            # condition to check when the height crosses k
            if ((y[i] > k && y[i-1] < k) || (y[i] < k && y[i-1] > k)){
                # get the x closer to k
                out = c(out, ifelse(which.min(abs(y[c(i,i-1)]-k))==1,
                    i, i-1))
                # -1 if lower interval, +1 if upper
                int = c(int, -sign(y[i] - y[i-1]))
                }
            # check if the height exactly equals k
            if (y[i] == k){
                out = c(out, i)
                # 0 if a maximum or minimum, -1 if lower, +1 if upper
                # y[i] can only be max or min if y[i] = k, so don't
                # check this condition for when height crosses k
                int = c(int, -sign(sign(y[i]-y[i-1]) +
                    sign(y[i+1]-y[i])))
                }
            }
        # ensure that first value is lower end and last is upper end
        if (is.null(int))
            return (NULL)
        if (int[1] == 1){
            int = c(-1, int)
            out = c(1, out)
            }
        if (int[length(int)] == -1){
            int = c(int, 1)
            out = c(out, length(y))
            }
        if (return.max.min)
            return (out)
        return (out[as.logical(int)])
        }
    # repeat until area under curve is <= specified prob
    # (note 14 jun: perhaps do some kind of iterative
    # convergence to reduce the number of iterations;
    # start in the middle i = floor(klength/2), and if
    # temp.prob is too low, set i=floor((i + 0)/2), else
    # set i = floor((i + klength)/2), repeat until value
    # is sufficiently close to prob. need to keep track of
    # previous "lower" and "upper" bounds
    while (c.prob > prob){
        temp.prob = 0
        int = zeros(dens$y, k[i])
        if (!is.null(int)){
            if (length(int) > 0){
                # sum the area in the intervals
                for (j in 1:(length(int)/2))
                    temp.prob = temp.prob + mean(x >= dens$x[
                         int[2*j-1]] & x <= dens$x[int[2*j]])
                # update (i think this always occurs)
                if (c.prob > temp.prob)
                    c.prob = temp.prob
                }
            }
        i = i + 1
        }
    return (dens$x[int])
    }
##### END DENSITY FUNCTIONS #####

##### OTHER FUNCTIONS #####
### bayes chi^2 goodness of fit given by V.E.Johnson (2004)
# so far assumes y is univariate, i.e. data = cbind(y, x)
# params are mcmc draws from the joint posterior
# cdf(data, params) is a function for the cdf of the likelihood
# K is the number of bins to use for the chi2 gof test
# alternatively, a can be given which contains the endpoints
# for a_0 to a_K
# returns the pvals
# (the null is that the model fits)
bayes.gof = function(data, params, cdf, K, a, every = 100){
    y = as.matrix(data)
    n = nrow(y)
    M = nrow(params)
    nparams = ncol(params)
    if (missing(K) && missing(a))
        K = round(n ^ 0.4)
    if (missing(a))
        a = seq(0, 1, length = K + 1)
    K = length(a) - 1 # for when only a is given
    p = diff(a) # need to check this for other values of a

    m = matrix(0, M, K)
    for (j in 1:M){
        if (floor(j/every) == j/every)
            cat("\rIteration:",j,"/",M)
        z = cdf(data, params[j,])
        temp.m = rep(K, n)
        for (k in 2:K)
            temp.m = temp.m - ifelse(z <= a[k], 1, 0)
        m[j,as.numeric(names(table(temp.m)))] = as.numeric(table(temp.m))
        if (j == M)
            cat("\n")
        }

    chi = double(M)
    for (l in 1:M)
        chi[l] = sum((m[l,] - n*p)^2 / (n*p))

    pvals = pchisq(chi, K - 1, lower.tail = FALSE)
    }
##### END OTHER FUNCTIONS #####
