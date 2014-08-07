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
#     minimum it can change by. For instance, is the current
#     sigma is 2, and in some window of mcmc iterations every
#     draw was accepted, the new candidate sigma would now
#     be 2*k, which should serve to reduce the acceptance rate.
#     On the other hand, if no draws are accepted, the new candidate
#     sigma would then be 2/k. I recommend k = window/50
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
# den has to be a density object
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
# trailing()
# truncate the numeric x "digits" places after the decimal
# useful for printing, especially number not much greater than 1
# or very, very close to 0.
trailing = function(x, digits = 4)
    formatC(x, digits=digits, format="f")

# nice.time()
# returns a string in nice time format given the number of seconds
# I have used this to estimate the remaining time in long simulations,
# in conjunction with as.numeric(Sys.time())
nice.time = function(seconds){
    # floor() or round() would work as well
    seconds = ceiling(seconds)
    days = seconds %/% (60*60*24)
    seconds = seconds %% (60*60*24)
    hours = seconds %/% (60*60)
    seconds = seconds %% (60*60)
    minutes = seconds %/% (60)
    seconds = seconds %% (60)
    return (paste0(days,"d ",hours,"h ",minutes,"m ",seconds, "s"))
    }

# mcmc_time()
# do - logical, should the function be run
# iter - value of 0 initializes the process
#   any other value increments the process
#   the function assumes iter is the same value
#   as the index of the parameter matrix
#   being modified. for instance, the m-h
#   loop often begins on iteration 2, and the
#   function assumes this
# every - how often should the results be written
#   to external text files. it should be the case
#   that (nmcmc + nburn) %% every == 0. for
#   computationally intensive loops, it's a good
#   idea to output the results periodically. if
#   the loop takes only a few minutes, then set
#   every = nburn + nmcmc
# params, accept, sigs - matrices and vectors 
#   containing parameter draws, binary indicator
#   of acceptances, and candidate sigmas 
# nburn - number of burnin draws
# nmcmc - number of mcmc iterations (does not
#   include burnin, i.e. total iterations =
#   nburn + nmcmc
# This function uses variables in global scope and
# so are reserved. These are:
#   mcmc_begin_time
#   mcmc_curr_time
# Requires function nice.time()
mcmc_time = function(do = TRUE, iter, every = 100, params, accept,
    sigs, nburn, nmcmc, dir = ".", prefix = "mcmc_"){
    if (!do)
        return (0)
    if (iter == 0){
        mcmc_begin_time <<- as.numeric(Sys.time())
    } else {
        mcmc_curr_time <<- as.numeric(Sys.time()) - mcmc_begin_time
        spi = mcmc_curr_time / (iter - 1)
        write.table(paste0("Iteration: ", iter, "\n",
            "  Burn-in: ", nburn, "\n",
            "    NMCMC: ", nmcmc, "\n\n",
            "  Elapsed: ", nice.time(mcmc_curr_time), "\n",
            "Remaining: ", nice.time(spi*(nmcmc+nburn-iter)), "\n\n",
            " Complete: ", trailing(100*iter/(nburn+nmcmc)),"%"),
            paste0(dir, "/", prefix, "status.txt"),
            col.names=F, row.names=F, quote=F)
        if (floor(iter/every) == iter/every){
            write.table(signif(params[(1+every*(ceiling(iter/every)
                -1)):iter,]), paste0(dir, "/", prefix, "params.txt"),
                row.names=F, col.names=F, quote=F, append=T)
            write.table(signif(accept[(1+every*(ceiling(iter/every)
                -1)):iter,]), paste0(dir, "/", prefix, "accept.txt"),
                row.names=F, col.names=F, quote=F, append=T)
            # print the latest parameter values and candidate sigmas
            # until burn in is complete, so the values may be used
            # as starting point for future runs
            if (iter <= nburn + 1){
                write.table(params[iter,],
                    paste0(dir, "/", prefix, "init_params.txt"),
                    row.names=F, col.names=F, quote=F)
                write.table(sigs,
                    paste0(dir, "/", prefix, "cand_sigmas.txt"),
                    row.names=F, col.names=F, quote=F)
                }
            }
        }
    }
##### END OTHER FUNCTIONS #####
