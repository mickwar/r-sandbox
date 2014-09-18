### needs to be put into bayes_functions
### compare this hpd.mult() with the one in bayes_functions.R
### and see the accuracy and speed of each
color.den = function(dens, from, to, col1 = 1, col2 = NULL){
    if (is.null(col2))
        col2 = col1
    index = which(dens$x > from & dens$x < to)
    polygon(c(dens$x[index][1], dens$x[index],
        dens$x[index][length(index)]), c(0, dens$y[index], 0),
        col=col1, border=col2)
    }
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
# x is 
hpd.mult = function(x, dens, prob = 0.95, tol, interactive = TRUE){
    max.k = max(dens$y)
    min.k = min(dens$y)
#   k = (max.k - min.k)/2
    k = runif(1, min.k, max.k)
    if (missing(tol))
        tol = max.k / 10000
    count = 0
    if (interactive){
        plot(dens)
        polygon(dens$x, dens$y, col='gray80')
        }
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
    while (max.k - min.k > tol){
        count = count + 1
        c.prob = 0
        int = zeros(dens$y, k)
        if (is.null(int)){
            int = c(1, length(dens$x))
            }
        # sum the area in the intervals
        for (j in 1:(length(int)/2))
            c.prob = c.prob + mean(x >= dens$x[
                 int[2*j-1]] & x <= dens$x[int[2*j]])
        if (interactive){
            abline(h=k)
            abline(h=c(max.k, min.k), col='blue')
            cat("Probability:",c.prob)
            readline()
            }
        # not a large enough region, lower k
        if (c.prob < prob)
            max.k = k
        # too much region, raise k
        if (c.prob > prob)
            min.k = k
        # right-e-o!
        if (c.prob == prob)
            max.k = min.k
        # pick new height to test at
        k = runif(1, min.k, max.k)
        }
    if (interactive)
        cat("Iterations:",count,"\n")
    return (dens$x[int])
    }


### bimodal
n = 10000
x = ifelse(runif(n) < 0.5, rnorm(n, 2.5, 1), rnorm(n, 6, 0.5))
dens = density(x, n = 20000)
plot(density(x))

#hpd.1 = hpd.mult(x, dens, visual = FALSE, print = FALSE)
hpd.1 = hpd.mult(x, dens)

plot(density(x), main="Mixture Normal", xlab="x", ylab="f(x)")
polygon(dens$x, dens$y, col='gray')
col1 = col.mult('gray', 'dodgerblue')
color.den(dens, hpd.1[1], hpd.1[2], col1 = 'dodgerblue')
color.den(dens, hpd.1[3], hpd.1[4], col1)
points(density(x), type='l')
abline(h=0)


#abline(v=hpd.1, lty=2, col='gray30', lwd=2)


#### more bimodal
#x = ifelse(runif(n) < 0.5, rnorm(n, -1, 1), rnorm(n, 6.5, 0.5))
#dens = density(x)
#plot(density(x))
#
#hpd.1 = hpd.mult(x, dens)
#abline(v=hpd.1, lty=2, col='gray')
