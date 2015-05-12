### Some color and plotting functions for nicer density plots.

### color.den()
# (by Arthur Lui, editted by Mickey)
# Colors a specified area under a density within an interval
# This should be used in conjunction with plot(dens, ...)
#
# Params: dens   - a density object
#         from   - numeric, the left endpoint of the area
#         to     - numeric, the right endpoint of the area
#         fill   - a color, the color that will fill in the area under the curve
#         border - a color, the outline of the area receives this color
#
# (note: there seems to be an issue with values to close
#        together, not sure exactly, but had errors)
color.den = function(dens, from, to, fill = 1, border = NULL){
    if (is.null(border))
        border = fill
    index = which(dens$x > from & dens$x < to)
    x0 = dens$x[index][1]
    x1 = dens$x[index][length(index)]
    x = dens$x[index]
    y = dens$y[index]
    
    # creating the shaded region
    polygon(c(x0, x, x1), c(0, y, 0), col = fill, border = border)
    }

# example 1
x = rnorm(1000)
dens = density(x)
plot(dens)
color.den(dens, -2, 2, "green3") # green3 is equivalent to color 3
lines(dens) # plot the density again since color.den() will draw over the original curve

# example 2
x = rnorm(1000)
dens = density(x)
plot(dens)
polygon(dens$x, dens$y, col = "dodgerblue")
color.den(dens, -2, 2, "mediumseagreen")
lines(dens)


### col.mult()
# Multiplies two colors together.
# If A and B are vectors in [0,1]^3 (i.e. {R, G, B}), then to multiply the
# colors A and B, do elementwise multiplication.
#
# Params: col1 and col2 are expected to be given R color names
# or integers of the form 0xRRGGBB (they can really be any integer, but I
# think having the colors in their hexadecimal form is easiest to see)
#
# requires int2rgb() if not using color names (e.g. "black")
#
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



init = array(0, dim=dim(dat))

library(jpeg)
path = "~/files/R/data/"
dat = readJPEG(paste0(path,"image.jpg"))

h = nrow(dat) # height in pixels
w = ncol(dat) # width in pixels
