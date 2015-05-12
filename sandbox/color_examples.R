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

# example 2 (making a density completely shaded with color and highlighting the middle)
x = rnorm(1000)
dens = density(x)
plot(dens)
polygon(dens$x, dens$y, col = "dodgerblue")
color.den(dens, -2, 2, "mediumseagreen")
lines(dens)

### int2rgb()
# Convert an integer/hexdecimal to an RGB. Returns the RGB values.
# 
# Params: x - an integer, between 0 and 16777215 = 256^3 - 1,
#             or between 0x000000 and 0xFFFFFF
int2rgb = function(x){
    hex = as.character(as.hexmode(x))
    hex = paste0("#", paste0(rep("0", 6-nchar(hex)), collapse=""), hex)
    col2rgb(hex)
    }

# example 1
int2rgb(5)
int2rgb(0x000005)

# example 2
int2rgb(1193046)
int2rgb(0x123456)

### col.mult()
# Multiplies two colors together.
# If A and B are vectors in [0,1]^3 (i.e. {R, G, B}), then to multiply the
# colors A and B, do elementwise multiplication.
#
# Params: col1 - a color, either in name or integer/hexidecimal
#         col2 - a color, either in name or integer/hexidecimal
#
# requires int2rgb() if using integer/hexidecimal form
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

# example 1
col.mult(0x123456, "black") # black times anything is always black

# example 2
col.mult(0x123456, "white") # multiplying by white doesn't change anything

# example 3
col.mult(0x123456, 0x888888) # multiplying by gray darkens the color equally in RGB

# example 4
x = rnorm(1000)
dens = density(x)
plot(dens)
polygon(dens$x, dens$y, col = "dodgerblue")

color.den(dens, -2, 2, "mediumseagreen")
lines(dens)

color.den(dens, -2, 2, col.mult("dodgerblue", "mediumseagreen"))
lines(dens)

color.den(dens, -2, 2, col.mult("dodgerblue", "gray50"))
lines(dens)

# By multiplying colors together for highlighting purposes, the plot looks
# less distracting by having appropriately related colors

# example 5
x = rnorm(1000)
dens = density(x)
plot(dens)
polygon(dens$x, dens$y, col = col.mult("dodgerblue", "gray50"))
color.den(dens, -2, 2, "dodgerblue")
lines(dens)

# It seems to work better when the darker region is the area of interest

### colgray()
# Convert a given RGB color to gray-scale
#
# Params: rgb    - a vector of length 3 where each element is between 0 and 1
#         method - either "yprime", "average", or "lightness"
#                  "yprime" is the method to account for human vision
col.gray = function(rgb, method = "yprime"){
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


# example 1
# library(jpeg)
# path = "~/files/R/data/"
# dat = readJPEG(paste0(path,"image.jpg"))
# 
# h = nrow(dat) # height in pixels
# w = ncol(dat) # width in pixels
