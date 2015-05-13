# Example usage of the functions in color_functions.R
source("../useful/color_functions.R")

### color.den examples
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


### int2rgb examples
# example 1
int2rgb(5)
int2rgb(0x000005)

# example 2
int2rgb(1193046)
int2rgb(0x123456)


### col.mult examples
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


### col.gray examples
# example 1
# library(jpeg)
# path = "~/files/R/data/"
# dat = readJPEG(paste0(path,"image.jpg"))
# 
# h = nrow(dat) # height in pixels
# w = ncol(dat) # width in pixels
