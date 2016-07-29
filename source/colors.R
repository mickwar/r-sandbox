### Some color and plotting functions for nicer density plots.

# get current list of objects in the console
TEMP_LIST_O = ls()

########## BEGIN



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

### fade()
# Make a color transparent, simple call to rgb()
#
# Params: color - either a string of a valid R color name or a vector of
#                 length 3 of RGB values in [0, 1]^3
#         alpha - the degree of transparency, 0 is invisible and 1 is no
#                 transparency
fade = function(color, alpha = 0.5){
    if (length(color) == 1){
        x = col2rgb(color) / 255
    } else {
        x = color
        }
    rgb(x[1], x[2], x[3], alpha)
    }

### make.phantom()
# Takes a character vector and can hide some components in that vector
# when making a title. Also allows different colors to be in the text.
# Adapted from Matt Heiner
#
# Params: text    - a vector of characters, each part that is desired to
#                   be either invisible or a different color should be
#                   in a separate element
#         display - numeric vector containing the index of which elements
#                   in text should be displayed, defaults to displaying all
#         colors  - either a single color or a vector of colors having length
#                   equal to the length of text so the colors match with
#                   the elements of text, defaults to black
#         sep     - a character string that separates the components in
#                   text, defaults to "" which is no separation
#         ...     - additional arguments used in title()
make.phantom = function(text, display, colors, sep = "", ...){
    # text: a character vector for 
    n = length(text)
    if (missing(display))
        display = 1:n
    if (missing(colors))
        colors = rep("black", n)
    if (length(colors))
        colors = rep(colors, n)
    for (i in display){
        if (i == 1)
            title(main = bquote( .(text[i]) * phantom(.(sep)) * phantom(.(paste0(text[-i], collapse = sep)))),
                col.main = colors[i], ...)
        if (i == n)
            title(main = bquote( phantom(.(paste0(text[-i], collapse = sep))) * phantom(.(sep)) * .(text[i])),
                col.main = colors[i], ...)
        if (i > 1 && i < n)
            title(main = bquote( phantom(.(paste0(text[1:(i-1)], collapse = sep))) * phantom(.(sep)) *
                .(text[i]) * phantom(.(sep)) * phantom(.(paste0(text[(i+1):n], collapse = sep)))),
                col.main = colors[i], ...)
        }
    }



########## END

# compare TEMP_LIST with ls() which now has additional functions
TEMP_LIST_N = ls()
TEMP_LIST_N = TEMP_LIST_N[-which(TEMP_LIST_N == "TEMP_LIST_O")] # remove TEMP_LIST_O
TEMP_LIST_N = TEMP_LIST_N[!(TEMP_LIST_N %in% TEMP_LIST_O)] # get only added functions

# print the sourced functions
if (length(TEMP_LIST_N) > 0){
    cat("Sourced functions:\n    ")
    cat(TEMP_LIST_N, sep="\n    ")
    }

rm(TEMP_LIST_O, TEMP_LIST_N) # remove temporary variables
