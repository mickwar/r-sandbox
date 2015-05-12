### Some color and plotting functions for nicer density plots.

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

int2rgb = function(x){
    hex = as.character(as.hexmode(x))
    hex = paste0("#", paste0(rep("0", 6-nchar(hex)), collapse=""), hex)
    col2rgb(hex)
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
