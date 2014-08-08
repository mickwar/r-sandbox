library(jpeg)
path = "~/files/R/data/"
dat = readJPEG(paste0(path,"image.jpg"))

h = nrow(dat) # height in pixels
w = ncol(dat) # width in pixels

col2gray = function(rgb, method="yprime"){
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
