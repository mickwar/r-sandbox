### todo
# default labeling
# less hardcoding some of the values for spacing (in
#   case i want to make adjustments)
# add options for vline and hline

textplot = function(x, y = NULL, pch="*",
    max.width = 150, max.height = 40,
    xlim = NULL, ylim = NULL,
    main = NULL, main.justify = NULL,
    xlab = NULL, ylab = NULL,
    xaxt = NULL, xaxt.nchar.max = 8,
    yaxt = NULL, yaxt.nchar.max = 8,
    xtick = TRUE, ytick = TRUE,
    box = TRUE){

    ### create functions
    # for repeating characters using "cat"
    fill = function(char, n)
        paste0(rep(char, n), collapse="")
    # to make cat default to sep="", like paste0
    cat0 = function(...)
        cat(..., sep="")
    sig = function(x)
        signif(x, 2)
    # to transform x from range(x) to [a, b]
    rescale = function(x, a, b)
        floor((x - min(x))/(max(x)-min(x)) * (b - a) + a)

    ### set up
    # process arguments
    if (any(is(x) == "density")){
        y = x$y
        x = x$x
        if (is.null(ylim))
            ylim = c(0, max(y))
        }
    if (is.null(y)){
        y = x
        x = seq(1, length(x), by = 1)
        }
    if (length(x) != length(y))
        stop("Lengths of x and y are unequal.")
    if (is.null(xlim))
        xlim = range(x)
    if (is.null(ylim))
        ylim = range(y)
    if (diff(range(xlim)) == 0)
        xlim = c(xlim[1] - 1, xlim[2] + 1)
    if (diff(range(ylim)) == 0)
        ylim = c(ylim[1] - 1, ylim[2] + 1)
    if (is.null(main))
        main = ""
    # main.justify is one of "c", "l", "r"
    if (is.null(main.justify) || length(main.justify) == 1)
        main.justify = rep("c", length(main))
    if (length(main.justify) != length(main))
        stop("Number of main titles given unequal to number of justifications.")
    if (is.null(xlab))
        xlab = "x"
    if (is.null(ylab))
        ylab = "y"
    # if xaxt and yaxt are provided, the given values will be
    # printed at the ticks, regardless of xlim and ylim
    if (xtick){
        if (is.null(xaxt)){
            # these values are taken from looking at plot() very closely
            f1 = log(1 + 5/9, 10)
            f2 = log(3 + 8/9, 10)
            f3 = log(7 + 7/9, 10)
            # mag - the digit place that the values should be changed at,
            # i.e. change at 10^mag
            mag = log(diff(range(xlim)), 10)
            prop.mag = mag - floor(mag)
            # check to see what spacing to use (1, 2, or 5)
            # spacing by 2
            if (prop.mag >= 0 && prop.mag < f1 + 1 - f3)
                width = 2*10^(floor(mag) - 1)
            # spacing by 5
            if (prop.mag >= f1 + 1 - f3 && prop.mag < f2 + 1 - f3)
                width = 5*10^(floor(mag) - 1)
            # spacing by 1
            if (prop.mag >= f2 + 1 - f3 && prop.mag < 1)
                width = 10^floor(mag)
            
            # creating spanning xaxt
            xaxt = seq(width*floor(xlim[1]/width -1), width*ceiling(xlim[2]/width+1), by = width)
    #       if (xlim[1] >= 0)
    #           xaxt = seq(0, xlim[2], by = width)
    #       if (xlim[2] <= 0)
    #           xaxt = seq(0, xlim[1], by = -width)
    #       if (xlim[1] < 0 && xlim[2] > 0){
    #           temp.1 = seq(0, xlim[1], by = -width)
    #           temp.2 = seq(0, xlim[2], by = width)
    #           xaxt = unique(sort(c(temp.1, temp.2)))
    #           }
            xaxt = xaxt[which(xaxt >= xlim[1] & xaxt <= xlim[2])]
            }
        xaxt = as.character(xaxt)
        xaxt = substr(xaxt, 1, xaxt.nchar.max)
        xaxt.len = length(xaxt)
        }
    if (ytick){
        if (is.null(yaxt)){
            # these values are taken from looking at plot() very closely
            f1 = log(1 + 5/9, 10)
            f2 = log(3 + 8/9, 10)
            f3 = log(7 + 7/9, 10)
            # mag - the digit place that the values should be changed at,
            # i.e. change at 10^mag
            mag = log(diff(range(ylim)), 10)
            prop.mag = mag - floor(mag)
            # check to see what spacing to use (1, 2, or 5)
            # spacing by 2
            if (prop.mag >= 0 && prop.mag < f1 + 1 - f3)
                width = 2*10^(floor(mag) - 1)
            # spacing by 5
            if (prop.mag >= f1 + 1 - f3 && prop.mag < f2 + 1 - f3)
                width = 5*10^(floor(mag) - 1)
            # spacing by 1
            if (prop.mag >= f2 + 1 - f3 && prop.mag < 1)
                width = 10^floor(mag)
            
            # creating spanning yaxt
            yaxt = seq(width*floor(ylim[1]/width -1), width*ceiling(ylim[2]/width+1), by = width)
    #       if (ylim[1] >= 0)
    #           yaxt = seq(width*floor(ylim[1]/width -1), ylim[2], by = width)
    #       if (ylim[2] <= 0)
    #           yaxt = seq(0, ylim[1], by = -width)
    #       if (ylim[1] < 0 && ylim[2] > 0){
    #           temp.1 = seq(0, ylim[1], by = -width)
    #           temp.2 = seq(0, ylim[2], by = width)
    #           yaxt = unique(sort(c(temp.1, temp.2)))
    #           }
            yaxt = yaxt[which(yaxt >= ylim[1] & yaxt <= ylim[2])]
            }
        yaxt = substr(yaxt, 1, yaxt.nchar.max)
        yaxt.len = length(yaxt)
        yaxt = as.character(yaxt[yaxt.len:1])
        }
    # ensure each pch is only one character long
    pch = substr(pch, 1, 1)
    # repeat pch character for the cbind later on
    if (length(pch) == 1)
        pch = rep(pch, length(x))
    # process pch when multiple characters are given
    if (length(pch) != length(x)){
        pch = rep(pch, ceiling(length(x)/length(pch)))
        pch = pch[1:length(x)]
        }

    # subset x and y to fit within xlim and ylim
    x = c(xlim[1], x, xlim[2])
    y = c(ylim[1], y, ylim[2])
    pch = c("?", pch, "?")
    keep.x = which(x >= xlim[1] & x <= xlim[2])
    x = x[keep.x]
    y = y[keep.x]
    pch = pch[keep.x]
    keep.y = which(y >= ylim[1] & y <= ylim[2])
    x = x[keep.y]
    y = y[keep.y]
    pch = pch[keep.y]
    pch = pch[-c(1, length(pch))]

    # calculations to fit within max.width and max.height
    # convert coordinates into the natural numbers
    # the +1 or +2 account for blank space
    left.margin = 0
    bot.margin = 0
    top.margin = 0
    if (nchar(main) != 0)
        top.margin = top.margin + length(main) + 1
    if (nchar(xlab) != 0)
        bot.margin = bot.margin + length(xlab) + 1
    if (nchar(ylab) != 0)
        left.margin = left.margin + nchar(ylab) + 1
    if (box){
        left.margin = left.margin + 1
        top.margin = top.margin + 1
        bot.margin = bot.margin + 1
        }
    if (xtick)
        bot.margin = bot.margin + 2
    if (ytick)
        left.margin = left.margin + max(nchar(yaxt)) + 2
    # the 8 accounts for the space on top and bottom of plot
    # h.range is the range of y values in the plot itself
    h.range = (max.height - (top.margin + bot.margin)):1
    w.range = round(seq(1, max.width - left.margin - 4, length=xaxt.len))
    # m.range is the range of x values in the plot itself
    m.range = 1:(max.width - left.margin - 4)
    # check xaxt characters, the +3 is for the space, tick mark and left side of the plot window
    if (left.margin + 3 - nchar(xaxt[1]) < 0)
        stop("Variable xaxt[1] contains too many characters. Consider rounding.")
    for (j in 2:xaxt.len){
        if (w.range[j] - w.range[j-1] - nchar(xaxt[j]) < 0)
            stop(paste0("Variable xaxt[",j,"] contains too many characters. Consider rounding."))
        }
    # convert x's and y's into the natural numbers for plot coordinates
    new.x = rescale(x, 1, max.width - 4 - left.margin)
    new.x = new.x[-c(1, length(new.x))]
    new.y = rescale(y, 1, max.height - 8)
    new.y = new.y[-c(1, length(new.y))]
    # matrix containing coordinates and which character to use for printing
    P = data.frame(new.x, new.y, pch)
    # remove repeat coordinates (keep the ones that appear first)
    P = P[as.numeric(rownames(unique(P[,1:2]))),]
    # (rownames conveniently returns row number since P is a new data.frame)
 
    # header
    # print main title (assuming center justify, change later)
    if (nchar(main) != 0){
        for (i in 1:length(main))
            cat0(fill(" ",max(0, floor((max.width - nchar(main[i]))/2))), main[i], "\n")
        # print blank line
        cat0("\n")
        }
    # print top line of frame
    if (box)
        cat0(fill(" ", left.margin - 1), "+", fill("-", max.width-4-left.margin), "+\n")


#       cat0(fill(" ", nchar(ylab)+1+max(nchar(yaxt))-nchar(yaxt[1])), yaxt[1],
#           " -+", fill("-", max.width-4-left.margin), "+\n")
    # body (and left margin)
    y.at = round(seq(h.range[1]+1, 0, length=length(yaxt)))
    mid = round((h.range[1]+1)/2)
    for (row in h.range){
        # row is on ylab location and a yaxt location
        if (row == mid && any(row == y.at))
            cat0(ylab, fill(" ", max(nchar(yaxt))-nchar(yaxt[which(row == y.at)])+1),
                yaxt[which(row == y.at)], " -|")
        # row is on ylab location but not on yaxt location
        if (row == mid && !any(row == y.at))
            cat0(ylab, fill(" ", max(nchar(yaxt))+1), "  |")
        # row is on yaxt and not mid
        if (row != mid && any(row == y.at))
            cat0(fill(" ", nchar(ylab)), fill(" ", max(nchar(yaxt))-
                nchar(yaxt[which(row == y.at)])+1), yaxt[which(row == y.at)], " -|")
        # row is on neither ylab or yaxt
        if (row != mid && !any(row == y.at))
            cat0(fill(" ", left.margin), "|")
        # print coordinates
        subP = P[P[,2] %in% row,]
        for (dot in m.range){
            get = which(dot == subP[,1])
            if (length(get) > 0){
                cat0(as.character(subP[get, 3]))
            } else {
                cat0(" ")
                }
            }
        cat0("|\n")
        }
    # print bottom line of frame
    cat0(fill(" ", nchar(ylab)+1+max(nchar(yaxt))-nchar(yaxt[yaxt.len])), yaxt[yaxt.len],
        " -+", fill("-", max.width-4-left.margin), "+\n")
    # print x-axis tick marks
    cat0(fill(" ", left.margin+3), "'")
    for (col in 2:length(w.range))
        cat0(fill(" ", w.range[col] - w.range[col - 1] - 1), "'")
    cat0("\n")
    # print x-values
    cat0(fill(" ", left.margin+4-nchar(xaxt[1])), xaxt[1])
    for (j in 2:xaxt.len)
        cat0(fill(" ", w.range[j] - w.range[j-1] - nchar(xaxt[j])), xaxt[j])
    cat("\n")
        
    # footer
    if (nchar(xlab) != 0){
        # print blank line
        cat0("\n")
        # print xlabel
        cat0(fill(" ",max(0, floor((max.width - nchar(xlab))/2))), xlab, "\n")
        }

    }
