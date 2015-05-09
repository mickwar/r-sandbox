do.test = function(p){
    chars = c("`", "1", "2", "3", "4", "5", "6", "7",
        "8", "9", "0", "-", "=", "q", "w", "e", "r",
        "t", "y", "u", "i", "o", "p", "[", "]", "\\",
        "a", "s", "d", "f", "g", "h", "j", "k", "l",
        ";", "'", "z", "x", "c", "v", "b", "n", "m",
        ",", ".", "/")

    n = length(chars)
    if (missing(p)) # number of reps for each character
        p = 10
    times = matrix(0, n, p)
    correct = times
    rand.order = sample(rep(1:n, p))
    readline("Hit enter to begin")
    for (i in 1:(n*p)){
        y = system.time(x <- readline(chars[rand.order[i]]))[3]
        k = which.min(times[rand.order[i],])
        correct[rand.order[i], k] = (x == chars[rand.order[i]])*1
        times[rand.order[i], k] = y
        }
    return (list("times"=times, "correct"=correct))
    }

x = do.test()

# overall correct %
mean(x$correct)

# mean time for each correct input
means = double(47)
for (i in 1:47){
    for (j in 1:10){
        means[i] = means[i] + x$times[i,j]*x$correct[i,j]
        }
    }
means = means / apply(x$correct, 1, sum)

# max/min time for correct inputs
maxs = double(47)
mins = double(47)+Inf
for (i in 1:47){
    for (j in 1:10){
        if (x$correct[i,j] == 1){
            maxs[i] = max(maxs[i], x$times[i,j])
            mins[i] = min(mins[i], x$times[i,j])
            }
        }
    }

plot(maxs, type='h', col='red', lwd=2, ylim = c(0, max(maxs)))
points(means, type='h', lwd=2)
points(mins, type='h', col='green', lwd=2)

# cooler plot
key.coors = matrix(
    c(  1, 4, "`",
        2, 4, "1",
        3, 4, "2",
        4, 4, "3",
        5, 4, "4",
        6, 4, "5",
        7, 4, "6",
        8, 4, "7",
        9, 4, "8",
        10, 4, "9",
        11, 4, "0",
        12, 4, "-",
        13, 4, "=",
        2.5, 3, "q",
        3.5, 3, "w",
        4.5, 3, "e",
        5.5, 3, "r",
        6.5, 3, "t",
        7.5, 3, "y",
        8.5, 3, "u",
        9.5, 3, "i",
        10.5, 3, "o",
        11.5, 3, "p",
        12.5, 3, "[",
        13.5, 3, "]",
        14.5, 3, "\\",
        2.8, 2, "a",
        3.8, 2, "s",
        4.8, 2, "d",
        5.8, 2, "f",
        6.8, 2, "g",
        7.8, 2, "h",
        8.8, 2, "j",
        9.8, 2, "k",
        10.8, 2, "l",
        11.8, 2, ";",
        12.8, 2, "'",
        3.1, 1, "z",
        4.1, 1, "x",
        5.1, 1, "c",
        6.1, 1, "v",
        7.1, 1, "b",
        8.1, 1, "n",
        9.1, 1, "m",
        10.1, 1, ",",
        11.1, 1, ".",
        12.1, 1, "/"),
    47, 3, byrow = TRUE)
keys = NULL
keys$x = as.numeric(key.coors[,1])
keys$y = as.numeric(key.coors[,2])
keys$c = key.coors[,3]

means = (means - min(means)) / diff(range(means))
maxs = (maxs - min(maxs)) / diff(range(maxs))
mins = (mins - min(mins)) / diff(range(mins))

plot.keys = function(what){
    m = min(what)
    M = max(what)
    what = (what - min(what)) / diff(range(what))
    make.square = function(x, y, fill = "white"){
        for (i in 1:length(x)){
            polygon(c(x[i]-0.5, x[i]-0.5, x[i]+0.5, x[i]+0.5),
                c(y[i]-0.5, y[i]+0.5, y[i]+0.5, y[i]-0.5),
                col = fill[i])
            }
        }
    plot(keys$x, keys$y, pch=keys$c, ylim=c(-1, 8), xlim=c(0, 15.5))
    make.square(keys$x, keys$y, rgb(what, 0.5, 1-what))
    points(keys$x, keys$y, pch=keys$c)
    points(seq(0, 15, length=100), rep(7, 100), pch=20,
        col=rgb(seq(0, 1, length=100), 0.5, seq(1, 0, length=100)))
    text(0.5, 6.5, "Fast")
    text(14.5, 6.5, "Slow")
    }


plot.keys(means)
plot.keys(maxs)
plot.keys(mins)



