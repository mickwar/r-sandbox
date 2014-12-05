movie.dat = read.table("~/files/R/data/movie_data.txt", header=T)
ben = read.table("~/files/R/data/movie_ben.txt", header=T)
ben[,1] = as.character(ben[,1])
movies = ben[,1]

# combine reponses into one category
new = double(length(movies))
new[apply(ben[,2:3], 1, function(x) all(x == c(1, 1)))] = 1 # liked and well executed
new[apply(ben[,2:3], 1, function(x) all(x == c(1, 0)))] = 2 # liked but not well executed
new[apply(ben[,2:3], 1, function(x) all(x == c(0, 1)))] = 3 # disliked, but well executed
new[apply(ben[,2:3], 1, function(x) all(x == c(0, 0)))] = 4 # disliked and not well executed
ben = ben[,-3]
ben[,2] = new

# remove the classes with few observations (i.e., liked but not well executed only has 2)
ben = ben[-which(ben[,2] == 2),]

movie.dat[,1] = as.character(movie.dat[,1])
movie.dat$MPAA = as.character(movie.dat$MPAA)

n.var = ncol(movie.dat) - 1
dat = data.frame(matrix(0, nrow(ben), n.var + ncol(ben)))
dat[,1:2] = ben[,2:1]
not_found = NULL
for (i in 1:nrow(ben)){
    j = which(ben[i,1] == movie.dat[,1])
    if (length(j) == 0){
        not_found[length(not_found)+1] = ben[i,1]
    } else {
        dat[i,3:(2+n.var)] = movie.dat[j,2:(n.var+1)]
        }
    }

names(dat) = c("ben", names(movie.dat))

# remove nas
dat = dat[-which(apply(dat, 1, function(x) any(is.na(x)))),]

dat = cbind(dat, "kids_I"=is.na(dat$kids_S)*1)
for (i in 3:5)
    dat[,i] = ifelse(is.na(dat[,i]), 0, dat[,i])

# functions
# if drop = TRUE, drop the last column
mpaa.matrix = function(x){
    level = c("PG13", "PG", "G")
    # NR is the dropped level
    n = length(x)
    k = length(level)
    out = matrix(0, n, k)
    for (i in 1:k)
        out[which(x == level[i]), i] = 1
    return (out)
    }
Y = dat[,1]
raw.X = as.matrix(cbind(dat[,3:5], mpaa.matrix(dat$MPAA),
    dat[,7:15]))
X = cbind(raw.X[,1:6], 100*(raw.X[,7]-raw.X[,9])/(raw.X[,7]+
    raw.X[,8]-raw.X[,9]-raw.X[,10]), 100*raw.X[,9]/(raw.X[,9]+ 
    raw.X[,10]), raw.X[,11:15])

# top critic rotten tomatoes indicator
X = cbind(X, "rot_I"=is.na(X[,8])*1)
X[,8] = ifelse(is.na(X[,8]), 0, X[,8])

var.names = c(names(dat)[3:5], as.character(unique(dat[,6])[1:3]),
    "rot_crit", "rot_top", names(dat)[11:15], "rot_I")
colnames(X) = var.names
rownames(X) = as.character(1:nrow(X))

X = X[,-c(8, 9, 11, 13, 14)]
