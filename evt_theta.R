set.seed(1)
n = 50000
Wn = 1/(-log(runif(n)))
theta = 0.10

Y = double(n)
Y[1] = Wn[1]/theta
for (i in 2:n)
    Y[i] = max( (1-theta)*Y[i-1], Wn[i] )

# plot(log(Y), type = 'l')
u = quantile(Y, 0.90)

ex.time = which(Y > u)



# N = 12
# #inter = c(1, 11, Inf, 2, 8, 1, 1, Inf)
# 
# ex.time = c(1, 2, 13, 33, 35, 43, 44, 45, 65, 70, 72, 74)


N = length(ex.time)
inter = diff(ex.time)

#inter = c(1, 11, 20, 2, 8, 1, 1, 20)

#theta = 0.5
theta.hat = min(1, (2*sum(inter - 1)^2) / ((N-1)*sum((inter - 1)*(inter - 2))))
C = floor(theta.hat*N)+1
#T_C = min(inter[inter >= C])
tmp = sort(inter, decreasing = TRUE)
T_C = tmp[C]
while (!(tmp[C-1] > T_C)){
    C = C - 1
    T_C = tmp[C]
    }



# sort(inter) > T_C

inter.Clust = inter[inter > T_C]
# inter
# 
# i_j = c(0, 5, 2, 3, 8, N)
# 
# inter > T_C

i_j = which(inter > T_C)
i_j = c(0, i_j, N)
ind.seq = rep(list(NULL), C)
intra.Clust = rep(list(NULL), C)
nice.S = rep(list(NULL), C)
nice.C = rep(list(NULL), C)

for (j in 2:(C+1)){
    ind.seq[[j-1]] = seq(i_j[j-1]+1, i_j[j]-1)
    if (i_j[j-1]+1 == i_j[j]){
#       nice.T[[j-1]] = NULL
    } else {
        intra.Clust[[j-1]] = inter[seq(i_j[j-1]+1, i_j[j]-1)]
        }
    nice.S[[j-1]] = ex.time[seq(i_j[j-1]+1, i_j[j])]
    nice.C[[j-1]] = Y[nice.S[[j-1]]]
    }

inter.Clust
intra.Clust
nice.S



### Bootstrap
theta.vec = double(20000)
for (i in 1:length(theta.vec)){

    samp.inter = sample(C-1, replace = TRUE)
    samp.intra = sample(C, replace = TRUE)

    tmp = c(inter.Clust[samp.inter], unlist(intra.Clust[samp.intra]))
    theta.vec[i] = min(1, (2*sum(tmp - 1)^2) / ((N-1)*sum((tmp - 1)*(tmp - 2))))
    }

plot(density(theta.vec))
abline(v = mean(theta.vec))
abline(v = min(1, (2*sum(inter - 1)^2) / ((N-1)*sum((inter - 1)*(inter - 2)))), lty = 2)

