source("functions.R")
field = matrix(scan("~/files/R/651/final_project/data/field.txt"), 9, 4, byrow=TRUE)
simul = matrix(scan("~/files/R/651/final_project/data/simul.txt"), 693, 13, byrow=TRUE)

### original data
field.x = field[,1]
field.y = field[,2:4]

simul.x = simul[,1]
simul.t = simul[,2:10]
simul.y = simul[,11:13]

# set constant variables
n = nrow(field)
m = nrow(simul)
px = 1
pt = ncol(simul.t)
p = px + pt
q = ncol(field.y)

### transform the data
fdat.y = logit(field.y)
sdat.y = logit(simul.y)
y.mean = apply(sdat.y, 2, mean)
y.sd = apply(sdat.y, 2, sd)

# standardize output based on simulation mean and sd
fdat.y = (fdat.y - matrix(rep(y.mean, n), n, q, byrow=TRUE)) /
    matrix(rep(y.sd, n), n, q, byrow=TRUE)
sdat.y = (sdat.y - matrix(rep(y.mean, m), m, q, byrow=TRUE)) /
    matrix(rep(y.sd, m), m, q, byrow=TRUE)

### ??
K.eta = svd(t(sdat.y))$u
apply(K.eta, 2, mean)
apply(K.eta, 2, sd) # 1/sqrt(m-1)
### ??

xmin = min(field.x)
xmax = max(field.x)
fdat.x = (field.x - xmin) / (xmax - xmin)
sdat.x = (simul.x - xmin) / (xmax - xmin)

tmin = apply(simul.t, 2, min)
tmax = apply(simul.t, 2, max)

sdat.t = simul.t - matrix(rep(tmin, m), m, pt, byrow=TRUE)
sdat.t = sdat.t / (matrix(rep(tmax, m), m, pt, byrow=TRUE) - matrix(rep(tmin, m), m, pt, byrow=TRUE))
sdat.xt = cbind(sdat.x, sdat.t)

