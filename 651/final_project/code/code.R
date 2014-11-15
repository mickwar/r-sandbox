source("functions.R")
field = matrix(scan("~/files/R/651/final_project/data/field.txt"), 9, 4, byrow=TRUE)
simul = matrix(scan("~/files/R/651/final_project/data/simul.txt"), 693, 13, byrow=TRUE)

# original data
field.x = field[,1]
field.y = field[,2:4]

simul.x = simul[,1]
simul.t = simul[,2:10]
simul.y = simul[,11:13]

n = nrow(field)
m = nrow(simul)
px = 1
pt = ncol(simul.t)
p = px + pt

# transform the data
fdat.y = logit(field.y)
sdat.y = logit(simul.y)


xmin = min(field.x)
xmax = max(field.x)
fdat.x = (field.x - xmin) / (xmax - xmin)
sdat.x = (simul.x - xmin) / (xmax - xmin)

tmin = apply(simul.t, 2, min)
tmax = apply(simul.t, 2, max)

sdat.t = simul.t - matrix(rep(tmin, m), m, pt, byrow=TRUE)
sdat.t = sdat.t / (matrix(rep(tmax, m), m, pt, byrow=TRUE) - matrix(rep(tmin, m), m, pt, byrow=TRUE))
sdat.xt = cbind(sdat.x, sdat.t)



0.001^y.dist[[1]]

