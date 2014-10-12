### bootstrap importance

# (take code from mp2)

### metropolis-hastings
autotune = function(accept, target = 0.25, k = 2.5)
    (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
        ceiling(accept-target))-1))^sign(accept-target)
window = 200

f.star = function(x)
    (1 + (x - 10)^2/3)^(-2)

xx = seq(0, 20, length=100)
plot(xx, f.star(xx), type='l')

nmcmc = 10000
nburn = 5000
params = double(nmcmc+nburn)
accept = double(nmcmc+nburn)
sig = 5.90

curr.post = log(f.star(params[1]))
cand.post = curr.post

for (i in 2:(nmcmc+nburn)){
    params[i] = params[i-1]
    cand = rnorm(1, params[i-1], sig)
    cand.post = log(f.star(cand))
    if (log(runif(1)) < cand.post - curr.post){
        curr.post = cand.post
        accept[i] = 1
        params[i] = cand
        }
    if (i <= nburn && floor(i/window) == i/window)
        sig = sig * autotune(mean(accept[(i-window+1):i]),
            target = 0.25, k = window/50)
    }

params = params[(nburn+1):(nburn+nmcmc)]
accept = accept[(nburn+1):(nburn+nmcmc)]

plot(params, type='l')
mean(accept)
sig

xx = seq(0, 20, length=100)
plot(xx, f.star(xx), type='l', lwd=3)
points(density(params,bw =((4*sd(params)^5)/(3*length(params)))^(1/5)),
    type='l', col='red', lwd=3)
