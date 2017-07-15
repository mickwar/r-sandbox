library(doMC)
library(foreach)
registerDoMC(4)


par_func = function(i)
    runif(10)^i

foreach(i = 0:5, .combine = rbind) %dopar% par_func(i)

