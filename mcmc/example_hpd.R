source("./density.R")

# draws from three distributions, x1 is bimodal, x2, x3 are unimodal
n = 200000
test.x1 = double(n)
for (i in 1:n)
    test.x1[i] = ifelse(runif(1) < 0.5, rnorm(1), rnorm(1, 6))
test.x2 = double(n)
for (i in 1:n){
    j = sample(3, 1)
    test.x2[i] = rnorm(1, c(-8,0,8)[j], c(2,0.5,1.5)[j])
    }
test.x3 = rbeta(n, 0.5, 0.5)
test.x4 = double(n)
for (i in 1:n)
    test.x4[i] = ifelse(runif(1) < 0.900, rnorm(1, -5), rnorm(1, 5))

# intervals from each function
uni1 = hpd.uni(test.x1)
uni2 = hpd.uni(test.x2)
uni3 = hpd.uni(test.x3)
uni4 = hpd.uni(test.x4)
mult1 = hpd.mult(test.x1)
mult2 = hpd.mult(test.x2)
mult3 = hpd.mult(test.x3)
mult4 = hpd.mult(test.x4)
# hpd.mult() should result in a probability region slightly
# less than that specified by the prob argument, but the area
# contained in the set is no greater than that of prob

# note that because of the density fuction, the intervals produced
# by hpd.mult() may be outside the support for the data. fix by
# imposing the constraints on the output of hpd.mult(), the
# probabilities computing in the function will remain the same.
mult3[1] = max(0, mult3[1])
mult3[length(mult3)] = min(1, mult3[length(mult3)])

# interval/set lengths
uni1[2] - uni1[1]
sum(mult1[seq(2, length(mult1), by=2)])-
    sum(mult1[seq(1, length(mult1), by=2)])

uni2[2] - uni2[1]
sum(mult2[seq(2, length(mult2), by=2)])-
    sum(mult2[seq(1, length(mult2), by=2)])

uni3[2] - uni3[1]
sum(mult3[seq(2, length(mult3), by=2)])-
    sum(mult3[seq(1, length(mult3), by=2)])

uni4[2] - uni4[1]
sum(mult4[seq(2, length(mult4), by=2)])-
    sum(mult4[seq(1, length(mult4), by=2)])

# probability within interval/set
mean(test.x1 >= uni1[1] & test.x1 <= uni1[2])
mean1 = 0
for (i in 1:(length(mult1)/2))
    mean1 = mean1 + mean(test.x1>=mult1[2*i-1] & test.x1<=mult1[2*i])
mean1

mean(test.x2 >= uni2[1] & test.x2 <= uni2[2])
mean2 = 0
for (i in 1:(length(mult2)/2))
    mean2 = mean2 + mean(test.x2>=mult2[2*i-1] & test.x2<=mult2[2*i])
mean2

mean(test.x3 >= uni3[1] & test.x3 <= uni3[2])
mean3 = 0
for (i in 1:(length(mult3)/2))
    mean3 = mean3 + mean(test.x3>=mult3[2*i-1] & test.x3<=mult3[2*i])
mean3

mean(test.x4 >= uni4[1] & test.x4 <= uni4[2])
mean4 = 0
for (i in 1:(length(mult4)/2))
    mean4 = mean4 + mean(test.x4>=mult4[2*i-1] & test.x4<=mult4[2*i])
mean4

# plots
par(mfrow=c(2,2), mar=c(3.6,2.6,2.1,0.6))
plot(density(test.x1), xlab="", ylab="")
abline(v=uni1, col='red')
abline(v=mult1, col='blue')

plot(density(test.x2), xlab="", ylab="")
abline(v=uni2, col='red')
abline(v=mult2, col='blue')

plot(density(test.x3), xlab="", ylab="")
abline(v=uni3, col='red')
abline(v=mult3, col='blue')

plot(density(test.x4), xlab="", ylab="")
abline(v=uni4, col='red')
abline(v=mult4, col='blue')
