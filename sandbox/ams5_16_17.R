# Coin flipping demo (discrete probability distribution example)
par(mfrow=c(3,2), mar = c(2.1, 2.1, 2.1, 1.1))
plot(0:1, dbinom(0:1, 1, 0.5), type='h', ylim = c(0, 1), lwd = 2)
plot(0:2, dbinom(0:2, 2, 0.5), type='h', ylim = c(0, 1), lwd = 2)
plot(0:3, dbinom(0:3, 3, 0.5), type='h', ylim = c(0, 1), lwd = 2)
plot(0:4, dbinom(0:4, 4, 0.5), type='h', ylim = c(0, 1), lwd = 2)
plot(0:5, dbinom(0:5, 5, 0.5), type='h', ylim = c(0, 1), lwd = 2)
plot(0:6, dbinom(0:6, 6, 0.5), type='h', ylim = c(0, 1), lwd = 2)


# "Biased coin"
par(mfrow=c(3,2), mar = c(2.1, 2.1, 2.1, 1.1))
plot(0:1, dbinom(0:1, 1, 0.2), type='h', ylim = c(0, 1), lwd = 2)
plot(0:2, dbinom(0:2, 2, 0.2), type='h', ylim = c(0, 1), lwd = 2)
plot(0:3, dbinom(0:3, 3, 0.2), type='h', ylim = c(0, 1), lwd = 2)
plot(0:4, dbinom(0:4, 4, 0.2), type='h', ylim = c(0, 1), lwd = 2)
plot(0:5, dbinom(0:5, 5, 0.2), type='h', ylim = c(0, 1), lwd = 2)
plot(0:6, dbinom(0:6, 6, 0.2), type='h', ylim = c(0, 1), lwd = 2)


### Chapter 16
# Number 1
dbinom(6000, 10000, 0.6)


# Number 4
par(mfrow=c(2,1), mar=c(4.1,3.1,3.1,1.1))

p = 1/6
n = 60
plot(0:n, dbinom(0:n, n, p), type='h')
abline(v=c(n*0.15, n*0.20)+c(-0.2,0.2), col = c("blue", "red"), lty = 2)
n = 600
plot(0:n, dbinom(0:n, n, p), type='h')
abline(v=c(n*0.15, n*0.20)+c(-0.2,0.2), col = c("blue", "red"), lty = 2)

# a
1-pbinom(60*0.2, 60, p)
1-pbinom(600*0.2, 600, p)

# b
1-pbinom(60*0.15, 60, p)
1-pbinom(600*0.15, 600, p)

# c
diff(pbinom(60*c(0.15, 0.2), 60, p))
diff(pbinom(600*c(0.15, 0.2), 600, p))

# d
dbinom(60*0.2, 60, p)
dbinom(600*0.2, 600, p)

# Number 9 (assuming p = 0.55)
1-pbinom(50, 100, 0.55)     # Probability of drawing 51 or more red (i.e. more blue)
1-pbinom(100, 200, 0.55)    # Probability of drawing 101 or more red (i.e. more blue)

### Chapter 17
# Number 1
x = c(1, 6, 7, 9, 9, 10)
# a
c(100 * min(x), 100 * max(x))

# b
sqrt(100) * sqrt(sum((x - mean(x))^2 / length(x))) 

# Number 6
45*1 + 2*23 + 3*32

x = c(1,1,2,3)
sqrt(100) * sqrt(sum((x - mean(x))^2 / length(x))) 

# Number 7
x = c(1,2,3,4,5,6)
sqrt(100) * sqrt(sum((x - mean(x))^2 / length(x))) 

(400 - 350) / (sqrt(100) * sqrt(sum((x - mean(x))^2 / length(x))))

# Number 9
# a
# expected winning for $1 bets
(ex = c(12/38*2 + (-1)*26/38, 1/38*35 + (-1)* 37/38))

# standard errors (for 1 bet)
p = 12/38
a = 2
b = -1
se1 = sqrt(a^2*p+b^2*(1-p) - (a*p + b*(1-p))^2)

p = 1/38
a = 35
b = -1
se2 = sqrt(a^2*p+b^2*(1-p) - (a*p + b*(1-p))^2)

# probabilities of coming out ahead
pnorm((0 - ex[1]) / se1, lower.tail = FALSE)
pnorm((0 - ex[2]) / se2, lower.tail = FALSE)



# b and c
# expected gain for $100 bet
ex = c(12/38*2 + (-1)*26/38, 1/38*35 + (-1)* 37/38) * 1000

x = c(rep(2, 12), rep(-1, 26))
se1 = sqrt(1000) * sqrt(sum((x - mean(x))^2 / length(x))) 

x = c(35, rep(-1, 37))
se2 = sqrt(1000) * sqrt(sum((x - mean(x))^2 / length(x))) 

# probability of winning more than $100
pnorm((100 - ex[1]) / se1, lower.tail = FALSE)
pnorm((100 - ex[2]) / se2, lower.tail = FALSE)

# proability of losing more than $100
pnorm((-100 - ex[1]) / se1)
pnorm((-100 - ex[2]) / se2)

# Number 12
x = 1:7
sqrt(100) * sqrt(sum((x - mean(x))^2 / length(x))) 
