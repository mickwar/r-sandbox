###########
### MPG example
library(date)
stdev = function(x) sqrt((length(x)-1)/length(x)) * sd(x)

# mpg data
dat = read.table("~/files/data/MPG_saturn.txt", header = TRUE)
y = dat$miles / dat$gallons
#y = y[-c(13, 20)] # removing the most extreme values

### Simple linear regression
t = as.numeric(as.date(gsub("-","", dat[,1])))
x = diff(t)
y = y[-1]
t = t[-1]

plot(as.date(t), y, type='b', pch = 20, ylab = "MPG", xlab = "Date")

plot(x, y, pch = 20, xlab = "Days since previous fill up", ylab = "MPG", cex.lab = 1.2)


ord = order(x)
mod = lm(y[ord] ~ x[ord]) #Some ad hoc grouping, I figure <10 was road trip
summary(mod)
# predict on data
plot(x, y, pch=20, main = "MPG", xlab = "Days since previous fill up", ylab = "MPG", cex.lab = 1.3); lines(x[ord], predict(mod), lwd = 3)

#coef(mod)
#cor(x, y) * sd2(y)
sqrt(1-cor(x, y)^2) * stdev(y)

sqrt(mean((y[ord] - predict(mod))^2))

# residuals
plot(rstudent(mod), pch = 20, main = "Standardized Residuals"); abline(h = 0, lwd = 2)

# fitted vs observed
#plot(predict(mod), y[ord], pch = 20, main = "Fitted vs Observed"); abline(0, 1)

### Natural splines (number of days between fill ups to predict MPG)
library(splines)
ord = order(x)
mod = lm(y[ord] ~ ns(x[ord], knots= c(10, 30, 60))) #Some ad hoc grouping, I figure <10 was road trip
deviance(mod)
summary(mod)

sqrt(mean((y[ord] - predict(mod))^2)*52/47)

# predict on data
plot(x, y, pch=20, main = "MPG", xlab = "Days since previous fill up", ylab = "MPG"); lines(x[ord], predict(mod), lwd = 3)

# residuals
plot(rstudent(mod), pch = 20, main = "Standardized Residuals"); abline(h = 0, lwd = 2)

# fitted vs observed
#plot(predict(mod), y[ord], pch = 20, main = "Fitted vs Observed"); abline(0, 1)

###########
# Ch 9 number 4
library(MASS)

men.mean = c(70, 144)
women.mean = c(64, 120)

set.seed(2)
men = mvrnorm(50, men.mean, matrix(c(3, 4.762, 4.762, 21), 2, 2))
women = mvrnorm(50, women.mean, matrix(c(3, 4.762, 4.762, 21), 2, 2))

# men points
plot(men.mean[1], men.mean[2], pch = 20, xlim = c(55, 75), ylim = c(100, 170), col = 'blue', cex = 4)
points(women.mean[1], women.mean[2], pch = 20, col = 'red', cex = 4)

abline(a=men.mean[2] - men.mean[1] * (0.6 * 21/3), b = (0.6 * 21/3),
    col = 'lightblue3', lty=2, lwd = 2)
abline(a=women.mean[2] - women.mean[1] * (0.6 * 21/3), b = (0.6 * 21/3),
    col = 'pink2', lty=2, lwd = 2)

points(men, col = 'blue', pch = 20)
points(women, col = 'red', pch = 20)


# correlation separately
cor(men)[1,2]
cor(women)[1,2]

# correlation together
cor(rbind(men, women))[1,2]


##########
# Ch 9 number 9
assistant = c(3.3, 2.9, 4.1, 3.3, 2.7, 3.4, 2.8, 2.1, 3.7, 3.2, 2.4)
course = c(3.5, 3.2, 3.1, 3.3, 2.8, 3.5, 3.6, 2.8, 2.8, 3.3, 3.3)
final = c(70, 64, 47, 63, 69, 69, 69, 63, 53, 65, 64)

cor(assistant, final)

cor(assistant, course)

cor(course, final)

plot(assistant, final, pch = 20)
plot(course, final, pch = 20)
plot(assistant, course, pch = 20)

mod = lm(final ~ assistant + course)
summary(mod)

