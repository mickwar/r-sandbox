# By Arthur Lui
# luiarthur@gmail.com

color.den <- function(den,from,to) {
  # Colors area under a density within an interval
  # den has to be a density object
  polygon(c(0, den$x[den$x> from & den$x < to], to), 
          c(0, den$y[den$x>=from & den$x <= to], 0),col="red")
}

color.fn <- function(f,from,to) {
  x <- seq(from,to,by=(to-from)/10000)
  polygon(c(0, x,    to), 
          c(0, f(x), 0),col="red")
}


color.emp <- function(x,y,from,to) {
  polygon(c(0, x[x> from & x < to],    to), 
          c(0, y[x> from & x < to], 0),col="red")
}

#Examples: ######################################################

## color.den
#  x <- rnorm(10000)
#  denx <- density(x)
#  plot(denx)
#  color.den(denx,0,5)
#
## color.fn
#  fn <- function(x) dnorm(x)
#  curve(dnorm(x),-5,5)
#  color.fn(fn,0,5)
#
## color.emp
#  x <- seq(-5,5,length=10)
#  y <- dnorm(x)
#  plot(x,y,type='l')
#  color.emp(x,y,0,5)
