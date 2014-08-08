mu = 15
sigma = 5
iterations = 10000
confidence = 0.99
intervals=matrix(0,iterations,2)
for (i in 1:iterations){
	N=1000
	test = rnorm(N, mu, sigma)
	n = ceiling(N/20)
	randoms = sample(N,n)
	
	test2=test[randoms]
	
	X = mean(test2)
	s = sd(test2)/sqrt(n)
	
	intervals[i,] = X + qt(c((1-confidence)/2, (1+confidence)/2), n-1) * s
	}

(sum((intervals - mu)[,2]<0) + sum((intervals - mu)[,1]>0)) / iterations
# Should be about 1 - confidence

int.lengths = intervals[,2]-intervals[,1]
