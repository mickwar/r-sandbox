iter.fact=function(numerator, denominator, out=1){
	if (missing(denominator))
		denominator=1
	while (numerator > 1 || denominator > 1){
		out = out * numerator / denominator
		if (numerator > 1)
			numerator = numerator - 1
		if (denominator > 1)
			denominator = denominator - 1
		}
	return (out)
	}

birthday=function(n){
	i=365
	p=1
	while (n > 1){
		p=p*(i-1)/365
		i=i-1
		n=n-1
		}
	return (1-p)
	}

same.simu=function(n, x, iter=1000, show.birthdays=FALSE){
	# n = number of people in group
	# x = number of people whose birthday to share (minimum)
	if (x < 2 || is.na(x))
		x=0
	x=floor(x)
	birthdays=list()
	out = integer(iter)
			# 0, no matching birthdays for that simulation
			# >1, number of matching birthday found
	for (i in 1:iter){
		birthdays[[i]] = sample(1:365, n, TRUE)
		count = 0
		for (person in 1:(n-1)){
			match = 1
			for (compare in (person+1):n){
				if (birthdays[[i]][person]==birthdays[[i]][compare])
					match = match + 1
				}
			if (match >= x)
				count = count + 1
			}
		out[i]=count
		}
#	if (show.birthdays)
#		print (birthdays)
	return (out)
	}

out=same.simu(13,2,10000, TRUE)
mean(out>0)
plot(density(out))
hist(out,col='gray',freq=F,breaks=0:10)
points(0:10,dpois(0:10,mean(out)),type='o',col='red')

### Looks like the probability of number of birthdays follows
### a Poisson distribution
