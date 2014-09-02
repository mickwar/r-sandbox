# for integer n, return each prime factor
factors=function(n){
	n=abs(n)
	out=NULL
	new=0
	if (floor(n) != ceiling(n))
		return (factors(floor(n)))
	for (i in 2:n){
		if (n/i == floor(n/i)){
			new=n/i
			break
			}
		}
	if (new==1){
		return (i)
		} else {
			return (c(i, factors(new)))
		} 
	}
factors(223092870)
factors(1073741824)

primes=NULL
upto=10000
for (i in 2:upto){
	if (factors(i)[1]==i)
		primes[length(primes)+1]=i
	}
