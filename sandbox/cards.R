### f simulates two shuffled decks of cards (of size n) that counts
### the number of matches. I.e. flip over the top card of each deck
### and count if both cards are the same. Do this for all cards in the deck.
### The result is a draw from a Poisson(1)
f = function(n){
    x = sample(n)
    y = sample(n)
    return (sum(x == y))
#   if (sum(x == y) == 0)
#       return (0)
#   return (1)
    }

n = 20
m = 100000
p = double(m)

for (i in 1:m)
    p[i] = f(n)

plot(table(p)/m, type='h')
xx = seq(0, max(p), by = 1)
points(xx, dpois(xx, 1), col = 'red')

table(p)/m - dpois(xx, 1)

c(mean(p), 1-exp(-1))

1 - (1^0) * exp(-0) / (gamma(0+1))
exp(0)

1 - (1^0) * exp(-1) / gamma(0+1)
1-exp(-1)

mean(p)*(1-mean(p))
