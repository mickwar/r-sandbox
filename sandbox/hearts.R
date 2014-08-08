cards = c(  "02c", "03c", "04c", "05c", "06c", "07c", "08c",
            "09c", "10c", "11c", "12c", "13c", "14c",
            "02d", "03d", "04d", "05d", "06d", "07d", "08d",
            "09d", "10d", "11d", "12d", "13d", "14d",
            "02s", "03s", "04s", "05s", "06s", "07s", "08s",
            "09s", "10s", "11s", "12s", "13s", "14s",
            "02h", "03h", "04h", "05h", "06h", "07h", "08h",
            "09h", "10h", "11h", "12h", "13h", "14h")

# 11 = Jack
# 12 = Queen
# 13 = King
# 14 = Ace

deal = function(cards){
    shuffle = sample(cards)
    H1 = shuffle[1:13]
    H2 = shuffle[14:26]
    H3 = shuffle[27:39]
    H4 = shuffle[40:52]
    return(list(H1, H2, H3, H4))
    }

trick = function(hands, starter, broken.hearts=FALSE){
    played = character(4)
    num.trick = sum(is.na(hands[[1]]))
    if (num.trick == 0){
        for (i in 1:4){
            if (any(hands[[i]] == "02c"))
                starter = i
            }
        }
    for (i in (1+((starter + (-1:2)) %% 4))){
        if (i == starter){
            if (num.trick == 0){
                valid.play = "02c"
            } else {
                if (broken.hearts){
                    # Able to start with any card
                    valid.play = hands[[i]][1:(13-num.trick)]
                } else {
                    # Start with non-hearts, unless only hearts available
                    # 
                    non.hearts = which(substring(hands[[i]], 3, 3) != "h")
                    if (length(non.hearts) > 0){
                        valid.play = hands[[i]][non.hearts]
                    } else {    # no hearts
                        valid.play = hands[[i]][1:(13-num.trick)]
                        }
                    }
                }
        } else {
            suit.played = substring(played[starter], 3, 3)
            suit.index = which(substring(hands[[i]], 3, 3) == suit.played)
            if (length(suit.index) > 0){
                valid.play = hands[[i]][suit.index]
            } else {    # no cards in leading suit
                valid.play = hands[[i]][1:(13-num.trick)]
                }
            # Can't break hearts or play Queen of Spades on first trick
            if (num.trick == 0){
                valid.play = valid.play[valid.play != "12s"]
                valid.play = valid.play[substring(valid.play, 3, 3) != "h"]
                }
            }
        # Play random card from valid cards to play
        played[i] = sample(valid.play, 1)
        # Replace position of played card with last non-NA and set
        # last non-NA to NA.
        at = which(hands[[i]] == played[i])
        hands[[i]][at] = hands[[i]][13-num.trick]
        hands[[i]][13-num.trick] = NA
        }
    return (list(played, hands))
    }

simulate.round = function(hands){
    trick.winner = 1
    broke.hearts = 0
    points = c(0, 0, 0, 0)
    tricks.won = c(0, 0, 0, 0)
    moon = c(0, 0, 0, 0)
    for (i in 1:4){
        if (any(hands[[i]] == "02c"))
            trick.winner = i
        }
    for (i in 1:13){
        out = trick(hands, trick.winner, broke.hearts)
        played = out[[1]]
        hands = out[[2]]
        # Get winner of trick
        # Get suit played
        suit.played = substring(played[trick.winner], 3, 3)
        # Subset played cards into same suit as suit played, then find maximum
        same.suit.cards = played[substring(played, 3, 3) == suit.played]
        winning.card = same.suit.cards[which.max(as.numeric(substring(same.suit.cards, 1, 2)))]
        trick.winner = which(played == winning.card)
        # Add up points, number of hearts, and queen of spades
        points[trick.winner] = points[trick.winner] + sum(substring(played, 3, 3) == "h") + 13*sum(played == "12s")
        tricks.won[trick.winner] = tricks.won[trick.winner] + 1
        }
    if (any(points == 26)){
        moon.shooter = which(points == 26)
        points = c(26, 26, 26, 26)
        points[moon.shooter] = 0
        moon[moon.shooter] = 1
        }
    return (list(points, tricks.won, moon))
    }

hands = deal(cards)
simulate.round(hands)[[1]]

##### Main simulation
### Data to collect
### Response variables
#   Points, tricks won, hand won, moon shot
### X variables
#   Max card value, min card value, sum of card value, has queen of spades,
#   has two of clubs, n clubs, n diamonds, n spades, n hearts


### Collect X variables
get.xs = function(hands){
    out = matrix(0, 4, 9)
    suit = c("c", "d", "s", "h")
    for (i in 1:4){
        out[i, 1] = max(as.numeric(substring(hands[[i]], 1, 2)))
        out[i, 2] = min(as.numeric(substring(hands[[i]], 1, 2)))
        out[i, 3] = sum(as.numeric(substring(hands[[i]], 1, 2)))
        out[i, 4] = sum(hands[[i]] == "12s")
        out[i, 5] = sum(hands[[i]] == "02c")
        for (j in 1:4)
            out[i, 5+j] = sum(substring(hands[[i]], 3, 3) == suit[j])
        }
    return (out)
    }

#library(foreach)
#library(doMC)
#registerDoMC(2)

#useParallel = TRUE

nhands = 10    # number of hands to use
nplays = 100    # number of times to play per hand

dat = matrix(0, 4*nhands, 12)
#engine = function(){
#    for (j in 1:nplays){
#        simulate.round(hands)
for (i in 1:nhands){
    hands = deal(cards)
    points = matrix(0, nplays, 4)
    tricks = matrix(0, nplays, 4)
    moons = matrix(0, nplays, 4)
    at = seq(1+4*(i-1), 4*i)
    dat[at, 4:12] = get.xs(hands)
    for (j in 1:nplays){
        out = simulate.round(hands)
        points[j, ] = out[[1]]
        tricks[j, ] = out[[2]]
        moons[j, ] = out[[3]]
        }
    dat[at, 1] = apply(points, 2, mean)
    dat[at, 2] = apply(tricks, 2, mean)
    dat[at, 3] = apply(moons, 2, mean)
    }
colnames(dat) = c("avg_point", "avg_trick", "avg_moon", "max", "min", "sum",
   "has_QS", "has_2C", "n_club", "n_diam", "n_spad", "n_hear") 
moons = dat[dat[,3]>0,]
apply(moons, 2, mean)

Y = dat[, 1]
X = dat[, 2:12]
#X = X[,-2]
n = length(Y)
k = ncol(X)

betahat = solve(t(X) %*% X) %*% t(X) %*% Y
s2 = 1/(n-k)*t(Y- X %*% betahat) %*% (Y - X %*% betahat)
err = sqrt(diag(as.vector(s2)*solve(t(X) %*% X)))
t.stat = betahat/err
p.val = 2*pt(-abs(t.stat), n-k)
est = cbind("mean"=as.vector(betahat), "std.err"=as.vector(err),
    "t.stat"=as.vector(t.stat), "p.val"=as.vector(p.val))
rownames(est) = rownames(betahat)

est

writeLines(as.character(Y), "./hearts_Y.txt")
write.table(X, "./hearts_X.txt")
write.table(est, "./hearts_estimate.txt")
