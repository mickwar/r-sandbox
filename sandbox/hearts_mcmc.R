logit = function(x)
    log(x/(1-x))
ilogit = function(x)
    1/(1+exp(-x))

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

trick = function(hands, starter, params, broken.hearts=FALSE){
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
        values = as.numeric(substring(valid.play, 1, 2))
        suit = substring(valid.play, 3, 3)
        suit.num = matrix(0, length(suit), 4)
        for (j in 1:length(suit)){
            if (suit[j] == "c")
                suit.num[j,1] = 1
            if (suit[j] == "d")
                suit.num[j,2] = 1
            if (suit[j] == "s")
                suit.num[j,3] = 1
            if (suit[j] == "h")
                suit.num[j,4] = 1
            }

        X = cbind(1, values, suit.num)

        played[i] = sample(valid.play, 1, prob = ilogit(X %*% params))
        # Replace position of played card with last non-NA and set
        # last non-NA to NA.
        at = which(hands[[i]] == played[i])
        hands[[i]][at] = hands[[i]][13-num.trick]
        hands[[i]][13-num.trick] = NA
        }
    return (list("played"=played, "hands"=hands))
    }


simulate.round = function(hands, beta){
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
        out = trick(hands, trick.winner, beta, broke.hearts)
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
    return (list("points"=points, "ntricks"=tricks.won, "moon"=moon))
    }


# init beta
# beta[1] is intercept (is this necessary?, probably not, keep for now)
# beta[2] is coefficient for point value on card
# beta[3:6] is coefficient for suit indicator

# this isn't really mcmc. to use the metropolis algorithm,
# need to actually construct a likelihood with priors
nmcmc = 5000
nparams = 6
params = matrix(0, nmcmc, nparams)
accept = matrix(0, nmcmc, nparams)
sig = rep(0.1, nparams)


hands = deal(cards)
for (i in 2:nmcmc){
    cat("\rIteration:",i,"/",nmcmc)
    params[i,] = params[i-1,]
    cand.params = params[i,]
    for (j in 1:nparams){
        cand.params[j] = rnorm(1, params[i,j], sig[j])
        if (which.min(simulate.round(hands, cand.params)$points) == 1){
            params[i,j] = cand.params[j]
            accept[i,j] = 1
        } else {
            cand.params[j] = params[i,j]
            }
        }
    if (i == nmcmc)
        cat("\n")
    }

for (i in 1:nparams){
    plot(density(params[,i]))
    abline(v=0)
    if (i != nparams)
        readline()
    }

apply(accept, 2, mean)
mean(accept)

pred = double(nmcmc)
for (i in 1:nmcmc)
    pred[i] = 1*(which.min(simulate.round(hands, params[i,])$points) == 1)

mean(pred)
mean(pred[4501:5000])
plot(pred, pch=20)
