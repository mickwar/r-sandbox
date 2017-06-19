dat = read.table("./basketball_scores.txt", sep = ",", header = FALSE)
n = NROW(dat)

# Starting rating for all players
starting_rating = 850

# minimum score
min_score = 600

# like it sounds (involved in calculation of observed and expected scores)
logistic_scale = 5

# For every point_diff rating points that team A is higher
# than team B, team A must score 1 more point to increase
# their rating (involved in calculation of expected score only)
point_diff = 80


# as.character(c(dat[1,1], dat[1,2]))
# as.character(dat[1,1])
# as.character(dat[1,2])

team.A = strsplit(as.character(dat[,1]), " ")
team.B = strsplit(as.character(dat[,2]), " ")

A.score = as.numeric(sapply(team.A, function(x) tail(x, 1)))
B.score = as.numeric(sapply(team.B, function(x) tail(x, 1)))

team.A = lapply(team.A, function(x) x[-length(x)])
team.B = lapply(team.B, function(x) x[-length(x)])

players = c("michael", "tony", "racer", "mickey", "trevor", "austen", "lai", "victor")
scores = matrix(0, n + 1, length(players))
ngames = matrix(0, n + 1, length(players))
nwin = matrix(0, n + 1, length(players))
nloss = matrix(0, n + 1, length(players))
scores[1,] = rep(starting_rating, length(players))
colnames(scores) = players

for (i in 1:n){
    # Carry over the scores
    scores[i+1, ] = scores[i, ]
    scores[i+1, ] = scores[i, ]

    # Carry over number of games
    ngames[i+1, ] = ngames[i, ]
    nwin[i+1, ] = nwin[i, ]
    nloss[i+1, ] = nloss[i, ]

    # Get index for winners and losers
    A.ind = which(players %in% team.A[[i]])
    B.ind = which(players %in% team.B[[i]])

    # Use sum of player scores to determine team score
    R_A = mean(scores[i, A.ind])
    R_B = mean(scores[i, B.ind])

    # Expected score for the winners
    E_A = 1 / (1 + exp((R_B - R_A) / (logistic_scale*point_diff)))

    # Expected score for the losers
    E_B = 1 - E_A
    
    # K-factor (adjusted by number of players on the team)
    K.fac = double(length(players))
    K.fac[c(A.ind, B.ind)] = 48
    for (j in A.ind){
        K.fac[j] = ifelse(scores[i, j] >= 900, 32, K.fac[j])
        K.fac[j] = ifelse(scores[i, j] >= 1200, 16, K.fac[j])
        K.fac[j] = ifelse(scores[i, j] >= 1500, 8, K.fac[j])
        }
    for (j in B.ind){
        K.fac[j] = ifelse(scores[i, j] >= 900, 32, K.fac[j])
        K.fac[j] = ifelse(scores[i, j] >= 1200, 16, K.fac[j])
        K.fac[j] = ifelse(scores[i, j] >= 1500, 8, K.fac[j])
        }

    # Reduce K-factor by size of the team
    K.fac[A.ind] = K.fac[A.ind] / sqrt(length(A.ind))
    K.fac[B.ind] = K.fac[B.ind] / sqrt(length(B.ind))

    # Divide a player's K-factor by 2 when that player has played more than 20 games
    K.fac[A.ind] = K.fac[A.ind] / ifelse(ngames[i, A.ind] > 20, 2, 1)
    K.fac[B.ind] = K.fac[B.ind] / ifelse(ngames[i, B.ind] > 20, 2, 1)

    # Calculate observed as using the score difference (max 10 point diff)
#   S_A = (min(win.score[i] - loss.score[i], 10) + 10) / 20
#   S_B = 1 - S_A
    S_A = 1 / (1 + exp((B.score[i] - A.score[i]) / logistic_scale))
    S_B = 1 - S_A

    # Adjust player scores
    # Team A won or tied (A.score >= B.score), but not by enough (S_A <= E_A),
    # but still don't decrease the score
    if ( (A.score[i] >= B.score[i]) && (S_A <= E_A) )
        K.fac[A.ind] = 0
    # Same for B
    if ( (B.score[i] >= A.score[i]) && (S_B <= E_B) )
        K.fac[B.ind] = 0

    scores[i+1, A.ind] = scores[i, A.ind] + K.fac[A.ind]*(S_A - E_A)
    scores[i+1, B.ind] = scores[i, B.ind] + K.fac[B.ind]*(S_B - E_B)

    # No player's score can drop below 600
    scores[i+1, scores[i+1,] < min_score] = min_score

    ngames[i+1, A.ind] = ngames[i, A.ind] + 1
    ngames[i+1, B.ind] = ngames[i, B.ind] + 1

    if (A.score[i] > B.score[i]){
        nwin[i+1, A.ind] = nwin[i, A.ind] + 1
        nloss[i+1, B.ind] = nloss[i, B.ind] + 1
        }
    if (B.score[i] > A.score[i]){
        nwin[i+1, B.ind] = nwin[i, B.ind] + 1
        nloss[i+1, A.ind] = nloss[i, A.ind] + 1
        }

    }

# Official scores
round(scores)

tail(round(scores), 8)

card = rbind(tail(round(scores), 1), tail(ngames, 1), tail(nwin, 1), tail(nloss, 1))
rownames(card) = c("ELO", "GP", "W", "L")
card = data.frame(card)

# The average score across all players over time (would like to see this increase)
apply(scores, 1, mean)

# Change in player scores over time
apply(scores[,-8], 2, diff)

# For every 20 rank difference, a win by 1 point is considered a draw
1 - (1 / (1 + exp(10 / 4.5)))  # obs
1 - (1 / (1 + exp(200 / (5*20))))   # exp

# 1 - (1 / (1 + exp(d / x)))  = p
# 1-p = (1 / (1 + exp(d / x)))
# 1 / (1-p) =  (1 + exp(d / x))
# 1/(1-p) - 1 =  exp(d / x)
# log(1/(1-p) - 1) = d / x
# x = d / log(1/(1-p)-1)

d = 5
p = 0.7
sc = d / log(1/(1-p)-1)


ss = 0:20
yy = 1 - (1 / (1 + exp(ss / sc)))
plot(ss, yy, type = 'l')
abline(v = ss[which.max(diff(yy) < 0.10) + 1])
abline(v = ss[which.max(diff(yy) < 0.05) + 1])
abline(v = ss[which.max(diff(yy) < 0.01) + 1])
abline(h = seq(0.5, 0.9, by = 0.1), lty = 2)

random.team = function(players, card){
    n = length(players)
    k = choose(n, n/2) / 2    # Number of possible teams

    require(gtools)
    teamA = head(combinations(n, n/2), k)
    teamB = tail(combinations(n, n/2), k)[k:1,]
    rownames(teamB) = NULL

    if (any(apply(cbind(teamA, teamB), 1, function(x) length(unique(x))) != n))
        stop("problem with getting team combinations")

    scA = apply(teamA, 1, function(y) mean(unlist(card[players[y]][1,])))
    scB = apply(teamB, 1, function(y) mean(unlist(card[players[y]][1,])))
    
    dd = abs(scA - scB)
    pp = 1 - (dd - min(dd)) / diff(range(dd))
    pp = pp / sum(pp)

    pick = sample(k, 1, prob = pp)

    return (list("teamA"=players[teamA[pick,]], "ratingA" = round(scA[pick]),
        "teamB"=players[teamB[pick,]], "ratingB" = round(scB[pick])))
    }

random.team(c("mickey", "tony", "racer", "michael", "trevor", "austen"), card)
