dat = read.table("./basketball_scores.txt", sep = ",", header = FALSE)
n = NROW(dat)

# Make a function that takes in team rating and number of games played
# by each player and compute the new ratings. Use inside the loop.
# Can also be used to calculate ratings for hypothetical games.


# Need to standardize game length (7 (6) or 14 (12) minutes only)
# Or just use the observed score as Score / Time (still need time)
# Or consider something like percent of points won by: 10 to 11 is different than 100 to 101

# Starting rating for all players
starting_rating = 1200

# minimum score
min_score = 600

# like it sounds (involved in calculation of observed and expected scores)
logistic_scale = 5

# For every point_diff rating points that team A is higher
# than team B, team A must score 1 more point to increase
# their rating (involved in calculation of expected score only)
point_diff = 40

# K-factor groups
K.group = c(600, 1000, 1400, 1800, 2200, 2600, Inf)
# K.val = c(120, 96, 72, 48, 24)*2
# K.val = c(120, 96, 72, 48, 24, 12)*2
K.val = c(6, 5, 4, 3, 2, 1) * 40

c(120, 96, 72, 48, 24)*2
c(120, 96, 72, 48, 24, 12)*2
c(6, 5, 4, 3, 2, 1) * 40

# Multiply K.factor by K.game(ngames) --- First few games should have volatile
# changes to ratin, but the impact grows less and elss over time
K.game = function(n){
    out = double(length(n))
    for (i in 1:length(n)){
        if (n[i] <= 10)
            out[i] = 4
        if (n[i] > 10 && n[i] <= 20)
            out[i] = 2
        if (n[i] > 20 && n[i] <= 30)
            out[i] = 1
        if (n[i] > 30)
            out[i] = 1/2
        }
    return (out)
    }


# K.val represents the largest possible change, within a certain
# tier, between two evenly matched opponents
# The smallest change (a 1-point victory) is roughly K.val/20
# These are assuming the first K.game[1] games have already been played

# Function for changing K.fac multiplayer in a team
team_func = function(n) 1/sqrt(n)




team.A = strsplit(as.character(dat[,1]), " ")
team.B = strsplit(as.character(dat[,2]), " ")

A.score = as.numeric(sapply(team.A, function(x) tail(x, 1)))
B.score = as.numeric(sapply(team.B, function(x) tail(x, 1)))

team.A = lapply(team.A, function(x) x[-length(x)])
team.B = lapply(team.B, function(x) x[-length(x)])


# Get players
players = unique(c(unlist(lapply(strsplit(as.character(dat[,1]), " "), function(x) x[-length(x)])),
    unlist(lapply(strsplit(as.character(dat[,2]), " "), function(x) x[-length(x)]))))
players = sort(players)



scores = matrix(0, n + 1, length(players))
ngames = matrix(0, n + 1, length(players))
nwin = matrix(0, n + 1, length(players))
nloss = matrix(0, n + 1, length(players))
scores[1,] = rep(starting_rating, length(players))
# scores[1,] = ending
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

    # Use mean of player scores to determine team score
    R_A = mean(scores[i, A.ind])
    R_B = mean(scores[i, B.ind])

    # Expected score for the winners
    E_A = 1 / (1 + exp((R_B - R_A) / (logistic_scale*point_diff)))

    # Expected score for the losers
    E_B = 1 - E_A
    
    # K-factor (adjusted by number of players on the team)
    K.fac = double(length(players))
#   K.fac[c(A.ind, B.ind)] = 48
    for (j in A.ind){
        for (k in 1:length(K.val))
            if (scores[i, j] >= K.group[k] & scores[i, j] < K.group[k+1])
                K.fac[j] = K.val[k]
        }
    for (j in B.ind){
        for (k in 1:length(K.val))
            if (scores[i, j] >= K.group[k] & scores[i, j] < K.group[k+1])
                K.fac[j] = K.val[k]
        }

    # Reduce K-factor by size of the team
    K.fac[A.ind] = K.fac[A.ind] * team_func(length(A.ind))
    K.fac[B.ind] = K.fac[B.ind] * team_func(length(B.ind))

    # Adjust K-factor by number of games played
    K.fac[A.ind] = K.fac[A.ind] * K.game(ngames[i, A.ind])
    K.fac[B.ind] = K.fac[B.ind] * K.game(ngames[i, B.ind])

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

# around = c("austen", "lai", "mickey", "racer", "tony", "trevor")
# ending = tail(scores, 1)
# ind = which(colnames(scores) %in% around)
# ind = 1:8
# scores = scores[,ind]
# ngames = ngames[,ind]
# nwin = nwin[,ind]
# nloss = nloss[,ind]

# Official scores
round(scores)

tail(round(scores), 10)
tail(scores, 8)

final_score = round(tail(scores, 1))
final_ngames = tail(ngames, 1)
final_nwin = tail(nwin, 1)
final_nloss = tail(nloss, 1)


card = rbind(final_score, final_ngames,
    final_nwin, final_nloss)

rownames(card) = c("ELO", "GP", "W", "L")
card = data.frame(card)

# The average score across all players over time (would like to see this increase)
apply(scores, 1, mean)

# Change in player scores over time
apply(scores, 2, diff)

# Player scores over time
matplot(apply(scores, 2, diff))

mm = apply(scores, 1, mean)
dd = diff(apply(scores, 1, mean))
cc = rep("black", n)
cc[which(dd > 0)] = "green"
cc[which(dd < 0)] = "red"

scores[scores == 1200] = NA
matplot(scores)
#lines(apply(scores, 1, mean))  # Average of all players over time
segments(x0 = 1:n, x1 = 2:(n+1), y0 = mm[-(n+1)], y1 = mm[-1], col = cc)
legend("bottomleft", legend = colnames(scores),
    pch = c(as.character(c(1:9, 0)), letters, LETTERS), bty = 'n', col = (1:ncol(scores)-1) %% 6 + 1)

# Pick a random team based on scores
random.team = function(players, card, lowest = FALSE){
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
    pp = 1 - (dd - min(dd)) / diff(range(dd)) + 0.5
    pp = pp / sum(pp)

    pick = ifelse(lowest, which.max(pp), sample(k, 1, prob = pp))

    return (list("teamA"=players[teamA[pick,]], "ratingA" = round(scA[pick]),
        "teamB"=players[teamB[pick,]], "ratingB" = round(scB[pick]),
        "diff" = round(abs(scA[pick] - scB[pick]))))
    }

team = c("mickey", "tony", "racer", "trevor", "seth", "lai")
random.team(team, card, FALSE)

