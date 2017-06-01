dat = read.table("./basketball_scores.txt", sep = ",", header = FALSE)
n = NROW(dat)

### Need to count number of games played

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
scores[1,] = rep(800, length(players))
colnames(scores) = players

for (i in 1:n){
    # Carry over the scores
    scores[i+1, ] = scores[i, ]
    scores[i+1, ] = scores[i, ]

    # Carry over number of games
    ngames[i+1, ] = ngames[i, ]
    ngames[i+1, ] = ngames[i, ]

    # Get index for winners and losers
    A.ind = which(players %in% team.A[[i]])
    B.ind = which(players %in% team.B[[i]])

    # Use sum of player scores to determine team score
    R_A = mean(scores[i, A.ind])
    R_B = mean(scores[i, B.ind])

    # Expected score for the winners
    E_A = 1 / (1 + exp((R_B - R_A) / 40))

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

    # Divide a player's K-factor by 2 when that player has played more than 30 games
    K.fac[A.ind] = K.fac[A.ind] / ifelse(ngames[i, A.ind] > 30, 2, 1)
    K.fac[B.ind] = K.fac[B.ind] / ifelse(ngames[i, B.ind] > 30, 2, 1)

    # Calculate observed as using the score difference (max 10 point diff)
#   S_A = (min(win.score[i] - loss.score[i], 10) + 10) / 20
#   S_B = 1 - S_A
    S_A = 1 / (1 + exp((B.score[i] - A.score[i]) / 2))
    S_B = 1 - S_A

    # Adjust player scores
    # Team A won or drawed (A.score >= B.score), but not by enough (S_A <= E_A),
    # but still don't decrease the score
    if ( (A.score[i] >= B.score[i]) && (S_A <= E_A) )
        K.fac[A.ind] = 0
    # Same for B
    if ( (B.score[i] >= A.score[i]) && (S_B <= E_B) )
        K.fac[B.ind] = 0

    scores[i+1, A.ind] = scores[i, A.ind] + K.fac[A.ind]*(S_A - E_A)
    scores[i+1, B.ind] = scores[i, B.ind] + K.fac[B.ind]*(S_B - E_B)

    # No player's score can drop below 600
    scores[i+1, scores[i+1,] < 600] = 600

    ngames[i+1, A.ind] = ngames[i, A.ind] + 1
    ngames[i+1, B.ind] = ngames[i, B.ind] + 1

    }

# Official scores
round(scores)

tail(round(scores), 1)

# The average score across all players over time (would like to see this increase)
apply(scores, 1, mean)

# Change in player scores over time
apply(scores, 2, diff)

# For every 20 rank difference, a win by 1 point is considered a draw
(1 / (1 + exp(-1 / 2)))
(1 / (1 + exp(-20 / 40)))
