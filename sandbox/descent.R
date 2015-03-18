roll.attack.dice = function(nblue = 1, nred = 0, nyellow = 0){
    hearts = 0
    range = 0
    surge = 0
    blue = list(c(-Inf, -Inf, -Inf), c(2, 2, 1), c(2, 4, 0),
        c(1, 6, 1), c(2, 3, 0), c(1, 5, 0))
    red = list(c(3, 0, 1), c(1, 0, 0), c(2, 0, 0),
        c(2, 0, 0), c(3, 0, 0), c(2, 0, 0))
    yellow = list(c(2, 0, 1), c(2, 0, 0), c(1, 1, 0),
        c(0, 1, 1), c(1, 0, 1), c(1, 2, 0))
    if (nblue >= 1L){
        for (i in 1L:nblue){
            x = sample(blue, 1L, replace = TRUE)
            hearts = hearts + x[[1L]][1L]
            range = range + x[[1L]][2L]
            surge = surge + x[[1L]][3L]
            }
        }
    if (nred >= 1L){
        for (i in 1L:nred){
            x = sample(red, 1L, replace = TRUE)
            hearts = hearts + x[[1L]][1L]
            range = range + x[[1L]][2L]
            surge = surge + x[[1L]][3L]
            }
        }
    if (nyellow >= 1L){
        for (i in 1L:nyellow){
            x = sample(yellow, 1L, replace = TRUE)
            hearts = hearts + x[[1L]][1L]
            range = range + x[[1L]][2L]
            surge = surge + x[[1L]][3L]
            }
        }
    if (hearts == -Inf)
        hearts <- range <- surge <- 0
    return (c(hearts, range, surge))
    }

roll.defend.dice = function(nblack, ngray, nbrown){
    shield = 0
    black = c(4, 0, 3, 2, 2, 2)
    gray = c(3, 2, 1, 1, 1, 0)
    brown = c(2, 1, 1, 0, 0, 0)
    if (nblack >= 1)
        shield = shield + sum(sample(black, nblack, replace = TRUE))
    if (ngray >= 1)
        shield = shield + sum(sample(gray, ngray, replace = TRUE))
    if (nbrown >= 1)
        shield = shield + sum(sample(brown, nbrown, replace = TRUE))
    return (shield)
    }


nrolls = 100000
attack.roll = matrix(0, nrolls, 3)
defend.roll = double(nrolls)
for (i in 1:nrolls){
    attack.roll[i,] = roll.attack.dice(1, 0, 1)
    defend.roll[i] = roll.defend.dice(0, 1, 0)
    }

mean(attack.roll[,1] > defend.roll)

