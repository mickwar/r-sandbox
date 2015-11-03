### Simulate monopoly moving

squares = c(
    "Go", "Mediterranean", "Community1", "Baltic", "Income", "RRReading",
    "Oriental", "Chance1", "Vermont", "Connecticut", "Jail", "StCharles",
    "Electrc", "States", "Virginia", "RRPenn", "StJames", "Community2",
    "Tennessee", "NewYork", "FreeParking", "Kentucky", "Chance2", "Indiana",
    "Illinois", "RRB&O", "Atlantic", "Ventnor", "Plumbing", "Marvin",
#   "GotoJail", # Omitted since you can't end a move on this piece
    "Pacific", "NorthCarolina", "Community3", "Pennsylvania", "RRShortline",
    "Chance3", "ParkPlace", "Luxury", "Boardwalk")

move = function()
    sample(6, 2, replace = TRUE)

chance = function(loc){
    random.chance = sample(c("3spaces", "jail", "boardwalk", "reading",
        "go", "illinois", "stcharles", "utility", "railroad", "other"),
        1, prob = c(1,1,1,1,1,1,1,1,1,15-9))
    switch(random.chance,
          "3spaces" = ifelse(loc - 3 >= 0, loc - 3, 40 + loc - 3),
            "other" = loc,
               "go" = 0,
          "reading" = 5,
             "jail" = 10,
        "stcharles" = 11,
         "illinois" = 24,
        "boardwalk" = 39,
          "utility" = ifelse(loc == 22, 28, 12),
         "railroad" = ifelse(loc == 7, 15, ifelse(loc == 22, 25, 5)))
    }

community = function(loc){
    random.chest = sample(c("jail", "go", "other"),
        1, prob = c(1,1,16-2))
    switch(random.chest,
            "other" = loc,
               "go" = 0,
             "jail" = 10)
    }

# n = Number of times to roll the dice (i.e. to land on a square)
get.turns = function(n){
    double.jail = 0
    doubles = 0
    turns = 0   # Number of actual turns (accounting for doubles)
    locs = double(n)
    for (i in 2:n){
        x = move()
        doubles = ifelse(x[1] == x[2], doubles + 1, 0)

        locs[i] = locs[i-1] + sum(x)

        # Count number of turns
        if (doubles == 0)
            turns = turns + 1

        # Rolled three doubles in a row
        if (doubles == 3){
            doubles <<- 0
            turns = turns + 1
            double.jail = double.jail + 1
            locs[i] = 10
            }

        # Passed Go
        if (locs[i] >= 40)
            locs[i] = locs[i] - 40

        # Landed on a Chance
        if ((locs[i] == 7) || (locs[i] == 22) || (locs[i] == 36))
            locs[i] = chance(locs[i])

        # Landed on a Community Chest
        if ((locs[i] == 2) || (locs[i] == 17) || (locs[i] == 33))
            locs[i] = community(locs[i])

        # Landed on Jail
        if (locs[i] == 30)
            locs[i] = 10
        }

    # The proportion of times a player lands on each location (as n goes to infinity)
    y = table(locs) / n
    landed = as.numeric(names(table(locs)))
    landed[landed < 30] = landed[landed < 30] + 1

    # Expand the vector so of length 39
    out = double(39)
    out[(1:39) %in% landed] = y


    names(out) = squares
    return (list("y" = out, "jail" = double.jail, "turns" = turns))
    }

B = 1e7     # Number of loops
n = 100     # Number if dice rolls

landings = matrix(0, B, 39)
num.doublejails = double(B)
num.turns = double(B)
for (i in 1:B){
    x = get.turns(n)
    landings[i,] = x$y
    num.doublejails[i] = x$jail
    num.turns[i] = x$turns
    }
b = landings
landings = head(landings, i-1)
num.doublejails = head(num.doublejails, i-1)
num.turns = head(num.turns, i-1)
colnames(landings) = squares

par(mfrow = c(2,1), mar = c(5.1, 4.1, 1.1, 1.1))
plot(table(num.doublejails) / B, type='h', ylab = "Proportion of runs",
    xlab = "Number of three double jails")
plot(table(num.turns) / B, type='h', ylab = "Proportion of runs",   
    xlab = "Number of total turns")
points(0:n, dbinom(0:n, n, mean(num.turns)/n), col = 'red', pch = 20, cex = 0.5)

mean(num.turns) / n


par(mfrow = c(2,2), mar = c(5.1, 4.1, 1.1, 1.1))
for (i in 1:39){
    plot(table(landings[,i]) / B, type = 'h', xlab = squares[i])
    abline(v = mean(landings[,i]), col = 'red')
    x = as.numeric(names(table(landings[,i])))
    xx = seq(min(x), max(x), by = 1/n)
#   points(xx, dpois(0:(length(xx)-1), num.turns[i]*mean(landings[,i])), col = 'red', pch = 20)
    readline()
    }

cbind(apply(landings, 2, mean), apply(landings, 2, sd))


# # Try as a discrete-space Markov chain?
# A = matrix(0, 40, 40)
# A[1,] = 
#c(1/36 * 1/16 + 6/36 * 1/16, 0, 1/36 * 14/16, 3/36, 4/36, 5/36, 6/36 * 
# for (i in 1:40){
#     (i+2):(i+12)
# c(1:6, 5:1) / 36

hist(y, breaks = 20, col = 'gray')
plot(locs[1:100], type='l')
points(locs[1:100], pch = 20, cex = 0.5)

plot(y, type='h')
par(mar = c(8.1, 4.1, 4.1, 2.1))
plot(y, type='h', axes = FALSE)
axis(side = 1, at = 1:39, labels = squares, las = 2)
axis(2)
