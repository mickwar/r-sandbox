attack = function(natt, ndef){
    if (ndef > natt)
        ndef = natt # cannot defend with more than the number of attackers
    out = 0
    # roll the two sets of dice
    attack_roll = sort(sample(6, natt, replace = TRUE), decreasing = TRUE)
    defend_roll = sort(sample(6, ndef, replace = TRUE), decreasing = TRUE)
    # count number of defenders that are killed
    for (i in 1:ndef)
        out = out + (attack_roll[i] > defend_roll[i])*1
    return(out)
    }

sim_battle = function(nattacker, ndefender, natt, ndef){
    nrolls = 0
    while (nattacker > 0 && ndefender > 0){
        nrolls = nrolls + 1
        x = attack(min(natt, nattacker), min(ndefender, ndef))
        nmax = min(natt, nattacker, ndef, ndefender)
        ndefender = ndefender - x
        nattacker = nattacker - (nmax - x)
        }
    return((c("natt"=nattacker, "ndef"=ndefender, "nrolls"=nrolls)))
    }

# run 10000 battles where you begin attacking with 8 and the
# opponent begins defending with 5, using the highest number
# of attacks and defenders possible each roll
x = matrix(0, 10000, 3)
for (i in 1:nrow(x))
    x[i,] = sim_battle(8, 5, 3, 2)

# estimated probabilities of the outcomes
table(x[,1]) / nrow(x)
table(x[,2]) / nrow(x)

# plots of the outcomes
par(mfrow=c(2,1))
plot(table(x[,1]) / nrow(x))
plot(table(x[,2]) / nrow(x))




####

mean(attack(1, 1, 100000)) # approximates 15/(6 * 6)
mean(attack(2, 1, 100000)) # approximates 125/(36 * 6)
mean(attack(3, 1, 100000)) # approximates 855/(216 * 6)

x11 = attack(1, 1, 100000)
x21 = attack(2, 1, 100000)
x31 = attack(3, 1, 100000)
x22 = attack(2, 2, 100000)
x32 = attack(3, 2, 100000)

table(x11) / length(x11)
table(x21) / length(x21)
table(x31) / length(x31)
table(x22) / length(x22)
table(x32) / length(x32)


temp = double(10000)
for (i in 1:length(temp))
    temp[i] = sum(attack(2, 1, 2))

table(temp)/length(temp)
table(x22) / length(x22)



temp2 = double(10000)
for (i in 1:length(temp2))
    temp2[i] = sum(attack(3, 1, 2))

table(temp2)/length(temp2)
table(x32) / length(x32)








### true probabilities (with fair dice)
# probability that attacker wins on 1 vs 1
15/(6 * 6)

# probability that attacker wins on 2 vs 1
125/(36 * 6)

# probability that attacker wins on 3 vs 1
855/(216 * 6)
