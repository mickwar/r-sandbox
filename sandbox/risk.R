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
x = matrix(0, 100000, 3)
for (i in 1:nrow(x))
    x[i,] = sim_battle(8, 5, 3, 2)

# estimated probabilities of the outcomes
p = table(x[,1]) / nrow(x)

# plots of the outcomes (for the attacker)
plot(p, ylab="", xlab="", ylim=c(0,1), cex.lab=2, cex.main = 2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=1.5)


#### battle rolls
x = matrix(0, 100000, 3) 
p = NULL

for (d in 1:8){
    cat("\rIter: ", d, "/", 8)
    for (i in d:(d+7)){
        for (j in 1:nrow(x)) x[j,] = sim_battle(i, d, 3, 2)
        p[[i-d+1]] = table(x[,1]) / nrow(x)
        }
    pdf(paste0("risk_d",d,".pdf"), height = 25, width = 20)
    par(mfrow=c(4,2), mar=c(4.1,4.1,2.1,2.1), oma=c(1,3,5,1))
    for (i in 1:8){
        plot(p[[i]], ylab="", xlab="", ylim=c(0,1), xlim=c(-0.05, length(p[[i]])-0.95))
#       text(0:(length(p[[i]])-1), p[[i]]+0.05, round(p[[i]], 3), cex=2.5)
        text(0:(length(p[[i]])-1), p[[i]]+0.05, round(p[[i]]*1000), cex=2.5)
        }
    mtext(paste0("Defender: ", d), 3, outer = TRUE, cex = 4)
    mtext("x-axis: remaining attackers", 1, outer = TRUE)
    dev.off()
    }

#### single roll probabilities
x11 <- x21 <- x31 <- x22 <- x32 <- double(1000000)
for (i in 1:length(x11)) x11[i] = attack(1, 1)
for (i in 1:length(x21)) x21[i] = attack(2, 1)
for (i in 1:length(x31)) x31[i] = attack(3, 1)
for (i in 1:length(x22)) x22[i] = attack(2, 2)
for (i in 1:length(x32)) x32[i] = attack(3, 2)

table(x11) / length(x11)
table(x21) / length(x21)
table(x31) / length(x31)
table(x22) / length(x22)
table(x32) / length(x32)

pdf("risk_roll.pdf", height = 25, width = 20)
par(mfrow=c(3,2), mar=c(4.1,5.1,2.1,2.1), oma=c(3,5,5,1))
p = table(x11)/ length(x11)
plot(p, ylab = 1, xlab="", main = 1, ylim=c(0,1), xlim=c(-0.05,1.05), cex.lab = 2, cex.main=2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=2.5)
plot(0, type="n", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
p = table(x21) / length(x21)
plot(p, ylab=2, xlab="", ylim=c(0,1), xlim=c(-0.05, 1.05), cex.lab=2, cex.main = 2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=2.5)
p = table(x22) / length(x22)
plot(p, ylab = "", xlab="", main = 2, ylim=c(0,1), xlim=c(-0.05, 2.05), cex.lab=2, cex.main = 2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=2.5)
p = table(x31) / length(x31)
plot(p, ylab=3, xlab="", ylim=c(0,1), xlim=c(-0.05, 1.05), cex.lab=2, cex.main = 2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=2.5)
p = table(x32) / length(x32)
plot(p, ylab="", xlab="", ylim=c(0,1), xlim=c(-0.05, 2.05), cex.lab=2, cex.main = 2)
text(0:(length(p)-1), p+0.05, round(p, 3), cex=2.5)
mtext("Number of Attackers", 2, outer = TRUE, cex = 3)
mtext("Number of Defenders", 3, outer = TRUE, cex = 3)
dev.off()




### true probabilities (with fair dice)
# probability that attacker wins on 1 vs 1
15/(6 * 6)

# probability that attacker wins on 2 vs 1
125/(36 * 6)

# probability that attacker wins on 3 vs 1
855/(216 * 6)

# probabilities for outcomes on 2 vs 2
# Pr(X1 > Y1 and X2 > Y2)
505/1296

# Pr(X1 > Y1 and Y2 > X2)

# Pr(Y1 > X1 and X2 > Y2)

# Pr(Y1 > X1 and Y2 > X2)
