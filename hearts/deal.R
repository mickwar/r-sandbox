# cards = c(  "02c", "03c", "04c", "05c", "06c", "07c", "08c",
#             "09c", "10c", "11c", "12c", "13c", "14c",
#             "02d", "03d", "04d", "05d", "06d", "07d", "08d",
#             "09d", "10d", "11d", "12d", "13d", "14d",
#             "02s", "03s", "04s", "05s", "06s", "07s", "08s",
#             "09s", "10s", "11s", "12s", "13s", "14s",
#             "02h", "03h", "04h", "05h", "06h", "07h", "08h",
#             "09h", "10h", "11h", "12h", "13h", "14h")
# 11 = Jack
# 12 = Queen
# 13 = King
# 14 = Ace

#cards = c(  "02c", "03c", "04c", "05c", "06c", "07c", "08c",
#            "09c", "10c", " Jc", " Qc", " Kc", " Ac",
#            "02d", "03d", "04d", "05d", "06d", "07d", "08d",
#            "09d", "10d", " Jd", " Qd", " Kd", " Ad",
#            "02s", "03s", "04s", "05s", "06s", "07s", "08s",
#            "09s", "10s", " Js", " Qs", " Ks", " As",
#            "02h", "03h", "04h", "05h", "06h", "07h", "08h",
#            "09h", "10h", " Jh", " Qh", " Kh", " Ah")

cards = c(  "02a", "03a", "04a", "05a", "06a", "07a", "08a",
            "09a", "10a", "11a", "12a", "13a", "14a",
            "02b", "03b", "04b", "05b", "06b", "07b", "08b",
            "09b", "10b", "11b", "12b", "13b", "14b",
            "02c", "03c", "04c", "05c", "06c", "07c", "08c",
            "09c", "10c", "11c", "12c", "13c", "14c",
            "02d", "03d", "04d", "05d", "06d", "07d", "08d",
            "09d", "10d", "11d", "12d", "13d", "14d")


deal = function(cards){
    H = matrix(sample(cards), 4, 13)
    # order: club, diamond, spade, heart, increasing
    for (i in 1:4){
        H[i,] = sort(H[i,])
        H[i,] = H[i,][order(substring(H[i,], 3, 3))]
        suits = substring(H[i,], 3, 3)
        suits = ifelse(suits == "d", "h", suits)
        suits = ifelse(suits == "c", "s", suits)
        suits = ifelse(suits == "b", "d", suits)
        suits = ifelse(suits == "a", "c", suits)
        substring(H[i,], 3, 3) = suits
        }
    return(H)
    }

write.table(deal(cards), "./hands.txt", quote = FALSE,
    col.names = FALSE, row.names = FALSE)

