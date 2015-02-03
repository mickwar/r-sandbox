library(RCurl)

olive = read.table(text = getURL(paste0("https://raw.githubusercontent.com/mickwar/",
    "data/master/666/oliver3b.txt")), header = TRUE)

dictionary = read.table(text = getURL(paste0("https://raw.githubusercontent.com/luiarthur/",
    "Boggle/master/dictionary.txt")))

# https://raw.githubusercontent.com/USER/REPOSITORY/BRANCH/FILE
