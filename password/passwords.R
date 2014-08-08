dict.path = "/usr/share/dict/"
file.name = "cracklib-small"

dictionary = read.table(paste0(dict.path, 
    file.name))[,1]
alternate = read.table("./alternatives.txt")

# can I find words that are in the dictionary
# even if some types something like "P@S$w0Rd"?
# or better, can i find the dictionary word
# in "AB1P@S$w0RdXY2"?
word = "P@S$w0Rd"
check.string = function(word){
    match.found = FALSE
    while (!match.found){

        }
    }

grep("zyg", dictionary)
