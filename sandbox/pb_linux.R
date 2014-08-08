get.width = function()
    min(100, as.numeric(system("tput cols", intern=TRUE)))
terminal.width = get.width()

pb.linux = function(current, max.iter){
    if (current == 1)
        cat(rep("#", terminal.width), "\n", sep="")
    ratio = max.iter / terminal.width
    if (floor(current %% max(1, ratio)) == 0){
        cat(rep("*", max(1, floor(1/ratio))), sep="")
        if ((floor(1/ratio*current)-
            floor(1/ratio)*current)-
            (floor(1/ratio*(current-1))-
            floor(1/ratio)*(current-1)) && ratio < 1) 
            cat("*", sep="")
        }
    if (current == max.iter)
        cat("\n")
    }

# when the iterations do not occur in order (such as for
# parallel processing)
pb.linux.nonseq.start = function()
    cat(rep("#", get.width()), "\n", sep="")

pb.linux.nonseq = function(at, current, max.iter){
    ratio = max.iter / get.width()
    if (at > current){
        current = at
        if (floor(current %% max(1, ratio)) == 0){
            cat(rep("*", max(1, floor(1/ratio))), sep="")
            if ((floor(1/ratio*current)-
                floor(1/ratio)*current)-
                (floor(1/ratio*(current-1))-
                floor(1/ratio)*(current-1)) && ratio < 1) 
                cat("*", sep="")
            }
        }
    if (current == max.iter)
        cat("\n")
    return (current)
    }
