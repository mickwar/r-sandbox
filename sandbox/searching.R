# grep in the command line is better than this
search.R=function(pattern,
    start.dir="C:/Users/Mickey/Desktop/Downloads/School/"){

    setwd(start.dir)

    files=list.files(getwd(),recursive=TRUE)
    
    # keep .r only files
    new=NULL
    for (i in 1:length(files)){
        if (tolower(substring(files[i],nchar(files[i])-1,
            nchar(files[i])))==".r")
            new[length(new)+1]=files[i]
        }

    # deparse the scripts and check for the pattern
    out=list(NULL)
    for (i in 1:length(new)){
        temp=grep(pattern,readLines(new[i],warn=FALSE))
        if (length(temp)>0){
            if (is.null(out[[1]]))
                out[[1]]=list(new[i],temp)
            else
                out[[length(out)+1]]=list(new[i],temp)
            }
        }

    return(out)
    }

(h=search.R("plot", "~/files/R"))
