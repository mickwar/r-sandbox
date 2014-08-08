# mw: written by Arthur Lui 2014, used to send mass emails,
# specifically home teaching assignments to an elders quorum
#Read and Organize Data: #######################################################
chr <- function(n) rawToChar(as.raw(n)) 

ht <- read.csv("HT.csv") # File containing the full names of HT-er and Teachee pairs
cn <- read.csv("contact.csv") # contact list with email and full names
cn[,1] <- as.character(cn[,1])
cn[,2] <- as.character(cn[,2])

teacher <- as.character(ht[,1])
student <- as.character(ht[,2])

label <- function(x) {
  lab <- 65
  x.lab <- NULL
  for (i in 1:(length(x)-1)) {
    if (x[i] == "" & x[i+1] != "") {
      lab <- lab + 1
    }  
    x.lab[i] <- ifelse(x[i]=="","",chr(lab))
  }
  x.lab[i+1] <- ifelse(x[i+1]=="","",chr(lab))
  x.lab
}

teacher.lab <- label(teacher)
student.lab <- label(student)

labs <- unique(teacher.lab)  
labs <- labs[labs != ""]

teacher.list <- lapply(as.list(labs),function(x) teacher[which(teacher.lab==x)])
student.list <- lapply(as.list(labs),function(x) student[which(student.lab==x)])

#Write Emails: #################################################################

write.email <- function(comps,studs,samp="sampEmail.txt",out="") {
  txt <- readLines(samp)
  n <- length(txt)
  comp.line <- which(txt=="<comps>")
  txt <- c(txt[1:comp.line-1],paste("-",comps),txt[(comp.line+1):n]) 
  n <- length(txt)
  stud.line <- which(txt=="<teachees>")
  txt <- c(txt[1:stud.line-1],paste("-",studs),txt[(stud.line+1):n]) 
  if (out!="") writeLines(txt,out)
  txt
}

emails <- lapply(as.list(1:length(labs)),function(x) write.email(teacher.list[[x]],student.list[[x]],out=paste("emails/email",x,".txt",sep="")))



send.mail <- function() {

  # sub functions:
    send <- function(email,address) {
      cmd <- paste("mail -s \"Home Teaching\"", address, "<", email)
      system(cmd)
    }

    mail.to.teachers <- function(gp) { # input is the group number
      tchrs.email <- cn[,2][which(cn[,1] %in% teacher.list[[gp]])]
      n <- length(tchrs.email)
      file <- paste("emails/email",gp,".txt",sep="")

      for (i in 1:n) {
        send(file,tchrs.email[i])
      }
    }
  
  # main:
    for (i in 1:length(teacher.list)) {
      mail.to.teachers(i)
    }
}


#send.mail()
