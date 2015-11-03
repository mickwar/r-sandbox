# Read in dataset
dat = read.csv("~/files/data/mustard_gas/data.csv")

# Remove the columns with only NAs
dat = dat[,-which(apply(dat, 2, function(x) all(is.na(x))))]

# Remove the rows which are all blank
dat = dat[-which(apply(dat, 1, function(x) all(x == ""))),]

# Get military branch
branch = as.character(dat$branch)

# Get birth and death years
birth = unlist(lapply(strsplit(as.character(dat$dob), "/"), function(x) ifelse(length(x)==0, "", x[length(x)])))
birth = as.numeric(ifelse(nchar(birth) == 4, birth, NA))

death = unlist(lapply(strsplit(as.character(dat$dod), "/"), function(x) ifelse(length(x)==0, "", x[length(x)])))
death = as.numeric(ifelse(nchar(death) == 4, death, NA))

# Remove weird birth years
rm = which(birth < 1895 | birth > 1932)

birth = birth[-rm]
death = death[-rm]
branch = branch[-rm]

# Remove those with no birth year
rm = which(is.na(birth))

birth = birth[-rm]
death = death[-rm]
branch = branch[-rm]

#cbind(is.na(birth), is.na(death))

# Remove branch that isn't Army or Navy (other categories have very little in them)
rm = which(!(branch == "Army" | branch == "Navy"))

birth = birth[-rm]
death = death[-rm]
branch = branch[-rm]

age = death - birth

hist(birth, col = 'gray')
hist(death, col = 'gray')

hist(age, col = 'gray')


