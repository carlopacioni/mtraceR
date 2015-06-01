# This small script has a few functions to test memory and speed of functions
# The general structure section is just a brief outline of what I set up to test
# a new function. I have (and assumed you would have) a bayesallfile one level
# up from the root of the package (i.e. "../") 



library(mtraceR)
library(lineprof)
library(pryr)

source(file = "./utilities/mem_lsos.R")
lsos()

#------------------------------------------------------------------------------#
# mtrace()

# link to large bayesallfile.gz (~600MB)
# https://www.dropbox.com/s/w7zoe714myp3iro/large_bayesallfile.gz?dl=0

# link to medium bayesallfile.gz (~300MB)
# https://www.dropbox.com/s/4a6mja85502b632/medium_bayesallfile.gz?dl=0

loc.path <- "C:/Users/30373314/Dropbox/mtraceR/TestFiles"
# loc.path <- "put here your path to test files"
large <- "large_bayesallfile.gz"
medium <- "medium_bayesallfile.gz"

# Assuming you are in the root directory of the package, this creates a folder 
# one level up where to store resutls 
dir.create("../tests", showWarnings=FALSE)
out <- "../tests"

profFUN <- lineprof(test<-mtrace(burn.in=0.25, dir.in=loc.path, dir.out=out,
                                 bayesallfile=medium), interval=0.1, torture=F)
lsos()
shine(profFUN)
profFUN

profsum<-(data.frame(prof="profFUN", time=sum(profFUN$time), 
                     alloc=sum(profFUN$alloc)))
profsum

time_mins <- profsum[, 2] / 60
time_mins


profFUN2 <- lineprof(test<-mtrace(burn.in=0.25, dir.in=loc.path, dir.out=out,
                                  bayesallfile=large), interval=0.1, torture=F)
lsos()
shine(profFUN2)

profsum<-rbind(profsum, data.frame(prof="profFUN2", time=sum(profFUN2$time), 
                                   alloc=sum(profFUN2$alloc)))
profsum

time_mins <- profsum[2, 2] / 60
time_mins

write.csv(profsum, file=paste0("./utilities/", "profsum", Sys.Date(), ".csv"))
write.csv(print(profFUN), file=paste0("./utilities/", "profmedium", Sys.Date(), ".csv"))
write.csv(print(profFUN2), file=paste0("./utilities/", "proflarge", Sys.Date(), ".csv"))

mem_used()
tt<-system.time(mchg<-mem_change(test1<-mtrace(burn.in=0.25, dir.in=loc.path, 
                                               dir.out=out,
                                               bayesallfile=medium)))
mem_used()
mchg
tt
