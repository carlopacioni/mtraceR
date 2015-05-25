# This small script has a few functions to test memory and speed of functions
# The general structure section is just a brief outline of what I set up to test
# a new function. I have (and assumed you would have) a bayesallfile one level
# up from the root of the package (i.e. "../") 
# In my case, mtrace works with default except for burn.in and npop. If you want 
# we can set them up programmatically (e.g. mtrace(burn.in = b, npop = np) and 
# have a line that source our own settings, something like:

# source(file = "../carlo.R")

# with carlo.R looking like this:

# b <- 0.25
# mp <- 2



library(lineprof)
library(pryr)

source(file = "./utilities/mem_lsos.R")
lsos()

###########################################################
# General structure
profFUN <- lineprof(test<-FUN, torture=F)
lsos()
shine(profFUN)


profsum<-(data.frame(prof="profFUN", time=sum(profFUN$time), 
                     alloc=sum(profFUN$alloc)))
profsum
time_mins <- profsum[, 2] / 60
time_mins

# in case case I do a second test with something different
profsum<-rbind(profsum, data.frame(prof="profFUN2", time=sum(profFUN2$time), 
                                   alloc=sum(profFUN2$alloc)))

mem_used()
tt<-system.time(mchg<-mem_change(test1<-FUN))
mem_used()
mchg
tt


##########################################################
# mtrace
profFUN <- lineprof(test<-mtrace(burn.in = 0.25, npop = 2, dir.in = "../"), 
                    torture=F)
lsos()
shine(profFUN)


profsum<-(data.frame(prof="profFUN", time=sum(profFUN$time), 
                     alloc=sum(profFUN$alloc)))
profsum

time_mins <- profsum[, 2] / 60
time_mins


profFUN2 <- lineprof(test<-mtrace(burn.in = 0.25, npop = 2, thin=10, dir.in = "../"), 
                    torture=F)
lsos()
shine(profFUN2)

profsum<-rbind(profsum, data.frame(prof="profFUN2", time=sum(profFUN2$time), 
                                   alloc=sum(profFUN2$alloc)))
profsum

time_mins <- profsum[2, 2] / 60
time_mins




mem_used()
tt<-system.time(mchg<-mem_change(test1<-mtrace(burn.in = 0.25, npop = 2,
                                               dir.in = "../")))
mem_used()
mchg
tt
