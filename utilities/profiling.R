# This small script has a few function to test memory and speed of functions
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

# I'm not using this for now, left it here just  in case
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


profsum<-rbind(profsum, data.frame(prof="profFUN2", time=sum(profFUN2$time), 
                                   alloc=sum(profFUN2$alloc)))

mem_used()
tt<-system.time(mchg<-mem_change(test1<-mtrace(burn.in = 0.25, npop = 2,
                                               dir.in = "../")))
mem_used()
mchg
tt
