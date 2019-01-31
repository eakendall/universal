# setwd("C:/Users/ekendal2/OneDrive - Johns Hopkins University/Research/universal regimen/universal")

date <- "20190109"
setting <- "SEA"

source("bpamz_cohort.R")


cohortsize <- 1e8 
# 1e7 is barely enough to capture all but the rarest combinations at least once. Really using this to get freqs, and would be easier to just calculate probs rather than creating a cohort but I've already done this so will go with it. Should use 1e8, but that req's >>4gb memory.

# make the cohort
cohort <- make.cohort(params=params, patientvars = patientvars, N=cohortsize)
c <- as.data.frame(table(cohort),stringsAsFactors = F)
c[rev(order(c$Freq)),]
c[] <- lapply(c, function(x) as.numeric(x))
distincts_with_freqs <- c[rev(order(c$Freq)),] # rather than modeling based on frequency in the population, model each patient type the same number of times and then bootstrap (if we increase the number of different types, may want to model the common types more than the rarer types)
c <- distincts_with_freqs %>% filter(MOXI==1|partialmoxi==0) # an imposible combo as coded
save(c, file=paste0("cohort.",date,".Rdata"))

sum(c$Freq)
tail(c %>% filter(INH==1|RIF==0)) # Rif monos essentially don't exist in India, so checking that other combos aren't missing
nrow(c %>% filter(Freq>0)) # number of patient types to be run, ~140 for India
c <- c %>% filter(Freq>0)
reps <- 5e3 # how many reps of each patient type to run? (will then sample Freq of each to recreate original cohort or subset thereof)
# 1e3 already gives a >1GB list for 3 scenarios, so may want to filter the cohort (e.g. to RR of FQ-R only) before the next step.  
# or, need to run sims and analyses on machine with more memory.

impact <- list()
impact$baseline <- modelcourse(scenario = "0", c, params, reps = reps)
impact$novelrr <- modelcourse(scenario = "1a", c, params, reps = reps)
impact$novelrrx <- modelcourse(scenario = "1x", c, params, reps = reps)
impact$novelpantb <- modelcourse(scenario = "3", c, params, reps = reps)
save(impact, file = paste0("impact.",date,".Rdata"))

moreimpact <- list()
moreimpact$all4 <- modelcourse(scenario = "2a", c, params, reps = reps)
moreimpact$all6 <- modelcourse(scenario = "2b", c, params, reps = reps)
save(moreimpact, file = paste0("moreimpact.",date,".Rdata"))
rm(moreimpact)

dst <- list()
dst$novelrrx <- impact$novelrrx
dst$noxxdr <- modelcourse(scenario = "3x", c, params, reps = reps)
dst$stepxxdr <- modelcourse(scenario = "4x", c, params, reps = reps)
dst$fullxxdr <- modelcourse(scenario = "5x", c, params, reps = reps)
save(dst, file = paste0("dst.",date,".Rdata"))
rm(dst)

# SA: delays Don't really want to advocate for getting rid of DST, and skipping xpert would have TB detection downsides that we aren't modeling. 
## ...so the main question of interest may be the elimination of delays associated with needing a separate DR regimen -- 
## ... i.e. what additional benefit from using the same drug combination for nearly all patients. For this, would want a regimen delay,
## ... and want to see how each of the scenarios is affected by it, so run the whole impact (and DST?) blocks again.
## ... Or maybe, just the impact blocks (although novelrr wouldn't change anything))
saveparams <- params
allparams <- read.csv("allparams.csv", header=T, stringsAsFactors = F)
# assuming SEA setting here
params <- as.numeric(allparams[,2]); names(params) <- allparams[,1]; params["DSTdelay"] <- 0
params["Tbdxtime_recurrence"] <- params["Tbdxtime_recurrenceratio"]*params["Tbdxtime"]
delays <- list()
delays$baseline <- modelcourse(scenario = "0", c, params, reps = reps)
delays$novelrr <- modelcourse(scenario = "1a", c, params, reps = reps)
delays$novelrrx <- modelcourse(scenario = "1x", c, params, reps = reps)
delays$novelpantb <- modelcourse(scenario = "3", c, params, reps = reps)
save(delays, file = paste0("delays.",date,".Rdata"))
params <- saveparams 
rm(delays)
