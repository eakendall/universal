# setwd("C:/Users/ekendal2/OneDrive - Johns Hopkins University/Research/universal regimen/universal")

date <- "20190226"
setting <- "SEA"
source("bpamz_cohort.R")

# cohortsize <- 1e8
# cohort <- make.cohort(params=params, patientvars = patientvars, N=cohortsize)
# c <- as.data.frame(table(cohort),stringsAsFactors = F)
# rm(cohort)
# c[rev(order(c$Freq)),]
# c[] <- lapply(c, function(x) as.numeric(x))
# distincts_with_freqs <- c[rev(order(c$Freq)),] # rather than modeling based on frequency in the population, model each patient type the same number of times and then bootstrap (if we increase the number of different types, may want to model the common types more than the rarer types)
# c <- distincts_with_freqs %>% filter(MOXI==1|partialmoxi==0) # an imposible combo as coded

cohort <- cohort.probs(params=params, patientvars = patientvars)

# go ahead and get the cohort freqs for sensitivity analysis also, so I can just run things once. 
setting <- "SAf"
source("bpamz_cohort.R")
SAfcohort <- cohort.probs(params=params, patientvars = patientvars)

minfreq <- 1/1e8
runcohort <- subset(cohort, cohort$Freq>minfreq|SAfcohort$Freq>minfreq)
nrow(runcohort)

## *** note that later cohort will be 2/3 pan-S HIV- new TB (1/3 sm+ and 1/3 sm-), and HIV+ will be common in SAf, so to minimize bias from a few reps, 
## I'm going to copy those 10x and split freq among them. 
## Then, 1e3 copies will give 1e4 actual copies of aything more common than ~5% of cohort.
## Thus, in a cohort desiredsize of 1e5, the largest segment will be 5e3, and 2e3 copies (4gb impact object) should be plenty. 


save(cohort, SAfcohort, minfreq, runcohort, file=paste0("cohortfreqs.",date,".Rdata"))
c <- runcohort

# and now set back to the main parametrs just to be safe
setting <- "SEA"
source("bpamz_cohort.R")

reps <- 5e3 # how many reps of each patient type to run? (will then sample Freq of each to recreate original cohort or subset thereof) -- 

impact <- list()
impact$baseline <- modelcourse(scenario = "0", c, params, reps = reps)
impact$novelrr <- modelcourse(scenario = "1a", c, params, reps = reps)
impact$novelrrx <- modelcourse(scenario = "1x", c, params, reps = reps)
impact$novelpantb <- modelcourse(scenario = "3", c, params, reps = reps)
save(impact, file = paste0("impact.",date,".Rdata"))
rm(impact)

# for supplement:
moreimpact <- list()
moreimpact$all4 <- modelcourse(scenario = "2a", c, params, reps = reps)
moreimpact$all6 <- modelcourse(scenario = "2b", c, params, reps = reps)
save(moreimpact, file = paste0("moreimpact.",date,".Rdata"))
rm(moreimpact)


# impact of xpert xdr: # will be looking at outcomes in MOXIR, Rif R and S, and results in cohort overall. All with high Xpert and different uses of fq dst. 
dst <- list()
dst$noxxdr <- modelcourse(scenario = "3x", c, params, reps = reps)
dst$stepxxdr <- modelcourse(scenario = "4x", c, params, reps = reps)
dst$fullxxdr <- modelcourse(scenario = "5x", c, params, reps = reps)
save(dst, file = paste0("dst.",date,".Rdata"))
rm(dst)




# # Sensitivity analysis: delays Don't really want to advocate for getting rid of DST, and skipping xpert would have TB detection downsides that we aren't modeling. 
# ## ...so the main question of interest may be the elimination of delays associated with needing a separate DR regimen -- 
# ## ... i.e. what additional benefit from using the same drug combination for nearly all patients. For this, would want a regimen delay,
# ## ... and want to see how each of the scenarios is affected by it, so run the whole impact (and DST?) blocks again.
# ## ... Or maybe, just the impact blocks (although novelrr wouldn't change anything))
# saveparams <- params
# allparams <- read.csv("allparams.csv", header=T, stringsAsFactors = F)
# # assuming SEA setting here
# params <- as.numeric(allparams[,2]); names(params) <- allparams[,1]; params["DSTdelay"] <- 0;
# params["Regimendelay"] <- 1 # 1 month regimendelay -- fixing but really I don't plan to use this except in sensis
# params["Tbdxtime_recurrence"] <- params["Tbdxtime_recurrenceratio"]*params["Tbdxtime"]
# delays <- list()
# delays$baseline <- modelcourse(scenario = "0", c, params, reps = reps)
# delays$novelrr <- modelcourse(scenario = "1a", c, params, reps = reps)
# delays$novelrrx <- modelcourse(scenario = "1x", c, params, reps = reps)
# delays$novelpantb <- modelcourse(scenario = "3", c, params, reps = reps)
# save(delays, file = paste0("delays.",date,".Rdata"))
# params <- saveparams 
# rm(delays)


# ADR sensitivity analysis:
params2 <- params
highparams <- as.numeric(allparams[,"UR.high"])
names(highparams) <- allparams[,1]
# use high end estimates of adr_other and adr_twodrugs
params2["adr_bpamz"] <- highparams["adr_bpamz"]
params2["adrfactor_other"] <- highparams["adrfactor_other"]
params2["adrfactor_z"] <- highparams["adrfactor_z"]
highres <- list()
highres$baseline <- modelcourse(scenario = "0", c, params2, reps = reps)
highres$novelrrx <- modelcourse(scenario = "1x", c, params2, reps = reps)
highres$noxxdr <- modelcourse(scenario = "3x", c, params2, reps = reps)
highres$fullxxdr <- modelcourse(scenario = "5x", c, params2, reps = reps)
save(highres, file = paste0("highres.",date,".Rdata"))
rm(highres)

# and then define a cohort with 2% BDQ res at baseline (noting symmetry for Pa, and worse if resistance to both esp if correlated), and assuming uncorrelated to other drugs (Worse if corerlated)
# will just model the BDQ+ cohort, then combine with original "combined" in analysis
bdqcohort <- c
bdqcohort[,"BDQ"] <- 1
bdqcohort[,"Freq"] <- c$Freq/48
combinedbdq <- list()
combinedbdq$baseline <- modelcourse(scenario = "0", bdqcohort, params, reps = reps)
combinedbdq$novelrrx <- modelcourse(scenario = "1x", bdqcohort, params, reps = reps)
combinedbdq$noxxdr <- modelcourse(scenario = "3x", bdqcohort, params, reps = reps)
combinedbdq$fullxxdr <- modelcourse(scenario = "5x", bdqcohort, params, reps = reps)
save(combinedbdq, file = paste0("combinedbdq.",date,".Rdata"))
rm(combinedbdq)




# sensitivity analysis: BPAMZ efficacy
params2 <- params
lowparams <- as.numeric(allparams[,"UR.low"])
names(lowparams) <- names(params)[1:length(lowparams)]
# actually need to also calculate a lower value for BPaMZ relative to HRZE, by using bpamhigh=T, c=equivalence at 5 of rather than 4
params2["BPaM_cxconv"] <- lowparams["BPaM_cxconv"]
params2["BPaZ_cxconv"] <- lowparams["BPaZ_cxconv"]
lowerefficacy <- list()
lowerefficacy$novelpantb <- modelcourse(scenario = "3", c, params2, reps = reps, bpamzhigh = TRUE)

# lowerefficacy$novelrrx <- modelcourse(scenario = "1x", c, params2, reps = reps, bpamzhigh = TRUE)
# lowerefficacy$noxxdr <- modelcourse(scenario = "3x", c, params2, reps = reps, bpamzhigh = TRUE)
# lowerefficacy$fullxxdr <- modelcourse(scenario = "5x", c, params2, reps = reps, bpamzhigh = TRUE)
save(lowerefficacy, file = paste0("lowerefficacy.",date,".Rdata"))
rm(lowerefficacy)
