# setwd("C:/Users/ekendal2/OneDrive - Johns Hopkins University/Research/universal regimen/bpamz")

# 0708: changing defaults for fq dst analysis (second manuscript) to 4 months for all who get bpamz (even if rr detected), with 6mo for rr as a sensitivity analysis

date <- "20190711"
set.seed(3423235)
setting <- "SAf"
source("bpamz_cohort.R")

setting <- "SAf"
SAfcohort <- cohort.probs(params=params, patientvars = patientvars)

# running analysis for first paper with South Africa's current Xpert and RIF-R,
# and will compare to lower xpert (1% for new, and half current for rerx), and/or higher RIF-R
# higher RIF-R will require resampling of cohort, say 3x higher ODDS in both new and retreatment (->9% and 21%), with same prev moxi and pza in rr and in rs as now, same correlation bewteen them, same inh in rr as now, and same prev inhmono.
increaseodds <- function(p, oddsfactor)
{
  odds <- p/(1-p)
  newodds <- odds*oddsfactor
  newp <- newodds/(1+newodds)
  return(newp)
}

# go ahead and get the cohort freqs for sensitivity analysis also, so I can just run things once. 
# note that I'll assume the same max Xpert coverage for both settings in the DST paper, so can use the same dst runs and weight by different cohorts
setting <- "SEA"
source("bpamz_cohort.R")
SEAcohort <- cohort.probs(params=params, patientvars = patientvars)

minfreq <- 1/1e8
runcohort <- SAfcohort[SAfcohort$Freq>minfreq|SEAcohort$Freq>minfreq , ]
nrow(runcohort)

## *** note that later cohort will be 2/3 pan-S HIV- new TB (1/3 sm+ and 1/3 sm-), and HIV+ will be common in SAf, so to minimize bias from a few reps, 
## I'm going to copy those 10x and split freq among them. 
## Then, 1e3 copies will give 1e4 actual copies of aything more common than ~5% of cohort.
## Thus, in a cohort desiredsize of 1e5, the largest segment will be 5e3, and 2e3 copies (4gb impact object) should be plenty. I'll double that to be sure.

save(SEAcohort, SAfcohort, minfreq, runcohort, file=paste0("cohortfreqs.",date,".Rdata"))

# reset setting just to be safe
setting <- "SAf"
source("bpamz_cohort.R")

reps <- 5e3 # how many reps of each patient type to run? (will then sample Freq of each to recreate original cohort or subset thereof) -- 

impact <- list()
impact$baseline <- modelcourse(scenario = "0", runcohort, params, reps = reps)
impact$baseline_lowx <- modelcourse(scenario = "0o", runcohort, params, reps = reps)
impact$novelrr <- modelcourse(scenario = "1", runcohort, params, reps = reps)
impact$novelrr_lowx <- modelcourse(scenario = "1o", runcohort, params, reps = reps)
impact$novelpantb <- modelcourse(scenario = "3", runcohort, params, reps = reps)
impact$novelpantb_lowx <- modelcourse(scenario = "3o", runcohort, params, reps = reps)
save(impact, file = paste0("impact.",date,".Rdata"))
rm(impact)

# # these next two are for second paper:
# # fixed durations, no dst
# moreimpact <- list()
# moreimpact$all4 <- modelcourse(scenario = "2a", runcohort, params, reps = reps)
# moreimpact$all6 <- modelcourse(scenario = "2b", runcohort, params, reps = reps)
# moreimpact$baseline_highx <- modelcourse(scenario = "0x", runcohort, params, reps = reps)
# save(moreimpact, file = paste0("moreimpact.",date,".Rdata"))
# rm(moreimpact)


# # Sensitivity analysis for first paper:

# ADR sensitivity analysis:
params2 <- params
highparams <- as.numeric(allparams[,"seahigh"])
if(setting=="SAf") highparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"safhigh"][!is.na(allparams[,"saf"])])
names(highparams) <- allparams[,1]
# use high end estimates of adr_other and adr_twodrugs
params2["adr_bpamz"] <- highparams["adr_bpamz"]
params2["adrfactor_other"] <- highparams["adrfactor_other"]
params2["adrfactor_z"] <- highparams["adrfactor_z"]
highres <- list()
highres$baseline <- modelcourse(scenario = "0", runcohort, params2, reps = reps)
highres$novelrr <- modelcourse(scenario = "1", runcohort, params2, reps = reps)
highres$novelpantb <- modelcourse(scenario = "3", runcohort, params2, reps = reps)
save(highres, file = paste0("highres.",date,".Rdata"))
rm(highres)

# and then define a cohort with 2% BDQ res at baseline (noting symmetry for Pa, and worse if resistance to both esp if correlated), and assuming uncorrelated to other drugs (Worse if corerlated)
# will just model the BDQ+ cohort, then combine with original "combined" in analysis
bdqcohort <- runcohort
bdqcohort[,"BDQ"] <- 1
bdqcohort[,"Freq"] <- runcohort$Freq/48 # should have beeen 49, but i'll adjust for this later -- and actually, it doesn't matter because we just use Freq for probs and set the desiredsize separately.
combinedbdq <- list()
combinedbdq$baseline <- modelcourse(scenario = "0", bdqcohort, params, reps = reps)
combinedbdq$novelrr <- modelcourse(scenario = "1", bdqcohort, params, reps = reps)
combinedbdq$novelpantb <- modelcourse(scenario = "3", bdqcohort, params, reps = reps)
save(combinedbdq, file = paste0("combinedbdq.",date,".Rdata"))
rm(combinedbdq)



# sensitivity analysis: BPAMZ efficacy
params2 <- params
lowparams <- as.numeric(allparams[,"sealow"])
if(setting=="SAf") lowparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"saflow"][!is.na(allparams[,"saf"])])
names(lowparams) <- names(params)[1:length(lowparams)]
# actually need to also calculate a lower value for BPaMZ relative to HRZE, by using bpamhigh=T, c=equivalence at 5 of rather than 4
params2["BPaM_cxconv"] <- lowparams["BPaM_cxconv"]
params2["BPaZ_cxconv"] <- lowparams["BPaZ_cxconv"]
lowerefficacy <- list()
lowerefficacy$novelrr <- modelcourse(scenario = "1", runcohort, params2, reps = reps, bpamzhigh = TRUE)
lowerefficacy$novelrr_lowx <- modelcourse(scenario = "1o", runcohort, params2, reps = reps, bpamzhigh = TRUE)
lowerefficacy$novelpantb <- modelcourse(scenario = "3", runcohort, params2, reps = reps, bpamzhigh = TRUE)
lowerefficacy$novelpantb_lowx <- modelcourse(scenario = "3o", runcohort, params2, reps = reps, bpamzhigh = TRUE)
save(lowerefficacy, file = paste0("lowerefficacy.",date,".Rdata"))
rm(lowerefficacy)


# sensitivity analysis: DR SOC
# improve to the point that DR Rx of FQ-S, if completed, has same cure prob as HRZE (for HR-S) if completed. But still lasts 18 months rather than 6. 
# and FQ-R carries same higher RR of fail/relapse
params2 <- params
params2["MDR_failrelapse_FQ-S"] <- (1+params2["Failures_per_recurrence"])*make.recurrence.matrix(params = params)["HR(ZE)","6"]
params2["MDR_failrelapse_FQ-R"] <- params2["MDR_failrelapse_FQ-S"]*(params["MDR_failrelapse_FQ-R"]/params["MDR_failrelapse_FQ-S"])
regimendurations <- c(6, 12, 4, 6, 6)
betterSOC <- list()
betterSOC$baseline <- modelcourse(scenario = "0", runcohort, params2, reps = reps)
betterSOC$baseline_lowx <- modelcourse(scenario = "0o", runcohort, params2, reps = reps)
betterSOC$novelrr <- modelcourse(scenario = "1", runcohort, params2, reps = reps)
betterSOC$novelrr_lowx <- modelcourse(scenario = "1o", runcohort, params2, reps = reps)
betterSOC$novelpantb <- modelcourse(scenario = "3", runcohort, params2, reps = reps)
betterSOC$novelpantb_lowx <- modelcourse(scenario = "3o", runcohort, params2, reps = reps)
save(betterSOC, file = paste0("betterSOC.",date,".Rdata"))
rm(betterSOC)



  # 


#### for dst paper #############
# addition of xpert xdr: # will be looking at outcomes in MOXIR, Rif R and S, and results in cohort overall. All with high Xpert and different uses of fq dst. 
dst <- list()
dst$noxxdr <- modelcourse(scenario = "3xa", runcohort, params, reps = reps)
dst$stepxxdr <- modelcourse(scenario = "4xa", runcohort, params, reps = reps)
dst$fullxxdr <- modelcourse(scenario = "5xa", runcohort, params, reps = reps)
dst$baseline <- modelcourse(scenario = "0x", runcohort, params, reps = reps)
save(dst, file = paste0("dst.",date,".Rdata"))
rm(dst)
dst6m <- list()
dst6m$noxxdr <- modelcourse(scenario = "3x", runcohort, params, reps = reps)
dst6m$stepxxdr <- modelcourse(scenario = "4x", runcohort, params, reps = reps)
dst6m$fullxxdr <- modelcourse(scenario = "5x", runcohort, params, reps = reps)
save(dst6m, file = paste0("dst6m.",date,".Rdata"))
rm(dst6m)


# Robustness and ADR sensivity analysis for DST paper:
params2 <- params
highparams <- as.numeric(allparams[,"seahigh"])
lowparams <- as.numeric(allparams[,"sealow"])
if(setting=="SAf") 
  { highparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"safhigh"][!is.na(allparams[,"saf"])])
    lowparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"saflow"][!is.na(allparams[,"saf"])]) }
names(highparams) <- allparams[,1]
names(lowparams) <- names(params)[1:length(lowparams)]
# use high end estimates of adr_other and adr_twodrugs (but not for baseline adr_bpamz anymore), and for efficacy BPaZ_cxconv and INH_multiplier
# params2["adr_bpamz"] <- highparams["adr_bpamz"]
params2["adrfactor_other"] <- highparams["adrfactor_other"]
params2["adrfactor_twodrugs"] <- highparams["adrfactor_twodrugs"]
# params2["adrfactor_z"] <- highparams["adrfactor_z"]
params2["BPaZ_cxconv"] <- lowparams["BPaZ_cxconv"]
params2["INH_multiplier"] <- highparams["INH_multiplier"]
robustdst <- list()
robustdst$noxxdr <- modelcourse(scenario = "3xa", runcohort, params2, reps = reps)
robustdst$stepxxdr <- modelcourse(scenario = "4xa", runcohort, params2, reps = reps)
robustdst$fullxxdr <- modelcourse(scenario = "5xa", runcohort, params2, reps = reps)
save(robustdst, file = paste0("robustdst.",date,".Rdata"))


params2 <- params
highparams <- as.numeric(allparams[,"seahigh"])
if(setting=="SAf") highparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"safhigh"][!is.na(allparams[,"saf"])])
names(highparams) <- allparams[,1]
# use high end estimates of adr_other and adr_twodrugs
params2["adr_bpamz"] <- highparams["adr_bpamz"]
params2["adrfactor_other"] <- highparams["adrfactor_other"]
params2["adrfactor_z"] <- highparams["adrfactor_z"]
highresdst <- list()
highresdst$noxxdr <- modelcourse(scenario = "3xa", runcohort, params2, reps = reps)
highresdst$stepxxdr <- modelcourse(scenario = "4xa", runcohort, params2, reps = reps)
highresdst$fullxxdr <- modelcourse(scenario = "5xa", runcohort, params2, reps = reps)
save(highresdst, file = paste0("highresdst.",date,".Rdata"))
rm(highresdst)


# hypermutator sensitivity analysis: Will combined the initially-resistant individuals from this population with the initially-pan-s from the main simulation.
params2 <- params
highparams <- as.numeric(allparams[,"seahigh"])
lowparams <- as.numeric(allparams[,"sealow"])
if(setting=="SAf") 
{ highparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"safhigh"][!is.na(allparams[,"saf"])])
lowparams[!is.na(allparams[,"saf"])] <- as.numeric(allparams[,"saflow"][!is.na(allparams[,"saf"])]) }
names(highparams) <- allparams[,1]
names(lowparams) <- names(params)[1:length(lowparams)]
params2["adrfactor_other"] <- highparams["adrfactor_other"]
params2["adrfactor_twodrugs"] <- highparams["adrfactor_twodrugs"]
params2["adrfactor_z"] <- highparams["adrfactor_z"]
params2["adrfactor_partialMFX"] <- highparams["adrfactor_partialMFX"]
params2["BPaZ_cxconv"] <- lowparams["BPaZ_cxconv"]
params2["INH_multiplier"] <- highparams["INH_multiplier"]
hypermutator <- list()
hypermutator$noxxdr <- modelcourse(scenario = "3xa", runcohort, params2, reps = reps)
hypermutator$stepxxdr <- modelcourse(scenario = "4xa", runcohort, params2, reps = reps)
hypermutator$fullxxdr <- modelcourse(scenario = "5xa", runcohort, params2, reps = reps)
save(hypermutator, file = paste0("hypermutator.",date,".Rdata"))


