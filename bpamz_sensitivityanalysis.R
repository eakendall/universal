# run for hi/low for each param for only 4 steps, baseline vs 3x vs 4x vs 5x, 
# and extract  cure1 for bpamz3x, baseline, and diff; and for 5x, 4x, and their diff; mean and sd for each

setwd("C:/Users/Administrator/OneDrive - Johns Hopkins University/Research/universal regimen/universal")

date <- "20190226"
setting <- "SEA"
require(plyr)
source("bpamz_cohort.R")
source("bpamz_result_utils.R")

sensitivities <- array(NA, dim=c(length(params), 32, 2))
dimnames(sensitivities) <- list("params"=names(params), 
                                "outcome"=paste0(
                                    rep(c("baselinecure1", "novelcure1", "regdiffcure1", "regratiocure1", "stepdstcure1", "fulldstcure1", "dstdiffcure1", "dstratiocure1",
                                          "baselinesuccess1", "novelsuccess1", "regdiffsuccess1", "regratiosuccess1", "stepdstsuccess1", "fulldstsuccess1", "dstdiffsuccess1", "dstratiosuccess1"), each=2),
                                    rep(c("_mean","_sd"), times=16)),
                                "value"=c("high","low"))

# for now, exclude those that are cohort characteristics (can pull from original impact, with different c)
varyparams <- !is.na(allparams$UR.low)
varyparams[1:18] <- FALSE

cohort <- cohort.probs(params=params, patientvars = patientvars)

# no need for the rare SAf rows here
minfreq <- 1/1e8
c <- subset(cohort, Freq>minfreq)

# but will need to use them later since that's the cohort used for dsts

reps <- 2e3 # how many reps of each patient type to run? (will then sample Freq of each to recreate original cohort or subset thereof) -- 


# run for hi/low for each param for only 4 steps, baseline vs 3x vs 4x vs 5x, 
for (p in (1:length(params))[varyparams]) # can use same cohort, need different impact
{
  paramshi <- paramslo <- params
  paramshi[p] <- allparams$UR.high[p]
  paramslo[p] <- allparams$UR.low[p]
  paramshi["Tbdxtime_recurrence"] <- paramshi["Tbdxtime"]*paramshi["Tbdxtime_recurrenceratio"]
  paramslo["Tbdxtime_recurrence"] <- paramslo["Tbdxtime"]*paramslo["Tbdxtime_recurrenceratio"]
  
  coursehi <- courselo <- list()
  coursehi$baseline <- shortmodelcourse(scenario = "0", c, paramshi, reps = reps)
  coursehi$noxxdr <- shortmodelcourse(scenario = "3x", c, paramshi, reps = reps)
  coursehi$stepxxdr <- shortmodelcourse(scenario = "4x", c, paramshi, reps = reps)
  coursehi$fullxxdr <- shortmodelcourse(scenario = "5x", c, paramshi, reps = reps)
  
  hicures <- loutcomeboot(individualoutcomefunction = step4cure,simoutput = coursehi, c = c)
  hisuccesses <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput = coursehi, c = c)
  
  sensitivities[p, paste0(rep(c("baselinecure1", "novelcure1", "stepdstcure1", "fulldstcure1"), each=2),rep(c("_mean","_sd"),times=4)), "high"] <- 
    unlist(lapply(hicures, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffcure1", "dstdiffcure1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply(hicures$noxxdr[,1,] - hicures$baseline[,1,], 2, mean)),
      sd(apply(hicures$noxxdr[,1,] - hicures$baseline[,1,], 2, mean)),
      mean(apply(hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,], 2, mean)),
      sd(apply(hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,], 2, mean)))

  sensitivities[p, paste0(rep(c("regratiocure1", "dstratiocure1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply((hicures$noxxdr[,1,] - hicures$baseline[,1,])/hicures$baseline[,1,], 2, mean)),
      sd(apply((hicures$noxxdr[,1,] - hicures$baseline[,1,])/hicures$baseline[,1,], 2, mean)),
      mean(apply((hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,])/hicures$stepxxdr[,1,], 2, mean)),
      sd(apply((hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,])/hicures$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("baselinesuccess1", "novelsuccess1", "stepdstsuccess1", "fulldstsuccess1"), each=2),rep(c("_mean","_sd"),times=4)), "high"] <- 
    unlist(lapply(hisuccesses, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffsuccess1", "dstdiffsuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply(hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,], 2, mean)),
      sd(apply(hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,], 2, mean)),
      mean(apply(hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply(hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,], 2, mean)))

  sensitivities[p, paste0(rep(c("regratiosuccess1", "dstratiosuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply((hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,])/hisuccesses$baseline[,1,], 2, mean)),
      sd(apply((hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,])/hisuccesses$baseline[,1,], 2, mean)),
      mean(apply((hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,])/hisuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply((hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,])/hisuccesses$stepxxdr[,1,], 2, mean)))
  
  rm(coursehi)
  rm(hicures)
  rm(hisuccesses)
  gc()
  
  courselo$baseline <- shortmodelcourse(scenario = "0", c, paramslo, reps = reps)
  courselo$noxxdr <- shortmodelcourse(scenario = "3x", c, paramslo, reps = reps)
  courselo$stepxxdr <- shortmodelcourse(scenario = "4x", c, paramslo, reps = reps)
  courselo$fullxxdr <- shortmodelcourse(scenario = "5x", c, paramslo, reps = reps)
  
  locures <- loutcomeboot(individualoutcomefunction = step4cure,simoutput = courselo, c = c)
  losuccesses <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput = courselo, c = c)
  
  sensitivities[p, paste0(rep(c("baselinecure1", "novelcure1", "stepdstcure1", "fulldstcure1"), each=2),rep(c("_mean","_sd"),times=4)), "low"] <- 
    unlist(lapply(locures, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffcure1", "dstdiffcure1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
    c(mean(apply(locures$noxxdr[,1,] - locures$baseline[,1,], 2, mean)),
      sd(apply(locures$noxxdr[,1,] - locures$baseline[,1,], 2, mean)),
      mean(apply(locures$fullxxdr[,1,] - locures$stepxxdr[,1,], 2, mean)),
      sd(apply(locures$fullxxdr[,1,] - locures$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("regratiocure1", "dstratiocure1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
    c(mean(apply((locures$noxxdr[,1,] - locures$baseline[,1,])/locures$baseline[,1,], 2, mean)),
      sd(apply((locures$noxxdr[,1,] - locures$baseline[,1,])/locures$baseline[,1,], 2, mean)),
      mean(apply((locures$fullxxdr[,1,] - locures$stepxxdr[,1,])/locures$stepxxdr[,1,], 2, mean)),
      sd(apply((locures$fullxxdr[,1,] - locures$stepxxdr[,1,])/locures$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("baselinesuccess1", "novelsuccess1", "stepdstsuccess1", "fulldstsuccess1"), each=2),rep(c("_mean","_sd"),times=4)), "low"] <- 
    unlist(lapply(losuccesses, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffsuccess1", "dstdiffsuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
    c(mean(apply(losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,], 2, mean)),
      sd(apply(losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,], 2, mean)),
      mean(apply(losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply(losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("regratiosuccess1", "dstratiosuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
    c(mean(apply((losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,])/losuccesses$baseline[,1,], 2, mean)),
      sd(apply((losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,])/losuccesses$baseline[,1,], 2, mean)),
      mean(apply((losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,])/losuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply((losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,])/losuccesses$stepxxdr[,1,], 2, mean)))
  rm(courselo)
  rm(locures)
  rm(losuccesses)
  gc()
  
  
  print(p)
  save(sensitivities, file=paste0("sensis.",date,".Rdata"))
}



load(paste0("dst.",date,".Rdata"))
smalldst <- lapply(dst, function(x) x[,,1:4,1:reps]); rm(dst)
load(paste0("impact.",date,".Rdata"))
smalldst$baseline <- impact$baseline[,,1:4,1:reps]
rm(impact)
smalldst <- smalldst[c("baseline","noxxdr","stepxxdr","fullxxdr")]
load(paste0("cohortfreqs.",date,".Rdata"))
indices <- cohort$Freq>minfreq|SAfcohort$Freq>minfreq
for (p in 1:18) # same impact but different sampling of cohort. will still need the same minfreq etc
{
  paramshi <- paramslo <- params
  paramshi[p] <- allparams$UR.high[p]
  paramslo[p] <- allparams$UR.low[p]
  paramshi["Tbdxtime_recurrence"] <- paramshi["Tbdxtime"]*paramshi["Tbdxtime_recurrenceratio"]
  paramslo["Tbdxtime_recurrence"] <- paramslo["Tbdxtime"]*paramslo["Tbdxtime_recurrenceratio"]
  
  chi <- cohort.probs(params = paramshi, patientvars = patientvars)
  chi <-chi[indices,]
  clo <- cohort.probs(params = paramslo, patientvars = patientvars)
  clo <-clo[indices,]
  
  hicures <- loutcomeboot(individualoutcomefunction = step4cure,simoutput = smalldst, c = chi)
  hisuccesses <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput = smalldst, c = chi)
  
  sensitivities[p, paste0(rep(c("baselinecure1", "novelcure1", "stepdstcure1", "fulldstcure1"), each=2),rep(c("_mean","_sd"),times=4)), "high"] <- 
    unlist(lapply(hicures, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffcure1", "dstdiffcure1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply(hicures$noxxdr[,1,] - hicures$baseline[,1,], 2, mean)),
      sd(apply(hicures$noxxdr[,1,] - hicures$baseline[,1,], 2, mean)),
      mean(apply(hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,], 2, mean)),
      sd(apply(hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("regratiocure1", "dstratiocure1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply((hicures$noxxdr[,1,] - hicures$baseline[,1,])/hicures$baseline[,1,], 2, mean)),
      sd(apply((hicures$noxxdr[,1,] - hicures$baseline[,1,])/hicures$baseline[,1,], 2, mean)),
      mean(apply((hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,])/hicures$stepxxdr[,1,], 2, mean)),
      sd(apply((hicures$fullxxdr[,1,] - hicures$stepxxdr[,1,])/hicures$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("baselinesuccess1", "novelsuccess1", "stepdstsuccess1", "fulldstsuccess1"), each=2),rep(c("_mean","_sd"),times=4)), "high"] <- 
    unlist(lapply(hisuccesses, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
  
  sensitivities[p, paste0(rep(c("regdiffsuccess1", "dstdiffsuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply(hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,], 2, mean)),
      sd(apply(hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,], 2, mean)),
      mean(apply(hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply(hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,], 2, mean)))
  
  sensitivities[p, paste0(rep(c("regratiosuccess1", "dstratiosuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "high"] <- 
    c(mean(apply((hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,])/hisuccesses$baseline[,1,], 2, mean)),
      sd(apply((hisuccesses$noxxdr[,1,] - hisuccesses$baseline[,1,])/hisuccesses$baseline[,1,], 2, mean)),
      mean(apply((hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,])/hisuccesses$stepxxdr[,1,], 2, mean)),
      sd(apply((hisuccesses$fullxxdr[,1,] - hisuccesses$stepxxdr[,1,])/hisuccesses$stepxxdr[,1,], 2, mean)))
  
  rm(hicures)
  rm(hisuccesses)
  gc()
  
  
    locures <- loutcomeboot(individualoutcomefunction = step4cure,simoutput = smalldst, c = clo)
    losuccesses <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput = smalldst, c = c)
    
    sensitivities[p, paste0(rep(c("baselinecure1", "novelcure1", "stepdstcure1", "fulldstcure1"), each=2),rep(c("_mean","_sd"),times=4)), "low"] <- 
      unlist(lapply(locures, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
    
    sensitivities[p, paste0(rep(c("regdiffcure1", "dstdiffcure1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
      c(mean(apply(locures$noxxdr[,1,] - locures$baseline[,1,], 2, mean)),
        sd(apply(locures$noxxdr[,1,] - locures$baseline[,1,], 2, mean)),
        mean(apply(locures$fullxxdr[,1,] - locures$stepxxdr[,1,], 2, mean)),
        sd(apply(locures$fullxxdr[,1,] - locures$stepxxdr[,1,], 2, mean)))
    
    sensitivities[p, paste0(rep(c("regratiocure1", "dstratiocure1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
      c(mean(apply((locures$noxxdr[,1,] - locures$baseline[,1,])/locures$baseline[,1,], 2, mean)),
        sd(apply((locures$noxxdr[,1,] - locures$baseline[,1,])/locures$baseline[,1,], 2, mean)),
        mean(apply((locures$fullxxdr[,1,] - locures$stepxxdr[,1,])/locures$stepxxdr[,1,], 2, mean)),
        sd(apply((locures$fullxxdr[,1,] - locures$stepxxdr[,1,])/locures$stepxxdr[,1,], 2, mean)))
    
    sensitivities[p, paste0(rep(c("baselinesuccess1", "novelsuccess1", "stepdstsuccess1", "fulldstsuccess1"), each=2),rep(c("_mean","_sd"),times=4)), "low"] <- 
      unlist(lapply(losuccesses, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))
    
    sensitivities[p, paste0(rep(c("regdiffsuccess1", "dstdiffsuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
      c(mean(apply(losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,], 2, mean)),
        sd(apply(losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,], 2, mean)),
        mean(apply(losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,], 2, mean)),
        sd(apply(losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,], 2, mean)))
    
    sensitivities[p, paste0(rep(c("regratiosuccess1", "dstratiosuccess1"), each=2),rep(c("_mean","_sd"),times=2)), "low"] <- 
      c(mean(apply((losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,])/losuccesses$baseline[,1,], 2, mean)),
        sd(apply((losuccesses$noxxdr[,1,] - losuccesses$baseline[,1,])/losuccesses$baseline[,1,], 2, mean)),
        mean(apply((losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,])/losuccesses$stepxxdr[,1,], 2, mean)),
        sd(apply((losuccesses$fullxxdr[,1,] - losuccesses$stepxxdr[,1,])/losuccesses$stepxxdr[,1,], 2, mean)))
    
    rm(locures)
    rm(losuccesses)
    gc()
    
    
  print(p)
  save(sensitivities, file=paste0("sensis.",date,".Rdata"))
}

write.table(sensitivities, file = paste0("sensitivities.",date,".csv"), sep=",")

# get baseline
cures <- loutcomeboot(individualoutcomefunction = step4cure,simoutput = smalldst, c = cohort[indices,])
successes <- loutcomeboot(individualoutcomefunction = treatmentsuccess,simoutput = smalldst, c = cohort[indices,])

baseline <- numeric(32)
names(baseline) <- dimnames(sensitivities)[[2]]

baseline[paste0(rep(c("baselinecure1", "novelcure1", "stepdstcure1", "fulldstcure1"), each=2),rep(c("_mean","_sd"),times=4))] <- 
  unlist(lapply(cures, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))

baseline[paste0(rep(c("regdiffcure1", "dstdiffcure1"), each=2),rep(c("_mean","_sd"),times=2))] <- 
  c(mean(apply(cures$noxxdr[,1,] - cures$baseline[,1,], 2, mean)),
    sd(apply(cures$noxxdr[,1,] - cures$baseline[,1,], 2, mean)),
    mean(apply(cures$fullxxdr[,1,] - cures$stepxxdr[,1,], 2, mean)),
    sd(apply(cures$fullxxdr[,1,] - cures$stepxxdr[,1,], 2, mean)))

baseline[paste0(rep(c("regratiocure1", "dstratiocure1"), each=2),rep(c("_mean","_sd"),times=2))] <- 
  c(mean(apply((cures$noxxdr[,1,] - cures$baseline[,1,])/cures$baseline[,1,], 2, mean)),
    sd(apply((cures$noxxdr[,1,] - cures$baseline[,1,])/cures$baseline[,1,], 2, mean)),
    mean(apply((cures$fullxxdr[,1,] - cures$stepxxdr[,1,])/cures$stepxxdr[,1,], 2, mean)),
    sd(apply((cures$fullxxdr[,1,] - cures$stepxxdr[,1,])/cures$stepxxdr[,1,], 2, mean)))

baseline[paste0(rep(c("baselinesuccess1", "novelsuccess1", "stepdstsuccess1", "fulldstsuccess1"), each=2),rep(c("_mean","_sd"),times=4))] <- 
  unlist(lapply(successes, function(x) c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean)))))

baseline[paste0(rep(c("regdiffsuccess1", "dstdiffsuccess1"), each=2),rep(c("_mean","_sd"),times=2))] <- 
  c(mean(apply(successes$noxxdr[,1,] - successes$baseline[,1,], 2, mean)),
    sd(apply(successes$noxxdr[,1,] - successes$baseline[,1,], 2, mean)),
    mean(apply(successes$fullxxdr[,1,] - successes$stepxxdr[,1,], 2, mean)),
    sd(apply(successes$fullxxdr[,1,] - successes$stepxxdr[,1,], 2, mean)))

baseline[paste0(rep(c("regratiosuccess1", "dstratiosuccess1"), each=2),rep(c("_mean","_sd"),times=2))] <- 
  c(mean(apply((successes$noxxdr[,1,] - successes$baseline[,1,])/successes$baseline[,1,], 2, mean)),
    sd(apply((successes$noxxdr[,1,] - successes$baseline[,1,])/successes$baseline[,1,], 2, mean)),
    mean(apply((successes$fullxxdr[,1,] - successes$stepxxdr[,1,])/successes$stepxxdr[,1,], 2, mean)),
    sd(apply((successes$fullxxdr[,1,] - successes$stepxxdr[,1,])/successes$stepxxdr[,1,], 2, mean)))


baseline <- t(rbind(baseline, baseline))

s  <- abind(sensitivities, baseline, along=1); 
dimnames(s) <- list(c(dimnames(sensitivities)[[1]],"baseline"), dimnames(sensitivities)[[2]], dimnames(sensitivities)[[3]])

save(s, file = paste0("sensisbaseline.",date,".Rdata"))
write.table(s, file = paste0("sensitivitieswithbaseline.",date,".csv"), sep=",")



# l4 <- loutcomeboot(individualoutcomefunction = step4cure, simoutput = smalldst, c = c)
# cure1 <- unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,],2,mean)),
#                                         sd(apply(y[,1,], 2, mean))), 3)))
# cure1  



######### sensitivity analysis, all params one-way
load("sensisbaseline.20190226.Rdata",verbose = T)

cbind(s[,"regdiffcure1_mean",],s["baseline","regdiffcure1_mean",1])
cbind(s[,"regdiffcure1_mean",],s["baseline","regdiffcure1_mean",1])

varied <- !is.na(s[,"regdiffcure1_mean",1]); varied[length(varied)] <- FALSE

par(mar=c(3,10,3,3))
barplot(s[varied,"baselinecure1_mean",1] - s["baseline","baselinecure1_mean",1] , 
        horiz = TRUE, col="white", xlim=c(-0.05,0.05), las=1, xaxt='n')
barplot(s[varied,"baselinecure1_mean",2] - s["baseline","baselinecure1_mean",1] , 
        horiz = TRUE, col="black", add=T, las=1, xaxt='n', yaxt='n')

par(mar=c(3,10,3,3))
barplot(s[varied,"novelcure1_mean",1] - s["baseline","novelcure1_mean",1] , 
        horiz = TRUE, col="white", xlim=c(-0.05,0.05), las=1, xaxt='n')
barplot(s[varied,"novelcure1_mean",2] - s["baseline","novelcure1_mean",1] , 
        horiz = TRUE, col="black", add=T, las=1, xaxt='n', yaxt='n')

barplot(s[varied,"novelcure1_mean",1] - s[varied,"baselinecure1_mean",1] - (s["baseline","novelcure1_mean",1] - s["baseline","baselinecure1_mean",1] ), 
        horiz = TRUE, col="white", xlim=c(-0.05,0.05), las=1, xaxt='n')
barplot(s[varied,"novelcure1_mean",2] - s[varied,"baselinecure1_mean",2] - (s["baseline","novelcure1_mean",1] - s["baseline","baselinecure1_mean",1] ), 
        horiz = TRUE, col="black", add=T, las=1, xaxt='n', yaxt='n')


par(mar=c(3,10,3,3))
barplot(s[varied,"regdiffcure1_mean",1] - s["baseline","regdiffcure1_mean",1] , 
        horiz = TRUE, col="white", xlim=c(-0.05,0.05), las=1, xaxt='n')
barplot(s[varied,"regdiffcure1_mean",2] - s["baseline","regdiffcure1_mean",2] , 
        horiz = TRUE, col="black", add=T, las=1, xaxt='n', yaxt='n')

par(mar=c(3,10,3,3))
barplot(s[varied,"dstdiffcure1_mean",1] - s["baseline","dstdiffcure1_mean",1] , 
        horiz = TRUE, col="white", xlim=c(-0.05,0.05), las=1, xaxt='n')
barplot(s[varied,"dstdiffcure1_mean",2] - s["baseline","dstdiffcure1_mean",2] , 
        horiz = TRUE, col="black", add=T, las=1, xaxt='n', yaxt='n')



