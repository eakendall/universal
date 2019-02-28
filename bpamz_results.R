setwd("C:/Users/Administrator/OneDrive - Johns Hopkins University/Research/universal regimen/universal")

date <- "20190226"
setting <- "SEA"
require(plyr)
source("bpamz_cohort.R")
source("bpamz_result_utils.R")

load(file = paste0("impact.",date,".Rdata"))

setting <- "SEA"
load(paste0("cohortfreqs.",date,".Rdata"))
c <- runcohort

# make bpamz outcome table
tableSreg <- function(impact, c, save)
{
# #  proportion cured after 1st and 2nd round:
  l4 <- loutcomeboot(individualoutcomefunction = step4812cure, simoutput = impact, c = c)
  if(save) save(l4, file=paste0("l4.",date,".Rdata"))
  # load(file=paste0("l4.",date,".Rdata"))
  # cured first round: overall, rif-s, rif-r, and fq-r:
  cure1 <- rbind(
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,],2,mean)),
                                 sd(apply(y[,1,], 2, mean))), 3))),
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum)),
                                 sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum))), 3))),
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum)),
                                          sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum))), 3)) ), 
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum)),
                                 sd(apply(y[,1,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum))), 3))) )
  rownames(cure1) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure1
  
  # cured second round: overall, rif-s, rif-r, and fq-r:
  cure2 <- rbind(
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,2,],2,mean)),
                                          sd(apply(y[,2,], 2, mean))), 3))),
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum)),
                                          sd(apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum))), 3))),
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum)),
                                          sd(apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum))), 3)) ), 
    unlist(lapply(l4, function(y) round(c(mean(apply(y[,2,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum)),
                                          sd(apply(y[,2,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum))), 3))) )
  rownames(cure2) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure2

  #cured by N=36 months
  # lc30 <- loutcomeboot(individualoutcomefunction = still.in.state, states=7, t=30, simoutput=impact, c=c)
  lc36 <- loutcomeboot(individualoutcomefunction = still.in.state, states=7, t=36, simoutput=impact, c=c)
  if(save) save(lc36, file=paste0("lc36.",date,".Rdata"))
  # load(paste0("lc36.",date,".Rdata"))
  
  cure36m <- rbind(
    unlist(lapply(lc36, function(y) round(c(mean(apply(y[,1,],2,mean)),
                                            sd(apply(y[,1,],2,mean))), 3))),
    unlist(lapply(lc36, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                             sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum))), 3))),
    unlist(lapply(lc36, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                            sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum))), 3))),
    unlist(lapply(lc36, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                         sd(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum))), 3))) )
  rownames(cure36m) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure36m
  
  # cure first round if treated
  lsuccess1 <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=impact, c=c)
  if(save) save(lsuccess1, file=paste0("lsuccess1.",date,".Rdata"))

  success <- rbind(
    unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,],2,sum)/apply(y[,1,],2,sum)),
                                            sd(apply(y[,2,],2,sum)/apply(y[,1,],2,sum))), 3))),
    unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,]*(y[,3,]==0),2,sum)/apply(y[,1,]*(y[,3,]==0),2,sum)),
                                            sd(apply(y[,2,]*(y[,3,]==0),2,sum)/apply(y[,1,]*(y[,3,]==0),2,sum))), 3))),
    unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,]*(y[,3,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1),2,sum)),
                                                 sd(apply(y[,2,]*(y[,3,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1),2,sum))), 3))),
    unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                 sd(apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum))), 3))) )
  rownames(success) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  success
  
  
  

  # total TB/infectious time:
  li <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                     simoutput=impact, c=c)
  if(save) save(li, file=paste0("li.",date,".Rdata"))
  # load(file=paste0("li.",date,".Rdata"))
  tbtime <- rbind(
    unlist(lapply(li, function(x) round(c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean))),3))), 
    unlist(lapply(li, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                    sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum))),3))),
    unlist(lapply(li, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                    sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum))),3))),
    unlist(lapply(li, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                         sd(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum))),3)))
  )
  rownames(tbtime) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  tbtime


  # time on treatment (of some kind) 
  followyears <- 20
  ltall <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(statetypes$treating, statetypes$treating_adr), cutofftime=12*followyears, carryforward=T, 
                        simoutput=impact, c=c)
  if(save) save(ltall, file=paste0("ltall.",date,".Rdata"))
  # load(file=paste0("ltall.",date,".Rdata"))
  
  rxtime <- rbind(
    unlist(lapply(ltall, function(x) round(c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean))),3))), 
    unlist(lapply(ltall, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                          sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum))),3))),
    unlist(lapply(ltall, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                          sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum))),3))),
    unlist(lapply(ltall, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                          sd(apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum))),3)))
  )
  rownames(rxtime) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  rxtime


  return(rbind(cure1, cure2, cure36m, success, tbtime, rxtime))
}


# make table showing outcomes, summary of full generation and recurrent matrix and adr matrix which will be included in appendix. 
rec <- make.recurrence.matrix()[c("HR(ZE)","R(ZE)","BPaMZ","BPamZ","BPaM","BPaZ", "BPam","BPa","MDR, FQ-S", "MDR, FQ-R"),c("4","6","18")]
# (rec <- round(rbind((rec[1,]*sum(c$Freq*(1-c$RIF)*(1-c$INH))+rec[2,]*sum(c$Freq*(1-c$RIF)*c$INH))/sum(c$Freq*(1-c$RIF)),
#       rec[3:nrow(rec),]),3)
# )
# to check: all displayed entries should be less than this, else need a min of the two:
round(rowSums(make.adr.matrix()[c("HR(ZE)","R(ZE)","BPaMZ","BPamZ","BPaM","BPaZ", "BPam","BPa","MDR, FQ-S", "MDR, FQ-R"),], na.rm=T)*100, 1)
write.table(round(rec,3), "recurrencetable.csv", sep=",", row.names = TRUE, col.names = TRUE)

# details for supplement:
write.table(rbind(unlist(set.scenario("0")), unlist(set.scenario("1a")), unlist(set.scenario("1x")),unlist(set.scenario("3")),
                  unlist(set.scenario("3x")),unlist(set.scenario("4x")), unlist(set.scenario("5x"))), file="scenarios.csv", sep=",")
drugs
make.active.regimen.matrix()[,1:(length(drugs)+length(regimens))] # describe, and refer to published code
a <- make.recurrence.matrix(); a[c(1:4,8:nrow(a)),c("9","12","18")] <- NA; round(a,3)
write.table(a, file="recurrencematrix.csv", sep=",")
a <- make.adr.matrix(); a[a==0]<-NA; round(a,3)
write.table(a, file="adrmatrix.csv", sep=",")




tSEA <- tableSreg(impact=impact, c=runcohort, save=T)
save(tSEA, file=paste0("tSEA.",date,".Rdata"))
write.table(tSEA, file=paste0("tSEA.",date,".csv"), sep = ",")
tSEA[1:16,] <- round(tSEA[1:16,]*100,1)
tSEA[17:24,] <- round(tSEA[17:24,],2)

tSEA2 <- cbind(paste0(tSEA[,1], " \u00b1 ",tSEA[,2]), paste0(tSEA[,3], " \u00b1 ",tSEA[,4]), paste0(tSEA[,5], " \u00b1 ",tSEA[,6]), paste0(tSEA[,7], " \u00b1 ",tSEA[,8]))
Encoding(tSEA2)<-"UTF-8"
tSEA2[1:16,] <- paste0(tSEA2[1:16,],"%")
tSEA2
write.table(tSEA2, file=paste0("tSEA2_success.",date,".csv"), sep = ",")
save(tSEA2, file=paste0("tSEA2.",date,".Rdata"))

# additional results for text
# % with initial regimen as hrze for rif-r in expanded xpert:
dxrr <- function(patiententry) {c(patiententry["Currentregimen",3]==4,patiententry["Currentregimen",3]>0,patiententry["RIF",1])} # treated rifr, treated at all, has rifr

lreg <- loutcomeboot(individualoutcomefunction = dxrr, simoutput = list("novelrr"=impact$novelrr, "novelrrx"=impact$novelrrx), c = c)
save(lreg, file=paste0("lreg.",date,".Rdata"))
# load("lreg.20190226.Rdata")
lapply(lreg, function(y) mean(apply(y[,1,]*y[,2,]*y[,3,],2,sum)/apply(y[,2,]*y[,3,],2,sum)))

# proportion of overall cohort cured in 3 rounds:
load(paste0("l4.",date,".Rdata"))
lapply(l4, function(y) mean(apply(y[,3,],2,mean)))
lapply(l4, function(y) sd(apply(y[,3,],2,mean)))
rm(l4)

# % reduction in time with TB
load("li.20190226.Rdata")
lapply(li, function(x) apply(x[,1,], 2, mean))
lapply(li, function(x) mean(apply(x[,1,]*(x[,which(patientvars=="RIF")+2,]==1), 2, sum)/apply((x[,which(patientvars=="RIF")+2,]==1), 2, sum)))
lapply(li, function(x) sd(apply(x[,1,]*(x[,which(patientvars=="RIF")+2,]==1), 2, sum)/apply((x[,which(patientvars=="RIF")+2,]==1), 2, sum)))

reduction <- lapply(li, function(x) apply((li$baseline[,1,]-x[,1,]), 2, sum)/apply(li$baseline[,1,],2,sum))
lapply(reduction, mean); lapply(reduction, sd)
rrreduction <- lapply(li, function(x) (apply(li$baseline[,1,]*(li$baseline[,which(patientvars=="RIF")+2,]==1),2,sum)-
                                              apply(x[,1,]*(x[,which(patientvars=="RIF")+2,]==1), 2, sum))/
                                                apply(li$baseline[,1,]*(li$baseline[,which(patientvars=="RIF")+2,]==1),2,sum))
univrrreduction <- lapply(li, function(x) (apply(li$novelrr[,1,]*(li$novelrr[,which(patientvars=="RIF")+2,]==1),2,sum)-
                                         apply(x[,1,]*(x[,which(patientvars=="RIF")+2,]==1), 2, sum))/
                        apply(li$novelrr[,1,]*(li$novelrr[,which(patientvars=="RIF")+2,]==1),2,sum))
rrreduction
lapply(rrreduction, mean); lapply(rrreduction, sd)
lapply(univrrreduction, mean); lapply(univrrreduction, sd)
rm(li)

# for supplement, South Africa
minfreq <- 1/1e8
SAfc <- SAfcohort[cohort$Freq>minfreq|SAfcohort$Freq>minfreq,]
nrow(SAfc)
## and now run the above to amke new table (without saving)
tSAf <- tableSreg(impact=impact, c=SAfc, save=F)
# save(tSAf, file=paste0("tSAf.",date,".Rdata"))
write.table(tSAf, file=paste0("tSAf.",date,".csv"), sep = ",")
save(tSAf, file=paste0("tSAf.",date,".Rdata"))
tSAf[1:16,] <- round(tSAf[1:16,]*100,1)
tSAf[17:24,] <- round(tSAf[17:24,],2)
tSAf2 <- cbind(paste0(tSAf[,1], " \u00b1 ",tSAf[,2]), paste0(tSAf[,3], " \u00b1 ",tSAf[,4]), paste0(tSAf[,5], " \u00b1 ",tSAf[,6]), paste0(tSAf[,7], " \u00b1 ",tSAf[,8]))
tSAf2[1:16,] <- paste0(tSAf2[1:16,],"%")
Encoding(tSAf2)<-"UTF-8"
write.table(tSAf2, file=paste0("tSAf2.",date,".csv"), sep = ",")


# for supplement, same duration for all
load(file=paste0("moreimpact.",date,".Rdata")) # 4 or 6 months for all
tmore <- tableSreg(impact=moreimpact, c=runcohort, save=F)
save(tmore, file=paste0("tmore.",date,".Rdata"))
load(paste0("tSEA.",date,".Rdata"))
tmore <- cbind(tSEA[,7:8], tmore)
write.table(tmore, file=paste0("tmore",date,".csv"), sep = ",")
tmore[1:16,] <- round(tmore[1:16,]*100,1)
tmore[17:24,] <- round(tmore[17:24,],2)
tmore2 <- cbind(paste0(tmore[,1], " \u00b1 ",tmore[,2]), paste0(tmore[,3], " \u00b1 ",tmore[,4]), paste0(tmore[,5], " \u00b1 ",tmore[,6]))
tmore2[1:16,] <- paste0(tmore2[1:16,],"%")
Encoding(tmore2)<-"UTF-8"
write.table(tmore2, file=paste0("tmore2.",date,".csv"), sep = ",")
save(tmore2, file=paste0("tmore2.",date,".Rdata"))

rm(moreimpact)


# for sensitivity analysis, change prevalence of moxi/PZA resistance: 
# Higher FQ R setting -- will be putting in supplement for part 1, but including in Fig for part 2 (and in text for part 3?)
params2 <- params
highparams <- as.numeric(allparams[,"UR.high"])
names(highparams) <- allparams[,1]
# use high end estimates of FQ-R for SEA
params2["MOXI-R-any_in_RIF-S"] <- highparams["MOXI-R-any_in_RIF-S"]
params2["MOXI-R-highlevel_in_RIF-S"] <- highparams["MOXI-R-highlevel_in_RIF-S"]
highFQScohort <- cohort.probs(params = params2, patientvars = patientvars)
# will use the same "impact" and "dst", but this different highFQcohort, restricted to same patient types as base model since those are the ones we ran.
highFQScohort <- highFQScohort[cohort$Freq>minfreq|SAfcohort$Freq>minfreq,]
# only care about some of the regimen impact outcomes:
#(and also, don't want to emphasize the hgher cure in rif-r due to more prior treatment and faster rediagnosis)
# l4 (proportion cured in 1 and 2 rounds), and li
# overall cohort, pantb,
# for higher RS MOXI-R, higher RS and RR MOXI-R, and higher MOXI-R and PZA+R in MOXI-R 
# so need three new cohorts:
# highFQScohort as defined above, and
params2["MOXI-R-any_in_RIF-R"] <- highparams["MOXI-R-any_in_RIF-R"]
params2["MOXI-R-highlevel_in_RIF-R"] <- highparams["MOXI-R-highlevel_in_RIF-R"]
highFQRcohort <- cohort.probs(params = params2, patientvars = patientvars)
highFQRcohort <- highFQRcohort[cohort$Freq>minfreq|SAfcohort$Freq>minfreq,]
params2["PZA-R_in_RIF-R"] <- highparams["PZA-R_in_RIF-R"]
highFQZRcohort <- cohort.probs(params = params2, patientvars = patientvars)
highFQZRcohort <- highFQZRcohort[cohort$Freq>minfreq|SAfcohort$Freq>minfreq,]



regdiffs <- function(impact, c)
{
  l4 <- loutcomeboot(individualoutcomefunction = step4812cure, simoutput = impact, c = c)
  # cured first round overall:
  cure1 <-t(array(unlist(lapply(l4, function(y) round(c(mean(apply(y[,1,],2,mean)),
                                          sd(apply(y[,1,], 2, mean))), 3)))[c(1:2,7:8)],dim=c(2,2)))
  cure1diff <- cbind(round(c(mean(apply(l4$novelpantb[,1,],2,mean) - apply(l4$baseline[,1,],2,mean)),
                                                     sd(apply(l4$novelpantb[,1,],2,mean) - apply(l4$baseline[,1,],2,mean))), 3),
                    round(c(mean((apply(l4$novelpantb[,1,],2,mean) - apply(l4$baseline[,1,],2,mean))/apply(l4$baseline[,1,],2,mean)),
                            sd((apply(l4$novelpantb[,1,],2,mean) - apply(l4$baseline[,1,],2,mean))/apply(l4$baseline[,1,],2,mean))), 3))
  
  li <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                     simoutput=impact, c=c)
  tbtime <- t(array(unlist(lapply(li, function(x) round(c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean))),3)))[c(1:2,7:8)],dim=c(2,2)))
  tbtimediff <- cbind(round(c(mean(apply(li$novelpantb[,1,],2,mean) - apply(li$baseline[,1,],2,mean)),
                        sd(apply(li$novelpantb[,1,],2,mean) - apply(li$baseline[,1,],2,mean))), 3),
                      round(c(mean((apply(li$novelpantb[,1,],2,mean) - apply(li$baseline[,1,],2,mean))/apply(li$baseline[,1,],2,mean)),
                              sd((apply(li$novelpantb[,1,],2,mean) - apply(li$baseline[,1,],2,mean))/apply(li$baseline[,1,],2,mean))), 3) )

  
    tbtimerr <- t(array(unlist(lapply(li, function(x) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                                              sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum))),3)))[c(1:2,7:8)],dim=c(2,2)))
 
   tbtimediffrr <- cbind(round(c(mean(apply(li$novelpantb[,1,]*(li$novelpantb[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(li$novelpantb[,which(patientvars=="RIF")+2,]==1,2,sum) - 
                                        apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(li$novelrrx[,which(patientvars=="RIF")+2,]==1,2,sum)),
                              sd(apply(li$novelpantb[,1,]*(li$novelpantb[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(li$novelpantb[,which(patientvars=="RIF")+2,]==1,2,sum) -
                                   apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(li$novelrrx[,which(patientvars=="RIF")+2,]==1,2,sum))), 3))
                      
  return(rbind(cure1, cure1diff, tbtime, tbtimediff , tbtimerr, tbtimediffrr)[c(1,5,2,6,3,4,7,8, 9,10),])
}
tresist <- cbind(regdiffs(impact, c), regdiffs(impact, highFQScohort), regdiffs(impact, highFQRcohort), regdiffs(impact, highFQZRcohort))
write.table(tresist, file=paste0("tresist",date,".csv"), sep = ",")
save(tresist, file=paste0("tresist.",date,".Rdata"))
tresist[c(1,3,5,6,8),] <- round(tresist[c(1,3,5,6,8),]*100,1)
tresist[c(2,4,7),] <- round(tresist[c(2,4,7),],2)
tresist2 <- cbind(paste0(tresist[,1], " \u00b1 ",tresist[,2]), paste0(tresist[,3], " \u00b1 ",tresist[,4]), paste0(tresist[,5], " \u00b1 ",tresist[,6]), paste0(tresist[,7], " \u00b1 ",tresist[,8]))
Encoding(tresist2)<-"UTF-8"
tresist2[c(1,3,5,6,8),] <- paste0(tresist2[c(1,3,5,6,8),],"%")
write.table(tresist2, file=paste0("tresist2",date,".csv"), sep = ",")


#
# for multiple scenarios, bar plots, :
yl <- loutcomeboot(function(x) x[c("eventtime","TBstate","RIF"),], impact, c, desiredsize=1e4)
zl <- lapply(yl, function(y) 
  abind(y[,seq(1,52,by=3),], y[,seq(2,53,by=3),], y[,seq(3,54,by=3),], 
        along = 4, new.names = list("patient"=1:1e4, "step"=1:18, "copy"=1:50, "a"=c("time","state","RIF"))))
zlp <- lapply(zl, function(z) aperm(z, perm = c(1,3,2,4)))
rm(zl)
save(yl, zlp, file=paste0("fig2statestatus.",date,".Rdata"))
rm(yl)

maxtime <- 30; times <- seq(0,maxtime,by=0.5)

statestatus <- function(z1, times=0:maxtime,rifs = c(1,0))
{ return(t(array(do.call('rbind', lapply(times, function(t) 
  tabulate(unlist( z1[,,"state"][
    cbind(1:dim(z1)[[1]], c(apply(z1[,,"time"], c(1), function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)))))][z1[,1,"RIF"] %in% rifs]), nbins = length(statetypes))
)), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))))
}


statestatusarray <- function(z2, times=0:30, rifs=c(0,1))
{
  aaply(z2, 2, function(z) statestatus(z,times, rifs))
}

fullstatenames <- c(" \nTB not yet diagnosed\n ", "[shouldn't show up]", "On treatment", " \nOn treatment &\nacquired resistance",
                    "Treatment failed", "Will relapse", "Cured", "Deceased", "Relapsed TB","Lost before (re)treatment")
rbind(unlist(statetypes), fullstatenames)
plotorder <- c(1,10,5,9,6,3,4,7,8)
pdf("Fig2.pdf", width=12, height=8)
layout(matrix(c(1,3,2,4,5,5), ncol=3), widths = c(5,5,2))
par(mar=c(3,3,1,5), oma=c(2,2,2,2))
require(RColorBrewer)
colors <- c("gray", palette(brewer.pal(8, name = "Accent")), "black")
colors <- c("gray", palette(brewer.pal(8, name = "Accent")), "black")
downshift <- 0.5
mapply( function(z, title)
{ 
  meanstates <- apply(statestatusarray(z[,1:5,,], times, rifs=1), c(2,3), mean)
  b <- barplot(xlab = "", ylab="", yaxt='n', xaxt='n', border = NA, space=0,
          meanstates[plotorder,], 
          col=colors, beside=F)
  mtext(title, side=3, font=2,line=0.5, cex=1)
  axis(side = 2, at = seq(0,sum(z[,1:5,1,"RIF"])/5,length=6),labels = seq(0,100,by=20), cex.axis=0.8, las=2)
  axis(side = 1, at = b[seq(1,length(b),length.out = maxtime/3+1)],labels = seq(0,maxtime,by=3), cex.axis=0.8)
  label <- meanstates[plotorder,ncol(meanstates)]/sum(z[,1,1,"RIF"]==1)*100 > 2
  # text(x = (b[length(b)]+c(0,3*ncol(meanstates)/maxtime,rep(0,sum(label)-2))), 
  text(x = rep(b[length(b)], sum(label)), 
       y = (cumsum(meanstates[plotorder,ncol(meanstates)])- meanstates[plotorder,ncol(meanstates)]/2)[label]-downshift,pos=4, xpd=NA,
       labels = paste0(round(meanstates[plotorder,ncol(meanstates)]/sum(z[,1,1,"RIF"]==1)*100,1)[label], "%"), cex=1)
  # segments(b[length(b)],(cumsum(meanstates[plotorder,ncol(meanstates)])[2]- meanstates[plotorder[2],ncol(meanstates)]/2) -1,
  #          b[length(b)]+3*ncol(meanstates)/maxtime,(cumsum(meanstates[plotorder,ncol(meanstates)])[2]- meanstates[plotorder[2],ncol(meanstates)]/2) -1, xpd=NA)
    }, 
          zlp, c("Current practice", "Novel regimen for RIF-R", "Novel RIF-R regimen + expanded RIF DST", "Novel regimen for all"))
mtext("Months elapsed since TB onset", 1, outer=T, cex=1, adj = .35)
mtext("% of RIF-R TB cohort", 2, outer=T, cex=1, line=-0.5)
# mtext("Status of RIF-R cohort over time", side=3, font=2, outer=T, line=0.5)
plot.new()
legend("center", legend = rev(fullstatenames[plotorder]), fill=rev(colors[1:9]), cex=1.2, xpd=NA,title = expression(bold("Status")))
dev.off()


# sensitivity analysis: delays
rm(impact)
load(file=paste0("delays.",date,".Rdata"))

li <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T,
                   simoutput=delays, c=c)
tbtime <- t(array(unlist(lapply(li, function(x) round(c(mean(apply(x[,1,],2,mean)), sd(apply(x[,1,],2,mean))),3)))[c(1:2,7:8)],dim=c(2,2)))
tbtimediff <- cbind(round(c(mean(apply(li$novelpantb[,1,],2,mean) - apply(li$novelrrx[,1,],2,mean)),
                            sd(apply(li$novelpantb[,1,],2,mean) - apply(li$novelrrx[,1,],2,mean))), 3),
                    round(c(mean((apply(li$novelpantb[,1,],2,mean) - apply(li$novelrrx[,1,],2,mean))/apply(li$novelrrx[,1,],2,mean)),
                            sd((apply(li$novelpantb[,1,],2,mean) - apply(li$novelrrx[,1,],2,mean))/apply(li$novelrrx[,1,],2,mean))), 3) )

tbtimerr <- t(array(unlist(lapply(li, function(y) round(c(mean(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                                          sd(apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum))),3)))[c(1:2,7:8)],dim=c(2,2)))

tbtimediffrr <- cbind(round(c(mean((apply(li$novelpantb[,1,]*(li$novelpantb[,which(patientvars=="RIF")+2,]==1),2,sum) - 
                                      apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum))/
                                     apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum)), 
                              sd((apply(li$novelpantb[,1,]*(li$novelpantb[,which(patientvars=="RIF")+2,]==1),2,sum) - 
                                    apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum))/
                                   apply(li$novelrrx[,1,]*(li$novelrrx[,which(patientvars=="RIF")+2,]==1),2,sum)) ), 3))


save(tbtime, tbtimediff, tbtimerr, tbtimediffrr, file=paste0("delayresults.",data,".Rdata"))

rm(delays)

# sensitivity analysis: BPAMZ efficacy
load(paste0("lowerefficacy.",date,".Rdata"))
tloeff <- tableSreg(impact=lowerefficacy, c=c, save=F)
write.table(tloeff, file=paste0("tloeff.",date,".csv"), sep = ",")
save(tloeff, file=paste0("tloeff.",date,".Rdata"))
tloeff[1:16,] <- round(tloeff[1:16,]*100,1)
tloeff[17:24,] <- round(tloeff[17:24,],2)
tloeff2 <- cbind(paste0(tloeff[,1], " \u00b1 ",tloeff[,2])) #, paste0(tloeff[,3], " \u00b1 ",tloeff[,4]), paste0(tloeff[,5], " \u00b1 ",tloeff[,6]))
tloeff2[1:12,] <- paste0(tloeff2[1:12,],"%")
load(paste0("tSEA2.",date,".Rdata"))
tloeff2 <- cbind(tSEA2[,4],tloeff2)
Encoding(tloeff2)<-"UTF-8"
write.table(tloeff2, file=paste0("tloeff2.",date,".csv"), sep = ",")
rm(lowerefficacy)



#############################################
# outcomes for moxi resistant and Xpert XDR

# let's tally total unsuccessful treatments (fail + relapse + death during treatment),
# and tally infectious years (starting from TB onset),
# for each of FQ-R RR and FQ-R RS (each weighted according to their prevalence in an overall TB cohort of 1e5),
# for each of FQ DST (simultaneous) and no,
# and for SEA plus a hypothetical setting with FQ-R in RS on the high end based on pakistan and Mumbai data (Assuming same prev PZA-R).
#( for this last, just need to define a new freq vector, order it like the original one, and use the same courses otherwise)
# and in the text, will report the increase in cure rates amogn RIF-R versus __ in RIF-S when FQ DST is added, 
# but that given RIF R only comprise __% of the overall cohort, the overall cure rate increases from __ to __ with stepwise and __ with simult.
# The number of Xpert XDR tests per additional cure is __ for stepwise and __ for simult. 
# Also give all results for higher prev setting. 

# moxi-r fractions
sum(c$Freq[c$MOXI==1])/sum(c$Freq)
sum(c$Freq[c$MOXI==1])/sum(c$Freq)*1e5
sum(c$Freq[c$MOXI==1&c$RIF==1])/sum(c$Freq)
sum(c$Freq[c$MOXI==1&c$RIF==0])/sum(c$Freq)
sum(c$Freq[c$MOXI==1&c$RIF==1])/sum(c$Freq[c$RIF==1])
sum(c$Freq[c$MOXI==1&c$RIF==0])/sum(c$Freq[c$RIF==0])

sum(c$Freq[c$MOXI==1&c$RIF==1])/sum(c$Freq[c$MOXI==0&c$RIF==1])/
sum(c$Freq[c$MOXI==1&c$RIF==0])/sum(c$Freq[c$MOXI==0&c$RIF==0])

sum(c$Freq[c$PZA==1&c$RIF==1])/sum(c$Freq[c$PZA==0&c$RIF==1])/
  sum(c$Freq[c$PZA==1&c$RIF==0])/sum(c$Freq[c$PZA==0&c$RIF==0])

sum(c$Freq[c$MOXI==1&c$RIF==1&c$PZA==0&c$partialmoxi==1])/sum(c$Freq[c$MOXI==1&c$RIF==1])
sum(c$Freq[c$MOXI==1&c$RIF==1&c$PZA==1&c$partialmoxi==1])/sum(c$Freq[c$MOXI==1&c$RIF==1])


# MOXI-R do worse:
load("lsuccess1.",date,".Rdata")
unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),
                                             sd(apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))), 3)))
unlist(lapply(lsuccess1, function(y) round(c(mean(apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum)),
                                             sd(apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))), 3)))
rm(l4)

# outcomes with DST approaches:
load(paste0("dst.",date,".Rdata"))

# cures with fq dst
l4dst <-  loutcomeboot(individualoutcomefunction = step4812cure, simoutput = dst, c = c)
save(l4dst, file=paste0("l4dst.",date,".Rdata"))

dstcure1RS <- lapply(l4dst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==0)*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply((y[,which(patientvars=="RIF")+3,]==0)*(y[,which(patientvars=="MOXI")+3,]==1),2,sum))
dstcure1RR <- lapply(l4dst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==1)*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply((y[,which(patientvars=="RIF")+3,]==1)*(y[,which(patientvars=="MOXI")+3,]==1),2,sum))
dstcure1all <- lapply(l4dst, function(y) apply(y[,1,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply((y[,which(patientvars=="MOXI")+3,]==1),2,sum))

lapply(dstcure1RS, mean); lapply(dstcure1RS, sd)
lapply(dstcure1RR, mean); lapply(dstcure1RR, sd)

rm(l4dst)

# unsuccessful treatments first round
unsuccessfulrx <- function(patiententry) # without, then with, prerxltfu; and append rif and moxi status, and partialmoxi and pza
  { return(c(patiententry["eventtype",3]==eventtypes$treatmentstart & 
             patiententry["TBstate",4] %in% c(statetypes$failed, statetypes$pendingrelapse, statetypes$deceased ),
        (patiententry["eventtype",2]==eventtypes$TBdiagnosis & 
           patiententry["TBstate",4] %in% c(statetypes$failed, statetypes$pendingrelapse, statetypes$deceased, 
                                            statetypes$pretreatmentlost ) ),
        patiententry["RIF",1],
        patiententry["MOXI",1],
        patiententry["partialmoxi",1],
        patiententry["PZA",1]) ) }
  
ludst <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = c)
save(ludst, file=paste0("ludst.",date,".Rdata"))
# load(file=paste0("ludst.",date,".Rdata"))
dstu1RS <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all <- lapply(ludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))

rm(ludst)

lidst <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                      simoutput=dst, c=c)
save(lidst, file=paste0("lidst.",date,".Rdata"))
# load(file=paste0("lidst.",date,".Rdata"))
dsttimeRS <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeRR <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeall <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))

rm(lidst)

# will boxplot these: 
# colors/shading for DST none vs step vs sim side by side, larger x axis divisions for rs vs rr vs overall, 
# and panels for u and i
require(RColorBrewer)
pdf("Fig3.pdf", width=6, height=10)
layout(c(1,2,3), heights = c(6,6,1))
# par(mfcol=c(2,1))
par(mar=c(4,4,1,1), oma=rep(0.5,4), lwd=0.3, cex=1)
cols <- brewer.pal(3, name = "Accent")
b <- boxplot(list("No FQ DST" = dstu1RS$noxxdr, "FQ DST in RIF-R" = dstu1RS$stepxxdr, "FQ DST for all" = dstu1RS$fullxxdr,
             "No FQ DST" = dstu1RR$noxxdr, "FQ DST in RIF-R" = dstu1RR$stepxxdr, "FQ DST for all" = dstu1RR$fullxxdr,
             "No FQ DST" = dstu1all$noxxdr, "FQ DST in RIF-R" = dstu1all$stepxxdr, "FQ DST for all" = dstu1all$fullxxdr),
        at=c(1,2,3,6,7,8, 11,12,13), 
        col=cols, notch=FALSE,
        ylab = "Unsuccessful treatments of MOXI-R TB",
        ylim=c(0,max(unlist(dstu1all))*1.05),
        xaxt ='n')
text(x=c(1,2,3, 6,7,8, 11,12,13, 16,17,18), y=b$stats[5,]-max(b$stats[5,])/100, pos=3, 
     labels = paste0(round(100*b$stats[3,]/rep(c(sum(c$Freq[c$MOXI==1&c$RIF==0]), sum(c$Freq[c$MOXI==1&c$RIF==1]), sum(c$Freq[c$MOXI==1])), each=3)/1e5, 1),"%"), cex=0.8)

axis(1, at = c(2,7,12), labels = c("Among RIF-S,\nMOXI-R TB", "Among RIF-R,\nMOXI-R TB", "Among all\nMOXI-R TB"), 
  lwd.ticks = F)
legend("topleft", "A", bty="n", cex = 1.5) 

b2 <- boxplot(list("No FQ DST" = dsttimeRS$noxxdr, "FQ DST in RIF-R" = dsttimeRS$stepxxdr, "FQ DST for all" = dsttimeRS$fullxxdr,
             "No FQ DST" = dsttimeRR$noxxdr, "FQ DST in RIF-R" = dsttimeRR$stepxxdr, "FQ DST for all" = dsttimeRR$fullxxdr,
             "No FQ DST" = dsttimeall$noxxdr, "FQ DST in RIF-R" = dsttimeall$stepxxdr, "FQ DST for all" = dsttimeall$fullxxdr),
        at=c(1,2,3,6,7,8, 11,12,13), 
        col=cols, notch=FALSE,
        ylab = "Person-months of active MOXI-R TB",
        ylim=c(0,max(unlist(dsttimeall))),
        xaxt ='n')
axis(1, at = c(2,7,12), labels = c("Among RIF-S,\nMOXI-R TB", "Among RIF-R,\nMOXI-R TB", "Among all\nMOXI-R TB"), 
     lwd.ticks = F)
legend("topleft", "B", bty="n", cex = 1.5) 
par(mar=c(1,1,1,1))
plot.new()
# par(mar=c(1,1,1,1))
legend("center", legend = c("None", "RIF-R", "All"), 
                 pch=15, col=cols, title = expression(bold("FQ DST performed for:")), xpd=NA, ncol=3, pt.cex = 3)
dev.off()


# sensitivity analysis plot, higher moxi and/or pza resistance:

ludsthighFQS <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highFQScohort)
dstu1allhighFQS <- lapply(ludsthighFQS, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
rm(ludsthighFQS)
ludsthighFQR <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highFQRcohort)
dstu1allhighFQR <- lapply(ludsthighFQR, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
rm(ludsthighFQR)
ludsthighFQZR <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highFQZRcohort)
dstu1allhighFQZR <- lapply(ludsthighFQZR, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
rm(ludsthighFQZR)

pdf("Fig3S.pdf", width=11, height=8)
layout(c(1,2), heights = c(6,1))
par(mar=c(4,4,1,1), oma=rep(0.5,4), lwd=0.3, cex=1)
cols <- brewer.pal(3, name = "Accent")
b <- boxplot(list("No FQ DST" = dstu1all$noxxdr, "FQ DST in RIF-R" = dstu1all$stepxxdr, "FQ DST for all" = dstu1all$fullxxdr,
                  "No FQ DST" = dstu1allhighFQS$noxxdr, "FQ DST in RIF-R" = dstu1allhighFQS$stepxxdr, "FQ DST for all" = dstu1allhighFQS$fullxxdr,
                  "No FQ DST" = dstu1allhighFQR$noxxdr, "FQ DST in RIF-R" = dstu1allhighFQR$stepxxdr, "FQ DST for all" = dstu1allhighFQR$fullxxdr,
                  "No FQ DST" = dstu1allhighFQZR$noxxdr, "FQ DST in RIF-R" = dstu1allhighFQZR$stepxxdr, "FQ DST for all" = dstu1allhighFQZR$fullxxdr),
             at=c(1,2,3, 5,6,7, 9,10,11, 13,14,15), 
             col=cols, notch=FALSE,
             ylab = "Unsuccessful treatments of MOXI-R TB",
             ylim=c(0,max(unlist(dstu1allhighFQZR))),
             xaxt ='n')
axis(1, at = c(2,6,10,14), labels = c("With original\nMOXI-R and PZA-R\n\n", "Increased MOXI-R\namong RIF-S\n\n", "Increased MOXI-R\namong RIF-S and RIF-R\n\n", "Increased MOXI-R\namong RIF-S and RIF-R,\nand increased PZA-R\namong RIF-R"), 
     lwd.ticks = F, padj=0.5)

plot.new()
# par(mar=c(1,1,1,1))
legend("center", legend = c("None", "RIF-R", "All"), 
       pch=15, col=cols, title = expression(bold("FQ DST performed for:")), xpd=NA, ncol=3, pt.cex = 3)
dev.off()


# Sensitivity alaysis, DST in SAf

load("dst.20190226.Rdata")
ludstSAf <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = SAfcohort)
save(ludstSAf, file=paste0("ludstSAf.",date,".Rdata"))
# load(file=paste0("ludst.",date,".Rdata"))
dstu1RS <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))

rm(ludstSAf)

lidstSAf <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                      simoutput=dst, c=SAfcohort)
save(lidstSAf, file=paste0("lidstSAf.",date,".Rdata"))
# load(file=paste0("lidst.",date,".Rdata"))
dsttimeRS <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeRR <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeall <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))

rm(lidstSAf)
rm(dst)

pdf("Fig3S2.pdf", width=6, height=10)
layout(c(1,2,3), heights = c(6,6,1))
# par(mfcol=c(2,1))
par(mar=c(4,4,1,1), oma=rep(0.5,4), lwd=0.3, cex=1)
cols <- brewer.pal(3, name = "Accent")
b <- boxplot(list("No FQ DST" = dstu1RS$noxxdr, "FQ DST in RIF-R" = dstu1RS$stepxxdr, "FQ DST for all" = dstu1RS$fullxxdr,
                  "No FQ DST" = dstu1RR$noxxdr, "FQ DST in RIF-R" = dstu1RR$stepxxdr, "FQ DST for all" = dstu1RR$fullxxdr,
                  "No FQ DST" = dstu1all$noxxdr, "FQ DST in RIF-R" = dstu1all$stepxxdr, "FQ DST for all" = dstu1all$fullxxdr),
             at=c(1,2,3,6,7,8, 11,12,13), 
             col=cols, notch=FALSE,
             ylab = "Unsuccessful treatments of MOXI-R TB",
             ylim=c(0,max(unlist(dstu1all))*1.05),
             xaxt ='n')
text(x=c(1,2,3, 6,7,8, 11,12,13, 16,17,18), y=b$stats[5,]-max(b$stats[5,])/100, pos=3, 
     labels = paste0(round(100*b$stats[3,]/rep(c(sum(c$Freq[c$MOXI==1&c$RIF==0]), sum(c$Freq[c$MOXI==1&c$RIF==1]), sum(c$Freq[c$MOXI==1])), each=3)/1e5, 1),"%"), cex=0.8)

axis(1, at = c(2,7,12), labels = c("Among RIF-S,\nMOXI-R TB", "Among RIF-R,\nMOXI-R TB", "Among all\nMOXI-R TB"), 
     lwd.ticks = F)
legend("topleft", "A", bty="n", cex = 1.5) 

b2 <- boxplot(list("No FQ DST" = dsttimeRS$noxxdr, "FQ DST in RIF-R" = dsttimeRS$stepxxdr, "FQ DST for all" = dsttimeRS$fullxxdr,
                   "No FQ DST" = dsttimeRR$noxxdr, "FQ DST in RIF-R" = dsttimeRR$stepxxdr, "FQ DST for all" = dsttimeRR$fullxxdr,
                   "No FQ DST" = dsttimeall$noxxdr, "FQ DST in RIF-R" = dsttimeall$stepxxdr, "FQ DST for all" = dsttimeall$fullxxdr),
              at=c(1,2,3,6,7,8, 11,12,13), 
              col=cols, notch=FALSE,
              ylab = "Person-months of active MOXI-R TB",
              ylim=c(0,max(unlist(dsttimeall))),
              xaxt ='n')
axis(1, at = c(2,7,12), labels = c("Among RIF-S,\nMOXI-R TB", "Among RIF-R,\nMOXI-R TB", "Among all\nMOXI-R TB"), 
     lwd.ticks = F)
legend("topleft", "B", bty="n", cex = 1.5) 
par(mar=c(1,1,1,1))
plot.new()
# par(mar=c(1,1,1,1))
legend("center", legend = c("None", "RIF-R", "All"), 
       pch=15, col=cols, title = expression(bold("FQ DST performed for:")), xpd=NA, ncol=3, pt.cex = 3)
dev.off()



# rr vs rs moxi-r:
sum((c$RIF==0)*c$Freq)/sum(c$Freq)
sum((c$RIF==0)*c$MOXI*c$Freq)/sum(c$MOXI*c$Freq)

# fold-diffs in unsuccessful outcomes averted, simultaneous vs stepwise
(dstu1all$fullxxdr-dstu1all$noxxdr)
(dstu1all$stepxxdr-dstu1all$noxxdr)
mean(dstu1all$stepxxdr-dstu1all$noxxdr)
median(dstu1all$stepxxdr-dstu1all$noxxdr)
mean((dstu1all$fullxxdr-dstu1all$noxxdr)/mean(dstu1all$stepxxdr-dstu1all$noxxdr))

(dsttimeall$fullxxdr-dsttimeall$noxxdr)
(dsttimeall$stepxxdr-dsttimeall$noxxdr)
mean((dsttimeall$fullxxdr-dsttimeall$noxxdr)/mean(dsttimeall$stepxxdr-dsttimeall$noxxdr))
mean((dsttimeall$fullxxdr-dsttimeall$noxxdr)/mean(dsttimeall$stepxxdr-dsttimeall$noxxdr))


# dsts performed:
dsts <- function(patiententry) # xperts, then xxdrs, each first round then overall; and append rif and moxi status, and partialmoxi and pza
{ return(c(floor(patiententry["DSTs",4]),
           round(10*(patiententry["DSTs",4]%%1)),
           floor(patiententry["DSTs",ncol(patiententry)]),
           round(10*(patiententry["DSTs",ncol(patiententry)]%%1)),
           patiententry["RIF",1],
           patiententry["MOXI",1],
           patiententry["partialmoxi",1],
           patiententry["PZA",1]) ) 
  }

dstdst <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=c)
save(dstdst, file=paste0("dstdst.",date,".Rdata"))

# total increased FQ DSTs from universal approach:
apply(dstdst$fullxxdr[,4,] - dstdst$stepxxdr[,4,] , 2, sum)
mean(apply(dstdst$fullxxdr[,2,] - dstdst$stepxxdr[,2,] , 2, sum)); sd(apply(dstdst$fullxxdr[,2,] - dstdst$stepxxdr[,2,] , 2, sum))

mean(apply(dstdst$fullxxdr[,2,] - dstdst$stepxxdr[,2,], 2, sum)/apply(dstdst$stepxxdr[,2,] - dstdst$noxxdr[,2,], 2, sum))
sd(apply(dstdst$fullxxdr[,2,] - dstdst$stepxxdr[,2,], 2, sum)/apply(dstdst$stepxxdr[,2,] - dstdst$noxxdr[,2,], 2, sum))


# dsts per added successful treatment (first treatment round only), comparing step to none and full to step:
(apply(dstdst$stepxxdr[,2,],2,sum) - apply(dstdst$fullxxdr[,2,],2,sum))/
  (dstu1all$fullxxdr - dstu1all$stepxxdr)
mean((apply(dstdst$stepxxdr[,2,],2,sum) - apply(dstdst$fullxxdr[,2,],2,sum))/
  (dstu1all$fullxxdr - dstu1all$stepxxdr));
sd((apply(dstdst$stepxxdr[,2,],2,sum) - apply(dstdst$fullxxdr[,2,],2,sum))/
     (dstu1all$fullxxdr - dstu1all$stepxxdr))

(apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum))/
  (dstu1all$stepxxdr - dstu1all$noxxdr)
median((apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum))/
         (dstu1all$stepxxdr - dstu1all$noxxdr))
apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum)
mean(dstu1all$stepxxdr - dstu1all$noxxdr)
median(dstu1all$stepxxdr - dstu1all$noxxdr)
mean((apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum))/
    median(dstu1all$stepxxdr - dstu1all$noxxdr))
mean((apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum))/
       mean(dstu1all$stepxxdr - dstu1all$noxxdr))
sd((apply(dstdst$noxxdr[,2,],2,sum) - apply(dstdst$stepxxdr[,2,],2,sum))/
       mean(dstu1all$stepxxdr - dstu1all$noxxdr))


rm(dstdst)



#####################################
# Acquired resistance

# load(paste0("dst.",date,".Rdata"))
dst <- dst[c("noxxdr","fullxxdr")]
load(paste0("impact.",date,".Rdata"))

combined <- list("baseline"=impact$baseline, "novelrrx"=impact$novelrrx, "noxxdr"=dst$noxxdr, "fullxxdr"=dst$fullxxdr)
save(combined, file=paste0("combined.",date,".Rdata"))
rm(impact)
rm(dst)

# # estimated reduction in infectious time and force of infection (not accounting for infectiousness during treatment, changes in a person's infectiousness over time, etc): 
# l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*10, simoutput=combined, c=c)
# # overall
# lapply(l, function(x) apply((l$baseline[,1,]-x[,1,]), 2, sum)/apply(l$baseline[,1,],2,sum))

# ## Tracking acquired resistance
# adrevents <- function(patiententry, drugs=c("RIF","INH","MOXI","BDQ","PA")) # total events, which drugs acquired, and initial resistance
# {
#   return(c(
#     sum( 
#       patiententry["eventtype",]==eventtypes$treatmentstart & 
#       patiententry["TBstate",]==statetypes$treating_adr & 
#       colSums(patiententry[drugs,] != cbind(patiententry[drugs,1],patiententry[drugs,1:(ncol(patiententry)-1)]))>=1 ),
#     rowSums(patiententry[drugs,] != cbind(patiententry[drugs,1],patiententry[drugs,1:(ncol(patiententry)-1)])),
#     patiententry[drugs,1] ))
# }
# 
# ladr <- loutcomeboot(individualoutcomefunction = adrevents, simoutput = combined, c = c)
# 
# lapply(ladr, function(y) apply(y[,1,],2,sum)) #total events
# lapply(ladr, function(y) apply(y[,1,],2,sum)/1e5) #total event rate
# lapply(ladr, function(y) apply(y[,1,]*y[,1+5+1,],2,sum)/apply(y[,1+5+1,],2,sum)) #rate of events among initial RIF-R
# lapply(ladr, function(y) apply(y[,1,]*y[,1+5+3,],2,sum)/apply(y[,1+5+3,],2,sum)) #rate of events among initial moxi-R
# lapply(ladr, function(y) apply(y[,1,]*y[,1+5+4,],2,sum)/apply(y[,1+5+1,],2,sum)) #rate of events among initial BDQ-R
# lapply(ladr, function(y) apply(y[,1,]*y[,1+1,],2,mean)) #acquisitions of (any) RIF-R
# lapply(ladr, function(y) apply(y[,1,]*y[,1+3,],2,mean)) #acquisitions of (any) MOXI-R ??
# lapply(ladr, function(y) apply(y[,1,]*y[,1+4,],2,mean)) #acquisitions of (any) BDQ-R ??
# 
#   

# infectious time with drug resistance: 
infectioustime <- function(patiententry, 
                           resdrugs=list(
                                         "rif"=list(c("RIF"), 1), 
                                         "moxi"=list(c("MOXI"), 1),
                                         "himoxi"=list(c("MOXI","partialmoxi"), c(1,0), "and"),
                                         "bdq"=list(c("BDQ"), 1), 
                                         "pa"=list(c("PA"), 1),
                                         "noveldrug"=list(c("BDQ","PA"), c(1,1), "or"),
                                         "none"=list(c())), 
                           includeinitial=TRUE, cutofftime=12*10) # follow up to 4th diagnosis ie step 13. exclude untreated (states 1 and 2) if includeinitial==F.
{
  elapsedtimes <- patiententry["eventtime",2:13] - patiententry["eventtime",1:12] 
  elapsedtimes[patiententry["eventtime",2:13]>cutofftime] <- cutofftime - (patiententry["eventtime",1:12])[patiententry["eventtime",2:13]>cutofftime]
  elapsedtimes[elapsedtimes<0] <- 0
  
  if (includeinitial) indices <-  patiententry["TBstate",1:12] %in% c(1,2,5,9,10) else 
    indices <-  patiententry["TBstate",1:12] %in% c(5,9,10)
  
  indexlist <- lapply(resdrugs, function(y) {
      if(length(y[[1]])==1) (indices & patiententry[y[[1]],1:12]==y[[2]]) else
      if(length(y[[1]])>1) {if(y[[3]]=="and") indices & apply(patiententry[y[[1]],1:12]==y[[2]], 2, all) else
                            if(y[[3]]=="or") indices & apply(patiententry[y[[1]],1:12]==y[[2]], 2, any) } else
      indices
    })
  
  t <- lapply(indexlist, function(y) sum(elapsedtimes*y))

  return(c(unlist(t), patiententry[patientvars,1]))
}

# names(combined) <- c("baseline", "novelrrx", "noxxdr", "fullxxdr")

lres <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=TRUE, simoutput = combined, c=c)
# lanytb <- loutcomeboot(individualoutcomefunction = infectioustime, resdrugs=list(c()), includeinitial=TRUE, simoutput = combined, c=c)
lrestreated <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = combined, c=c)
# lanytbtreated <- loutcomeboot(individualoutcomefunction = infectioustime, resdrugs=list(c()), includeinitial=FALSE, simoutput = combined, c=c)

save(lres, lrestreated, file=paste0("lres_etc.",date,".Rdata"))

rm(combined)

# note that if n% of the cohort is BDQ-R, then average outcomes will be a weighted average of outcomes in the original cohort, and those in the bdq-r version of the cohort.
# but a weighted average won't capture the variability among that 5%. 
# so if I want to who what happens with 4% BDQ-R (i.e. same as current rif-r), then I should run loutcomeboot with a desiredsize of 4/(100-4)*1e5, using combinedbdq,
# and then do a weighted average of each result with lres etc. Can just use the same c for the freq's, and I'll weight it later.
load(paste0("combinedbdq.",date,".Rdata"))
lbdqres <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=TRUE, simoutput = combinedbdq, c=c, desiredsize = 1e5*4/(100-4))
lbdqrestreated <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = combinedbdq, c=c, desiredsize = 1e5*4/(100-4))
save(lbdqres, lbdqrestreated, file=paste0("lbdqres_etc.",date,".Rdata"))

# rifr_bdq <- mapply(function(a,b) apply(a[,1,],2,mean)*.96 + apply(b[,1,],2,mean)*.04, a=lres, b=lbdqres)
# rifrtreated_bdq <- mapply(function(a,b) mean(apply(a[,1,],2,mean)*.96 + apply(b[,1,],2,mean)*.04), a=lrestreated, b=lbdqrestreated)

# boxplot(cbind(mapply(function(x,y) apply(y[,4,],2,mean), x=1, y=lres), mapply(function(a,b) apply(a[,4,],2,mean)*.96 + apply(b[,4,],2,mean)*.04, a=lres, b=lbdqres)), ylim=c(0,1))

# will mention bdq sensitiivty analysis in text,  plan to put a figure in appendix (with a focux on post-treatment time, and what that means for the next round of trnamission)

rm(combinedbdq)

load("highres.20190226.Rdata")
lhighres <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=TRUE, simoutput = highres, c=c)
lhighrestreated <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = highres, c=c)
save(lhighres, lhighrestreated, file=paste0("lhighres_etc.",date,".Rdata"))

boxplot(lapply(lrestreated, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lhighrestreated, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))

rm(highres)

###### FIGURE 4 #############
load(paste0("lres_etc.",date,".Rdata"))
load(paste0("lhighres_etc.",date,".Rdata"))
names(lhighres) <- names(lhighrestreated) <- c( "baseline", "novelrrx", "noxxdr", "fullxxdr")

a1 <- data.frame(lapply(lrestreated, function(y) apply(y[,1,],2,mean)/apply(y[,7,],2,mean)))
a2 <- data.frame(lapply(lrestreated, function(y) apply(y[,2,],2,mean)/apply(y[,7,],2,mean)))
a3 <- data.frame(lapply(lrestreated, function(y) apply(y[,4,],2,mean)/apply(y[,7,],2,mean)))
a4 <- data.frame(lapply(lrestreated, function(y) apply(y[,5,],2,mean)/apply(y[,7,],2,mean)))
a1$drug <- "Rifampin"
a2$drug <- "Moxifloxacin"
a3$drug <- "Bedaquiline"
a4$drug <- "Pretomanid"
a <- bind_rows(a1,a2,a3,a4)
am <- melt(a, id.vars=c("drug"))
am$params <- "Low ADR"

c1 <- data.frame(lapply(lres, function(y) apply(y[,1,],2,mean)/apply(y[,7,],2,mean)))
c2 <- data.frame(lapply(lres, function(y) apply(y[,2,],2,mean)/apply(y[,7,],2,mean)))
c3 <- data.frame(lapply(lres, function(y) apply(y[,4,],2,mean)/apply(y[,7,],2,mean)))
c4 <- data.frame(lapply(lres, function(y) apply(y[,5,],2,mean)/apply(y[,7,],2,mean)))
c1$drug <- "Rifampin"
c2$drug <- "Moxifloxacin"
c3$drug <- "Bedaquiline"
c4$drug <- "Pretomanid"
c <- bind_rows(c1,c2,c3,c4)
cm <- melt(a, id.vars=c("drug"))
cm$params <- "Low ADR"


b1 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,1,],2,mean)/apply(y[,7,],2,mean)))
b2 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,2,],2,mean)/apply(y[,7,],2,mean)))
b3 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,4,],2,mean)/apply(y[,7,],2,mean)))
b4 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,5,],2,mean)/apply(y[,7,],2,mean)))
b1$drug <- "Rifampin"
b2$drug <- "Moxifloxacin"
b3$drug <- "Bedaquiline"
b4$drug <- "Pretomanid"
b <- bind_rows(b1,b2,b3,b4)
bm <- melt(b, id.vars=c("drug"))
bm$params <- "High ADR"
require(dplyr)
require(data.table)
mm <- bind_rows(am,bm)
colnames(mm) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mm$Parameters <- factor(mm$Parameters, levels=c("Low ADR", "High ADR"))
levels(mm$Parameters) <- c("Low", "High")
# levels(mm$Parameters) <- c("Post-Treatment TB, Low BPaMZ ADR Risk", "Post-Treatment TB, High BPaMZ ADR Risk")

d1 <- data.frame(lapply(lhighres, function(y) apply(y[,1,],2,mean)/apply(y[,7,],2,mean)))
d2 <- data.frame(lapply(lhighres, function(y) apply(y[,2,],2,mean)/apply(y[,7,],2,mean)))
d3 <- data.frame(lapply(lhighres, function(y) apply(y[,4,],2,mean)/apply(y[,7,],2,mean)))
d4 <- data.frame(lapply(lhighres, function(y) apply(y[,5,],2,mean)/apply(y[,7,],2,mean)))
d1$drug <- "Rifampin"
d2$drug <- "Moxifloxacin"
d3$drug <- "Bedaquiline"
d4$drug <- "Pretomanid"
d <- bind_rows(d1,d2,d3,d4)
dm <- melt(d, id.vars=c("drug"))
dm$params <- "High ADR"
require(dplyr)
require(data.table)
mm2 <- bind_rows(cm,dm)
colnames(mm2) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mm2$Parameters <- factor(mm2$Parameters, levels=c("Low ADR", "High ADR"))
levels(mm2$Parameters) <- c("Low", "High")

mm$tbtime <- "Post-treatment time with TB"
mm2$tbtime <- "All time with TB"

mm3 <- rbind(mm,mm2)
mm3$Drug <- factor(mm3$Drug)
mm3$Scenario <- factor(mm3$Scenario, levels=rev(levels(mm3$Scenario)))
levels(mm3$Scenario) <- rev(c("Current Practice",
                          "BPaMZ for RIF-R, no FQ DST",
                          "BPaMZ for all, no FQ DST",
                          "BPaMZ for all, universal FQ DST"))

save(am, bm, cm, dm, mm, mm2, mm3, file="Fig4data.Rdata")


# total time, not proportion:

e1 <- data.frame(lapply(lrestreated, function(y) apply(y[,1,],2,sum)))
e2 <- data.frame(lapply(lrestreated, function(y) apply(y[,2,],2,sum)))
e3 <- data.frame(lapply(lrestreated, function(y) apply(y[,4,],2,sum)))
e4 <- data.frame(lapply(lrestreated, function(y) apply(y[,5,],2,sum)))
e1$drug <- "Rifampin"
e2$drug <- "Moxifloxacin"
e3$drug <- "Bedaquiline"
e4$drug <- "Pretomanid"
e <- bind_rows(e1,e2,e3,e4)
em <- melt(e, id.vars=c("drug"))
em$params <- "Low"

g1 <- data.frame(lapply(lres, function(y) apply(y[,1,],2,sum)))
g2 <- data.frame(lapply(lres, function(y) apply(y[,2,],2,sum)))
g3 <- data.frame(lapply(lres, function(y) apply(y[,4,],2,sum)))
g4 <- data.frame(lapply(lres, function(y) apply(y[,5,],2,sum)))
g1$drug <- "Rifampin"
g2$drug <- "Moxifloxacin"
g3$drug <- "Bedaquiline"
g4$drug <- "Pretomanid"
g <- bind_rows(g1,g2,g3,g4)
gm <- melt(g, id.vars=c("drug"))
gm$params <- "Low"


names(lhighres) <- names(lhighrestreated) <- c( "baseline", "novelrrx", "noxxdr", "fullxxdr")

f1 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,1,],2,sum)))
f2 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,2,],2,sum)))
f3 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,4,],2,sum)))
f4 <- data.frame(lapply(lhighrestreated, function(y) apply(y[,5,],2,sum)))
f1$drug <- "Rifampin"
f2$drug <- "Moxifloxacin"
f3$drug <- "Bedaquiline"
f4$drug <- "Pretomanid"
f <- bind_rows(f1,f2,f3,f4)
fm <- melt(f, id.vars=c("drug"))
fm$params <- "High"
require(dplyr)
require(data.table)
mmt <- bind_rows(em,fm)
colnames(mmt) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mmt$Parameters <- factor(mmt$Parameters, levels=c("Low", "High"))

h1 <- data.frame(lapply(lhighres, function(y) apply(y[,1,],2,sum)))
h2 <- data.frame(lapply(lhighres, function(y) apply(y[,2,],2,sum)))
h3 <- data.frame(lapply(lhighres, function(y) apply(y[,4,],2,sum)))
h4 <- data.frame(lapply(lhighres, function(y) apply(y[,5,],2,sum)))
h1$drug <- "Rifampin"
h2$drug <- "Moxifloxacin"
h3$drug <- "Bedaquiline"
h4$drug <- "Pretomanid"
h <- bind_rows(h1,h2,h3,h4)
hm <- melt(h, id.vars=c("drug"))
hm$params <- "High"
require(dplyr)
require(data.table)
mm2t <- bind_rows(gm,hm)
colnames(mm2t) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mm2t$Parameters <- factor(mm2t$Parameters, levels=c("Low", "High"))

mmt$tbtime <- "Post-treatment time with TB"
mm2t$tbtime <- "All time with TB"

mm3t <- rbind(mmt,mm2t)
mm3t$Drug <- factor(mm3t$Drug)
mm3t$Scenario <- factor(mm3t$Scenario, levels=rev(levels(mm3t$Scenario)))
levels(mm3t$Scenario) <- rev(c("Current Practice",
                              "BPaMZ for RIF-R, no FQ DST",
                              "BPaMZ for all, no FQ DST",
                              "BPaMZ for all, universal FQ DST"))

save(em, fm, gm, hm, mmt, mm2t, mm3t, file="Fig4data_t.Rdata")


require(ggplot2)

load("Fig4data_t.Rdata")
pdf("Fig4fixed.pdf", width=8, height=10)
# better would be total time (person-months in cohort, not %) with dr-tb
ggplot(data = transform(mm3t, Drug=factor(Drug, levels=c("Rifampin","Moxifloxacin","Bedaquiline","Pretomanid"))),
       aes(x=Scenario, y=Infectious_Time)) + 
  geom_col(aes(fill=Parameters), position = position_dodge(), width=0.8) + 
  facet_grid(Drug ~ tbtime) + 
  coord_flip() +
  # ggtitle("Potential for drug resistance transmission") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("BPaMZ implementation scenario") + 
  ylab("Person-months with active TB resistant to the specified drug, per 100,000 TB cases") +
 theme(legend.position="bottom") + guides(fill=guide_legend(ncol=2)) + labs(fill = "Risk of acquiring resistance during BPaMZ:") 
dev.off()

# ggplot(data = transform(mm, Drug=factor(Drug, levels=c("Moxifloxacin","Rifampin","Bedaquiline","Pretomanid"))),
#                         aes(x=Scenario, y=Infectious_Time)) + 
#   geom_col(aes(fill=Parameters), position = position_dodge(0.8), width=0.5) + 
#   facet_wrap( ~ Drug, ncol=2) + 
#   coord_flip() +
#   ggtitle("Potential for drug resistance transmission") + theme(plot.title = element_text(hjust = 0.5)) + 
#   ylab("TB in cohort after first treatment attempt:\nFraction resistant to the specified drug") + 
#   theme(legend.position="bottom") + guides(fill=guide_legend(ncol=2)) + labs(fill = "Fraction of TB\nresistant to drug") + 
#   geom_col(data=transform(mm2, Drug=factor(Drug, levels=c("Moxifloxacin","Rifampin","Bedaquiline","Pretomanid"))),
#                               aes(x=as.numeric(Scenario)+0.1, fill=Parameters), position= position_dodge(0.8), width=0.5) +
#   scale_fill_manual(values=c("black", "darkblue", "gray","lightblue")) 
# 
# dev.off()


### Interpretation and results text:
# RIF-S, MOXI-R TB:
load("combined.20190226.Rdata")
lrsmrtreated <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = combined, c=c)
mapply(function(x,y) apply(x[,1,],2,sum)/apply(y[,1,],2,sum), x=lrsmrtreated$noxxdr, y=lrsmrtreated$baseline)

 # Because most time with active TB occurred prior to diagnosis and treatment, BPaMZ had little effect on total time with active RIF-R or MOXI-R TB. 
# With our initial low estimates of the risk of acquiring resistance to MOXI, BDQ, or PA during BPaMZ treatment, the amount of BDQ and PA resistance remained low even among TB that remained active or recurred after attempted treatment. Some acquisition of BDQ and PA did occur with universal use of BPaMZ, but the combined person-years of BDQ-R and PA-R TB were outweighed by the reduction in time spend with RIF-R TB. 
# With high-end estimates for the risk of acquired resistance, however, BDQ-R and PA-R TB within the cohort could each exceed RIF-R TB after attempted BPaMZ treatment.
# FQ DST reduced preventable (i.e., occurring after TB diagnosis and treatment) time with MOXI-R TB by approximately
univ <- mm3t %>% filter(Drug=="Moxifloxacin" & Parameters=="Low" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, universal FQ DST") %>% select(Infectious_Time)
none <- mm3t %>% filter(Drug=="Moxifloxacin" & Parameters=="Low" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, no FQ DST") %>% select(Infectious_Time)
mean(1-unlist(univ/none)); sd(1-unlist(univ/none))

univ <- mm3t %>% filter(Drug=="Moxifloxacin" & Parameters=="High" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, universal FQ DST") %>% select(Infectious_Time)
none <- mm3t %>% filter(Drug=="Moxifloxacin" & Parameters=="High" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, no FQ DST") %>% select(Infectious_Time)
mean(1-unlist(univ/none)); sd(1-unlist(univ/none))

univ <- mm3t %>% filter(Drug=="Bedaquiline" & Parameters=="High" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, universal FQ DST") %>% select(Infectious_Time)
none <- mm3t %>% filter(Drug=="Bedaquiline" & Parameters=="High" & tbtime=="Post-treatment time with TB" & Scenario=="BPaMZ for all, no FQ DST") %>% select(Infectious_Time)
mean(1-unlist(univ/none)); sd(1-unlist(univ/none))

## In figure legend: With the "Low" acquired resistance parameters, the risk that pan-susceptible TB would acquire resistance to any drug when treated with BPaMZ was comparable to risk of acquiring RIF-R when treated with HRZE. With the "high" parameters, this risk increased to 1% for each of MOXI, BDQ, and PA, and we also used high-end estimates for how much this risk increased when MOXI or PZA resistance was present at baseline.


# Figure for BDQ (appendix):
# same as above, but with bdq (weighted with lres) instead of high adq param.
str(lbdqres); str(lres)
bdqprev <- 1 # change only if adjusting from the initial (4%) prevalence used in freq
i1 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,1,],2,sum) + apply(y[,1,],2,sum), x=lbdqrestreated, y=lrestreated))
i2 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,2,],2,sum) + apply(y[,2,],2,sum), x=lbdqrestreated, y=lrestreated))
i3 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,4,],2,sum) + apply(y[,4,],2,sum), x=lbdqrestreated, y=lrestreated))
i4 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,5,],2,sum) + apply(y[,5,],2,sum), x=lbdqrestreated, y=lrestreated))
i1$drug <- "Rifampin"
i2$drug <- "Moxifloxacin"
i3$drug <- "Bedaquiline"
i4$drug <- "Pretomanid"
i <- bind_rows(i1,i2,i3,i4)
im <- melt(i, id.vars=c("drug"))
im$params <- "High BDQ"
require(dplyr)
require(data.table)
mmb <- bind_rows(em,im)
colnames(mmb) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mmb$Parameters <- factor(mmb$Parameters, levels=c("Low", "High BDQ"))

j1 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,1,],2,sum) + apply(y[,1,],2,sum), x=lbdqres, y=lres))
j2 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,2,],2,sum) + apply(y[,2,],2,sum), x=lbdqres, y=lres))
j3 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,4,],2,sum) + apply(y[,4,],2,sum), x=lbdqres, y=lres))
j4 <- data.frame(mapply(function(x,y) bdqprev*apply(x[,5,],2,sum) + apply(y[,5,],2,sum), x=lbdqres, y=lres))
j1$drug <- "Rifampin"
j2$drug <- "Moxifloxacin"
j3$drug <- "Bedaquiline"
j4$drug <- "Pretomanid"
j <- bind_rows(j1,j2,j3,j4)
jm <- melt(j, id.vars=c("drug"))
jm$params <- "High BDQ"
mm2b <- bind_rows(gm,jm)
colnames(mm2b) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
mm2b$Parameters <- factor(mm2b$Parameters, levels=c("Low", "High BDQ"))

mmb$tbtime <- "Post-treatment time with TB"
mm2b$tbtime <- "All time with TB"

mm3b <- rbind(mmb,mm2b)
mm3b$Drug <- factor(mm3b$Drug)
mm3b$Scenario <- factor(mm3b$Scenario, levels=rev(levels(mm3b$Scenario)))
levels(mm3b$Scenario) <- rev(c("Current Practice",
                               "BPaMZ for RIF-R, no FQ DST",
                               "BPaMZ for all, no FQ DST",
                               "BPaMZ for all, universal FQ DST"))

save(em, gm, im, jm, mmb, mm2b, mm3b, file="Fig4data_b.Rdata")

pdf("Fig4S.pdf", width=8, height=12)
# better would be total time (person-months in cohort, not %) with dr-tb
ggplot(data = transform(mm3b, Drug=factor(Drug, levels=c("Rifampin","Moxifloxacin","Bedaquiline","Pretomanid"))),
       aes(x=Scenario, y=Infectious_Time)) +
  geom_col(aes(fill=Parameters), position = position_dodge(), width=0.8) +
  facet_grid(Drug ~ tbtime) +
  coord_flip() +
  # ggtitle("Potential for drug resistance transmission") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Person-time with  active TB resistant to the specified drug") +
  theme(legend.position="bottom") + guides(fill=guide_legend(ncol=2)) + labs(fill = "BDQ resistance prevalence at baseline:")
dev.off()


# 
# # tally up time alive: # need to cut all scenarios off at the same point, say 3 years - no, need a longer time window, see below:
# 
# l <- loutcomeboot(individualoutcomefunction =time.in.state, states=1:7, cutofftime=12*followyears, carryforward=T, 
#                   simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
# lapply(l, function(x) apply(x[,1,],2,summary))
# # why is survival staying same or going down with novelrr for RIF==1 and cutofftime=36?? (same for MXOI==0) need to recheck with larger iniital sim **
# ## one possibility: I model low mortality during treatment, so longer treatment courses allow less opportunity for TB death. Esp with a cutoff at 3 years, there's no opportunity for death after unsuccessful conventional MDR-TB treament, but some do die of TB in the period after novel RR TB treatement (even though ultimately their risk is lower)
# # but I'm seeing this to some extent even for very long cutofftimes -- novelrr doesn't reduce mortality. maybe bc most get hrze, so most mortality is before mdr treament starts, and of those who get to mdr treatment, there's still the efect that longer regimen defers relapse/failure mortality
# # for MOXI==1, the lack of benefit of pantb is expected. 
# #(there was also a problem with remaining diagnosed at the last step, which I've hopefuly fixed 12/30)
# # for RIF==0 at baseline, there's a slight mortality increase with novel RR in 20181218 -- just chance? **
# 
# # incremental months of life per patient: 
# lapply(l, function(x) apply(x[,1,] - l$baseline[,1,],2,mean))
# # plot years of life gained within first [3] years:
# par(mfrow=c(1,1), mar=c(4,4,3,1), oma=c(1,1,1,1))
# boxplot( lapply(l, function(x) apply(x[,1,] - l$baseline[,1,],2,mean)), main=paste0("Incremental months of life gained per patient,\nover ",followyears," years after TB onset"))
# 
# 
# # time deceased i.e. YLL (compare to time alive, should see inverse)
# # but don't trust this because of tail of life expectancy
# l <- loutcomeboot(individualoutcomefunction =time.in.state, states=8, cutofftime=12*followyears, carryforward=T, 
#                   simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
# lapply(l, function(x) apply(x[,1,],2,mean))
# lapply(l, function(x) apply(x[,1,]- l$baseline[,1,],2,mean)/12)
# 
# 
# # cures, as opposed to death or still on treatment at xx months:
# t <- 48; 
# l <- loutcomeboot(individualoutcomefunction =still.in.state, states=7, time=t, 
#                   simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
# boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main=paste0("Proportion with cure (vs ongoing TB or death) ",t/12," years after TB onset"), col=colors, ylim=c(0,1))
# 
# 
# # deaths within 5 years of TB onset:
# t <- 120
# l <- loutcomeboot(individualoutcomefunction =still.in.state, states=8, time=t, 
#                   simoutput=impact, c=c, include=(c$RIF==1&c$RxHist==0), copies = 3, desiredsize = 1e5)
# boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main=paste0("Mortality within ",t/12," years of TB onset, proportion of cohort"), col=colors, ylim=c(0,1))
# ## ** need to figure out why no mortality benefit for novelRR (here, RxHist and RIF do worse, and for new RIF there's no difference)
# ## are outcomes wrong? make.recurrence.matrix()[c("MDR, FQ-S", "MDR, FQ-R", "(ZE)","BPaMZ", "BPaM", "BPaZ"),c("6","18")] looks okay
# ## we saw above that we don't reduce YLL either.
# ## don't expect to reduce time to RR diagnosis (which in many cases will be long, with one or more (ZE) treatments first and taking 2 yrs or more),
# ## and once diagnosed as RR, expect TB mortality risk among ~10% after 18 months, vs among ~3% after 6 months --> 3% + 10%*3% over ~24 months,
# ## so by 24 months we should be seeing novelrr mortality benefit.
# ## it's possible the problem was different maxtimes between the scenarios, since I was setting mortality limit for each internally based on the max time recorded, 
# ## meaning longer maxtime and more opportunity for natural mortality in the one with longer time steps -- but that should cause more mortality for baseline. 
# # with larger samples it's clear the mortality reduction is small but ~2%. which would make sense of >50% of mortality is before 1st treatment, another 50%+ is after improper treatment, and after MDR treatment the new regimen reduces relapse/failure from 10% to 2% (10% reduction in recurrence prevents death in 5% of this 25%?)
# # 5x greater impact for novelrr. will want to do sensitivity analysis around xpert coverage params, make sure this responds as expected. 
# 
# # time to cure, if cured (where cure is assumed to happen when treatment stops):
# par(mar=c(3,3,3,1))
# l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
#                    simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
# lapply(l, function(x) apply(x[,1,],2,summary))
# boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRR-TB patients only")
# l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
#                    simoutput=impact, c=c, include=(c$RIF==0), copies = 3, desiredsize = 1e5)
# boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRS-TB patients only")
# 
# 
# # and for dsT:
# l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, 
#                   simoutput=dst, c=c, include=(c$MOXI==1), copies = 3, desiredsize = 1e4)
# lapply(l, function(x) apply(x[,1,],2,summary))
# lapply(l, function(x) apply(x[,1,],2,mean)) 
# boxplot(lapply(l, function(x) x[,1,1]), main="Total infectious time, RR-TB cases")
# 

# reduction <- lapply(impact, function(x)
#   (apply(outcomes$baseline, 4, function(y) (time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"])) - 
#      apply(x, 4, function(y) (time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"]))) /
#     apply(outcomes$baseline, 4, function(y) time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"]))
# b <- boxplot(reduction, main="Proportion of total infectious time averted, compared to baseline", col=colors[2:6])
# text(1:5,unlist(lapply(reduction,median))+0.008, paste0(round(unlist(lapply(reduction,median))*100,1),"%"))


# ld30 <- loutcomeboot(individualoutcomefunction = still.in.state, states=c(1,2,3,4,5,6,7,9), t=30, simoutput=impact, c=c, include=(c$RIF==1), copies=100, desiredsize = 1e5)
# aliveat30m <- rbind(unlist(lapply(ld30, function(l) mean(apply(l[,1,],1,mean)))),
#                     unlist(lapply(ld30, function(l) sd(apply(l[,1,],1,mean)))))
# 
# lu30 <- loutcomeboot(individualoutcomefunction = still.in.state, states=c(1), t=30, simoutput=impact, c=c, include=(c$RIF==1), copies=100, desiredsize = 1e5)
# undxat30m <- rbind(unlist(lapply(lu30, function(l) mean(apply(l[,1,],1,mean)))),
#                    unlist(lapply(lu30, function(l) sd(apply(l[,1,],1,mean)))))


