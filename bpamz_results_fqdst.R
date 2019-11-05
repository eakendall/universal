# setwd("C:/Users/Administrator/OneDrive - Johns Hopkins University/Research/universal regimen/bpamz/fq dst manuscript")

date <- "20190711" # fa dst runs, assumes 4mo as bpamz duration for all, except in sens analysis
setting <- "SAf"
require(plyr)
source("bpamz_cohort.R")
source("bpamz_result_utils.R")

load(paste0("cohortfreqs.",date,".Rdata"), verbose = T)
minfreq <- 1/1e8
SAfc <- SAfcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]
SEAc <- SEAcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]


# some bookkeeping
# making some alternative cohorts that we'll use for sensitivity analyses:

highrmparams <- params;
highrmparams["MOXI-R-any_in_RIF-S"] <- increaseodds(params["MOXI-R-any_in_RIF-S"], 1/5)
highrmparams["MOXI-R-any_in_RIF-R"] <- increaseodds(params["MOXI-R-any_in_RIF-R"], 1/5)
highrmparams["MOXI-R-highlevel_in_RIF-S"] <- increaseodds(params["MOXI-R-highlevel_in_RIF-S"], 1/5)
highrmparams["MOXI-R-highlevel_in_RIF-R"] <- increaseodds(params["MOXI-R-highlevel_in_RIF-R"], 1/5)
highrmparams["RIF-R_in_New"] <- increaseodds(params["RIF-R_in_New"], 1/10)
highrmparams["RIF-R_in_Retreatment"] <- increaseodds(params["RIF-R_in_Retreatment"], 1/10)
highrmcohort <- cohort.probs(params=highrmparams, patientvars = patientvars)
highrmc <- highrmcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]


highrmzparams <- highrmparams
highrmzparams["PZA-R_in_RIF-S"] <- increaseodds(highrmparams["PZA-R_in_RIF-S"], 1/5)
highrmzparams["PZA-R_in_RIF-R"] <- increaseodds(highrmparams["PZA-R_in_RIF-R"], 1/5)
highrmzcohort <- cohort.probs(params=highrmzparams, patientvars = patientvars)
highrmzc <- highrmzcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]


highrzparams <- highrmzparams
highrzparams["MOXI-R-any_in_RIF-S"] <- params["MOXI-R-any_in_RIF-S"]
highrzparams["MOXI-R-any_in_RIF-R"] <- params["MOXI-R-any_in_RIF-R"]
highrzparams["MOXI-R-highlevel_in_RIF-S"] <- params["MOXI-R-highlevel_in_RIF-S"]
highrzparams["MOXI-R-highlevel_in_RIF-R"] <-  params["MOXI-R-highlevel_in_RIF-R"]
highrzcohort <- cohort.probs(params=highrzparams, patientvars = patientvars)
highrzc <- highrzcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]



lowmzparams <- params; 
lowmzparams["PZA-R_in_RIF-R"] <- params["PZA-R_in_RIF-S"]
lowmzparams["MOXI-R-any_in_RIF-R"] <- params["MOXI-R-any_in_RIF-S"]
lowmzparams["MOXI-R-highlevel_in_RIF-R"] <- params["MOXI-R-highlevel_in_RIF-S"]
lowmzcohort <- cohort.probs(params=lowmzparams, patientvars = patientvars)
lowmzc <- lowmzcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]
weighted.mean(lowmzc$RIF, w = lowmzc$Freq)
weighted.mean(lowmzc$PZA, w = lowmzc$Freq)
weighted.mean(lowmzc$PZA[lowmzc$RIF==1], w = lowmzc$Freq[lowmzc$RIF==1])
weighted.mean(lowmzc$PZA[lowmzc$RIF==0], w = lowmzc$Freq[lowmzc$RIF==0])
weighted.mean(lowmzc$MOXI, w = lowmzc$Freq)
weighted.mean(lowmzc$MOXI[lowmzc$RIF==1], w = lowmzc$Freq[lowmzc$RIF==1])
weighted.mean(lowmzc$MOXI[lowmzc$RIF==0], w = lowmzc$Freq[lowmzc$RIF==0])

norzcorparams <- params; 
norzcorparams["PZA-R_in_RIF-R"] <- norzcorparams["PZA-R_in_RIF-S"] <- 
  (params["PZA-R_in_RIF-S"]*
     (params["Retreatment_in_All"]*(1-params["RIF-R_in_Retreatment"]) + (1-params["Retreatment_in_All"])*(1-params["RIF-R_in_New"])) + 
   params["PZA-R_in_RIF-R"]*
     (params["Retreatment_in_All"]*(params["RIF-R_in_Retreatment"]) + (1-params["Retreatment_in_All"])*(params["RIF-R_in_New"])))
norzcorcohort <- cohort.probs(params=norzcorparams, patientvars = patientvars)
norzcorc <- norzcorcohort[SEAcohort$Freq>minfreq|SAfcohort$Freq>minfreq,]


# make table showing outcomes, summary of full generation and recurrent matrix and adr matrix which will be included in appendix. 
rec <- make.recurrence.matrix(params)[c("HR(ZE)","R(ZE)","BPaMZ","BPamZ","BPaM","BPaZ", "BPam","BPa","MDR, FQ-S", "MDR, FQ-R"),c("4","6","18")]
# (rec <- round(rbind((rec[1,]*sum(c$Freq*(1-c$RIF)*(1-c$INH))+rec[2,]*sum(c$Freq*(1-c$RIF)*c$INH))/sum(c$Freq*(1-c$RIF)),
#       rec[3:nrow(rec),]),3)
# )
# to check: all displayed entries should be less than this, else need a min of the two:
round(rowSums(make.adr.matrix(params)[c("HR(ZE)","R(ZE)","BPaMZ","BPamZ","BPaM","BPaZ", "BPam","BPa","MDR, FQ-S", "MDR, FQ-R"),], na.rm=T)*100, 1)
write.table(round(rec,3), "recurrencetable.csv", sep=",", row.names = TRUE, col.names = TRUE)

# details for bpamz supplement:
drugs
make.active.regimen.matrix()[,1:(length(drugs)+length(regimens))] # describe, and refer to published code
a <- make.recurrence.matrix(params); a[c(1:4,8:nrow(a)),c("9","12","18")] <- NA; round(a,3)
write.table(a, file="recurrencematrix.csv", sep=",")
a <- make.adr.matrix(params); a[a==0]<-NA; round(a,3)
write.table(a, file="adrmatrix.csv", sep=",")


# make bpamz outcome table
tableSreg <- function(impact, c, save)
{
# #  proportion cured after 1st and 2nd round:
  l4 <- loutcomeboot(individualoutcomefunction = step4812cure, simoutput = impact, c = c)
  if(save) save(l4, file=paste0("l4.",date,".Rdata"))
  # cured first round: overall, rif-s, rif-r, and fq-r:
  cure1 <- rbind(
    unlist(lapply(l4, function(y) round(quantile((apply(y[,1,],2,mean)),c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(l4, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum)),
                                                 c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(l4, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum)),
                                                 c(0.5, 0.25,0.75)), 3)) ), 
    unlist(lapply(l4, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum)),
                                                 c(0.5, 0.25,0.75)), 3))) )
  rownames(cure1) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure1
  
  # cured second round: overall, rif-s, rif-r, and fq-r:
  cure2 <- rbind(
    unlist(lapply(l4, function(y) round(quantile((apply(y[,2,],2,mean)),
                                                 c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(l4, function(y) round(quantile((apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==0,2,sum)),
                                                 c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(l4, function(y) round(quantile((apply(y[,2,]*(y[,which(patientvars=="RIF")+3,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+3,]==1,2,sum)),
                                        c(0.5, 0.25,0.75)), 3)) ), 
    unlist(lapply(l4, function(y) round(quantile((apply(y[,2,]*(y[,which(patientvars=="MOXI")+3,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+3,]==1,2,sum)),
                                        c(0.5, 0.25,0.75)), 3))) )
  rownames(cure2) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure2

  #cured by N=36 months
  lc36 <- loutcomeboot(individualoutcomefunction = still.in.state, states=7, t=36, simoutput=impact, c=c)
  if(save) save(lc36, file=paste0("lc36.",date,".Rdata"))
  # load(paste0("lc36.",date,".Rdata"))
  
  cure36m <- rbind(
    unlist(lapply(lc36, function(y) round(quantile((apply(y[,1,],2,mean)),
                                                   c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lc36, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                                   c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lc36, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                                   c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lc36, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                                   c(0.5, 0.25,0.75)), 3))) )
  rownames(cure36m) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  cure36m
  
  # cure first round if treated
  lsuccess1 <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=impact, c=c)
  if(save) save(lsuccess1, file=paste0("lsuccess1.",date,".Rdata"))

  success <- rbind(
    unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)),
                                                        c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==0),2,sum)/apply(y[,1,]*(y[,3,]==0),2,sum)),
                                                        c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1),2,sum)),
                                                        c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                        c(0.5, 0.25,0.75)), 3))) )
  rownames(success) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  success
  
  # total TB/infectious time:
  li <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                     simoutput=impact, c=c)
  if(save) save(li, file=paste0("li.",date,".Rdata"))
  # load(file=paste0("li.",date,".Rdata"))
  tbtime <- rbind(
    unlist(lapply(li, function(x) round(quantile((apply(x[,1,],2,mean)), c(0.5, 0.25,0.75)), 3))), 
    unlist(lapply(li, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                          c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(li, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                          c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(li, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                          c(0.5, 0.25,0.75)), 3))))

  # time on treatment (of some kind) 
  followyears <- 20
  ltall <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(statetypes$treating, statetypes$treating_adr), cutofftime=12*followyears, carryforward=T, 
                        simoutput=impact, c=c)
  if(save) save(ltall, file=paste0("ltall.",date,".Rdata"))
  # load(file=paste0("ltall.",date,".Rdata"))
  
  rxtime <- rbind(
    unlist(lapply(ltall, function(x) round(quantile((apply(x[,1,],2,mean)), c(0.5, 0.25,0.75)), 3))), 
    unlist(lapply(ltall, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==0,2,sum)),
                                                    c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(ltall, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1),2,sum)/apply(y[,which(patientvars=="RIF")+2,]==1,2,sum)),
                                                    c(0.5, 0.25,0.75)), 3))),
    unlist(lapply(ltall, function(y) round(quantile((apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum)/apply(y[,which(patientvars=="MOXI")+2,]==1,2,sum)),
                                                    c(0.5, 0.25,0.75)), 3)))
  )
  rownames(rxtime) <- c("All", "RIF-S", "RIF-R", "MOXI-R")
  rxtime


  return(rbind(cure1, cure2, cure36m, success, tbtime, rxtime))
}




######################################

# for second manuscript:

# cohort compositions:
cs <- cbind(
c(
  "% RR",
  "% MR",
  "% MR RR",
  "% MR RS",
  "% ZR",
  "% ZR MR RR",
  "% ZR MR RS"),
paste0(round(100*c(
    sum(subset(SAfc, RIF==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, MOXI==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, RIF==1 & MOXI==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, RIF==0 & MOXI==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, PZA==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(SAfc$Freq),
    sum(subset(SAfc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(SAfc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(SEAc, RIF==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, MOXI==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, RIF==1 & MOXI==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, RIF==0 & MOXI==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, PZA==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(SEAc$Freq),
  sum(subset(SEAc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(SEAc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(highrmzc, RIF==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, MOXI==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, RIF==1 & MOXI==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, RIF==0 & MOXI==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, PZA==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(highrmzc$Freq),
  sum(subset(highrmzc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(highrmzc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(highrmc, RIF==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, MOXI==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, RIF==1 & MOXI==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, RIF==0 & MOXI==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, PZA==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(highrmc$Freq),
  sum(subset(highrmc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(highrmc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(highrzc, RIF==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, MOXI==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, RIF==1 & MOXI==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, RIF==0 & MOXI==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, PZA==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(highrzc$Freq),
  sum(subset(highrzc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(highrzc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(lowmzc, RIF==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, MOXI==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, RIF==1 & MOXI==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, RIF==0 & MOXI==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, PZA==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(lowmzc$Freq),
  sum(subset(lowmzc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(lowmzc$Freq)),1),"%"),
paste0(round(100*c(
  sum(subset(norzcorc, RIF==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, MOXI==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, RIF==1 & MOXI==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, RIF==0 & MOXI==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, PZA==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, RIF==0 & MOXI==1 & PZA==1)$Freq)/sum(norzcorc$Freq),
  sum(subset(norzcorc, RIF==1 & MOXI==1 & PZA==1 )$Freq)/sum(norzcorc$Freq)),2),"%")

)
cs
write.csv(cs, file="cs.csv")

sum(subset(SAfc, RIF==1)$Freq)/sum(SAfc$Freq)
sum(subset(SEAc, RIF==1)$Freq)/sum(SEAc$Freq)
sum(subset(SAfc, RIF==1 & MOXI==1)$Freq)/sum(SAfc$Freq)
sum(subset(SEAc, RIF==1 & MOXI==1)$Freq)/sum(SEAc$Freq)
sum(subset(SEAc, RIF==0 & MOXI==1)$Freq)/sum(subset(SEAc, RIF==0)$Freq)
sum(subset(SAfc, PZA==1)$Freq)/sum(SAfc$Freq)
sum(subset(SEAc, PZA==1)$Freq)/sum(SEAc$Freq)
sum(subset(SAfc, MOXI==1)$Freq)/sum(SAfc$Freq)
sum(subset(SEAc, MOXI==1)$Freq)/sum(SEAc$Freq)

sum(subset(SEAc, RIF==1 & MOXI==1)$Freq)/sum(subset(SEAc, RIF==1)$Freq)/
  (sum(subset(SEAc, RIF==0 & MOXI==1)$Freq)/sum(subset(SEAc, RIF==0)$Freq))
sum(subset(SAfc, RIF==1 & MOXI==1)$Freq)/sum(subset(SAfc, RIF==1)$Freq)/
  (sum(subset(SAfc, RIF==0 & MOXI==1)$Freq)/sum(subset(SAfc, RIF==0)$Freq))

# proportion of moxi-r (and pza-r) occurring in rs:
sum(subset(SEAc, RIF==0 & MOXI==1)$Freq)/sum(subset(SEAc, MOXI==1)$Freq)
sum(subset(SEAc, RIF==0 & PZA==1)$Freq)/sum(subset(SEAc, PZA==1)$Freq)

sum(subset(SAfc, RIF==0 & MOXI==1)$Freq)/sum(subset(SAfc, MOXI==1)$Freq)
sum(subset(SAfc, RIF==0 & PZA==1)$Freq)/sum(subset(SAfc, PZA==1)$Freq)

sum(subset(SEAc, partialmoxi==0 & MOXI==1)$Freq)/sum(subset(SEAc, MOXI==1)$Freq)
sum(subset(SAfc, partialmoxi==0 & MOXI==1)$Freq)/sum(subset(SAfc, MOXI==1)$Freq)




# outcomes with different FQ DST appraoches:

load(paste0("dst.",date,".Rdata"))

# dst outcomes in table format for each setting:
tSAfdst <- tableSreg(impact = dst, c = SAfc, save = T)
save(tSAfdst, file=paste0("tSAfdst.",date,".Rdata"))
temp <- tSAfdst
temp[1:16,] <- round(temp[1:16,]*100,1)
temp[17:24,] <- round(temp[17:24,],2)
temp[1:16, c(1,3,4,6,7,9,10,12)] <- paste0(temp[1:16, c(1,3,4,6,7,9,10,12)],"%")
tSAf2dst <- cbind(paste0(temp[,1], " (",temp[,2],"-",temp[,3],")"), 
                  paste0(temp[,4], " (",temp[,5],"-",temp[,6],")"), 
                  paste0(temp[,7], " (",temp[,8],"-",temp[,9],")"),
                  paste0(temp[,10], " (",temp[,11],"-",temp[,12],")"))
colnames(tSAf2dst) <- c("None",
                        "Stepwise (For RIF-R TB)", "Universal", "")
write.table(tSAf2dst, file=paste0("tSAf2dst.",date,".csv"), sep = ",")
save(tSAf2dst, file=paste0("tSAf2dst.",date,".Rdata"))
tSAf2dst


tSEAdst <- tableSreg(impact = dst, c = SEAc, save = F)
save(tSEAdst, file=paste0("tSEAdst.",date,".Rdata"))
temp <- tSEAdst
temp[1:16,] <- round(temp[1:16,]*100,1)
temp[17:24,] <- round(temp[17:24,],2)
temp[1:16, c(1,3,4,6,7,9,10,12)] <- paste0(temp[1:16, c(1,3,4,6,7,9,10,12)],"%")
tSEA2dst <- cbind(paste0(temp[,1], " (",temp[,2],"-",temp[,3],")"), 
                  paste0(temp[,4], " (",temp[,5],"-",temp[,6],")"), 
                  paste0(temp[,7], " (",temp[,8],"-",temp[,9],")"),
                  paste0(temp[,10], " (",temp[,11],"-",temp[,12],")"))
colnames(tSEA2dst) <- c("None",
                        "Stepwise (For RIF-R TB)", "Universal", "")
write.table(tSEA2dst, file=paste0("tSEA2dst.",date,".csv"), sep = ",")
save(tSEA2dst, file=paste0("tSEA2dst.",date,".Rdata"))
tSEA2dst



#############################################
# outcomes for moxi resistant and Xpert XDR

lsuccess1SAf <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=SAfc)
lsuccess1SEA <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=SEAc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
save(lsuccess1SEA, lsuccess1SAf, file=paste0("lsuccess.dst.",date,".Rdata"))

lsuccess1 <- lsuccess1SAf # can toggle between these
lsuccess1 <- lsuccess1SEA

# treated, cured, rif, moxi, partialmoxi, pza
# all rR
unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1),2,sum)),
                                             c(0.25,0.5,0.75)), 3)))
# RR MR
unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),
                                             c(0.25,0.5,0.75)), 3)))
# RS MR
unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum)),
                                             sd(apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))), 3)))
# all RS
unlist(lapply(lsuccess1, function(y) round(quantile((apply(y[,2,]*(y[,3,]==0),2,sum)/apply(y[,1,]*(y[,3,]==0),2,sum)),
                                             sd(apply(y[,2,]*(y[,3,]==0),2,sum)/apply(y[,1,]*(y[,3,]==0),2,sum))), 3)))

rm(lsuccess1)


# cure among MFX-R RR TB detected by Xpert - use only for none vs stepwise comparison
lknownRR <- loutcomeboot(individualoutcomefunction = knownRRsuccess, simoutput = dst, c=SAfc)
lapply(lknownRR, function(y) quantile((apply(y[,1,],2,sum)/apply(y[,2,],2,sum)), c(0.25,0.5,0.75)))


# outcomes with DST approaches:

  
ludst <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = SEAc)
save(ludst, file=paste0("ludst.",date,".Rdata"))
# number of poor outcomes in moxi-R:
dstu1RS <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all <- lapply(ludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))

rm(ludst)

# ludst will be first round treatments, but li will be infectious time over all four cycles
lidst <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                      simoutput=dst, c=SEAc)
save(lidst, file=paste0("lidst.",date,".Rdata"))
# load(file=paste0("lidst.",date,".Rdata"))
dsttimeRS <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeRR <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeall <-  lapply(lidst, function(y) apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))

rm(lidst)

# and SAf dst results for figure:
ludstSAf <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = SAfc)
save(ludstSAf, file=paste0("ludstSAf.",date,".Rdata"))
# load(file=paste0("ludstSAf.",date,".Rdata"))
dstu1RSSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RRSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1allSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))

rm(ludstSAf)

lidstSAf <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9,10), cutofftime=12*20, carryforward=T, 
                         simoutput=dst, c=SAfc)
save(lidstSAf, file=paste0("lidstSAf.",date,".Rdata"))
# load(file=paste0("lidstSAf.",date,".Rdata"))
dsttimeRSSAf <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==0)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeRRSAf <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="RIF")+2,]==1)*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))
dsttimeallSAf <-  lapply(lidstSAf, function(y) apply(y[,1,]*(y[,which(patientvars=="MOXI")+2,]==1),2,sum))

rm(lidstSAf)

# for later sensitivity analyses:
lsuccess1RMZ <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=highrmzc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
lsuccess1lowMZ <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=lowmzc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
lsuccess1RM <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=highrmc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
lsuccess1RZ <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=highrzc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
save(lsuccess1RMZ, lsuccess1lowMZ, file=paste0("lsuccess.dst.sensis.",date,".Rdata"))


ludstRMZ <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highrmzc)
dstu1RS_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
save(ludstRMZ, file=paste0("ludstRMZ.",date,".Rdata"))
rm(ludstRMZ)

ludstlowMZ <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = lowmzc)
dstu1RS_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
save(ludstlowMZ, file=paste0("ludstlowMZ.",date,".Rdata"))
rm(ludstlowMZ)

ludstRM <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highrmc)
dstu1RS_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
save(ludstRM, file=paste0("ludstRM.",date,".Rdata"))
rm(ludstRM)

ludstRZ <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highrzc)
dstu1RS_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
save(ludstRZ, file=paste0("ludstRZ.",date,".Rdata"))
rm(ludstRZ)


# proportion of unsuccessful treatments that are in moxi-R:
load(paste0("ludstSAf.",date,".Rdata"))
lapply(lapply(lapply(ludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum)/apply(y[,1,],2,sum)), quantile, c(0.5,0.25,0.75)),round,3)
load(paste0("ludst.",date,".Rdata"))
lapply(lapply(lapply(ludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum)/apply(y[,1,],2,sum)), quantile, c(0.5,0.25,0.75)),round,3)
# and proportion of population:
weighted.mean(SAfc$MOXI, w = SAfc$Freq)
weighted.mean(SEAc$MOXI, w = SEAc$Freq)

# # text results: cure rates, MFX RR/S TB:
# lsuccess1SAf <- loutcomeboot(individualoutcomefunction = treatmentsuccess, simoutput=dst, c=SAfc) # 1 event 2 state 3 rif 4 moxi 5 partial moxi
# unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum)),
#                                             c(0.25,0.5,0.75), , 3)))
# unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),
# c(0.25,0.5,0.75), , 3)))


# outcomes among MFX-R RR-TB that was missed by Xpert: 
# need a function that looks for success among those receiving a particular regimen in the previous step 
# -- particularly, for our text we want those who were on bpamz4 (and had RIF res), otherwise function can look like treatmentsuccess **
lmissedRR <- loutcomeboot(individualoutcomefunction = missedRRsuccess, simoutput = dst, c=SAfc)
lapply(lmissedRR, function(y) apply(y[,1,],2,sum)/apply(y[,2,],2,sum))
lapply(lmissedRR, function(y) quantile((apply(y[,1,],2,sum)/apply(y[,2,],2,sum)), c(0.25,0.5,0.75)))
lmissedRR_SEA <- loutcomeboot(individualoutcomefunction = missedRRsuccess, simoutput = dst, c=SEAc)
lapply(lmissedRR_SEA, function(y) quantile((apply(y[,1,],2,sum)/apply(y[,2,],2,sum)), c(0.25,0.5,0.75)))

# # cure among all MFX-R TB: see tables
# 
# # cure among all MFX-R RR TB: as above -- (still using SAf lsuccess1 here)
# unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),
#                                              c(0.25,0.5,0.75), 3)))


# single combined SEA and SAf DST figure:
require(RColorBrewer)
# pdf("dstFig2cures.pdf", width=11, height=8)
tiff("dstFig2cures.tif", width=9, height=7, units='in', res=300)
layout(array(c(1,3,2,3), dim=c(2,2)), heights = c(6,1))
par(mar=c(4,3,1,1), oma=rep(c(0.5,2,0.5,0.5)), lwd=0.3, cex=1)
# cols <- brewer.pal(4, name = "Dark2")
cols <- brewer.pal(6, name = "Greys")[c(1,3,5)]
labcols <- brewer.pal(9, name = "Greys")[6:8]
b3 <- boxplot(list("No FQ DST" = dstu1RSSAf$noxxdr, "FQ DST in RIF-R" = dstu1RSSAf$stepxxdr, "FQ DST for all" = dstu1RSSAf$fullxxdr,
                   "No FQ DST" = dstu1RRSAf$noxxdr, "FQ DST in RIF-R" = dstu1RRSAf$stepxxdr, "FQ DST for all" = dstu1RRSAf$fullxxdr,
                   "No FQ DST" = dstu1allSAf$noxxdr, "FQ DST in RIF-R" = dstu1allSAf$stepxxdr, "FQ DST for all" = dstu1allSAf$fullxxdr),
              at=c(1,2,3,6,7,8, 11,12,13), 
              col=cols, notch=FALSE, outline = F, range=0, pars=list(lwd=1),
              # ylab = "Unsuccessful treatments of MOXI-R TB,\namong a cohort of 100,000 individuals with TB",
              ylim=c(0,max(unlist(dstu1all))*1.05),
              xaxt ='n',
              main="South Africa")
text(x=c(1,2,3, 6,7,8, 11,12,13)-0.5 + c(-0.2,-0.1,0.1), y=b3$stats[5,]-max(b3$stats[5,])/100+30 + c(20,2.5,-15), pos=4, 
     labels = paste0(round(100*b3$stats[3,]/rep(c(sum(SAfc$Freq[SAfc$MOXI==1&SAfc$RIF==0]), sum(SAfc$Freq[SAfc$MOXI==1&SAfc$RIF==1]), sum(SAfc$Freq[SAfc$MOXI==1])), each=3)/1e5, 1),"%"), cex=0.9,
     col=labcols, xpd=NA, font=2)
# segments(x0 = c(1,2,3)+0.6, y0 = b3$stats[5,1:3]-max(b3$stats[5,])/100+60, 
#          x1 = c(2,2,2)+0.7, y1 = b3$stats[5,c(2,2,2)]-max(b3$stats[5,])/100+140, 
#          col=cols[1:3])
# text(x=2.5, y=b3$stats[5,c(2,2,2)]-max(b3$stats[5,])/100+150, labels = "Median proportion unsuccessfully\ntreated within subset",
#      pos = 4, srt=90, cex=0.9)

axis(1, at = c(2,7,12), labels = c(paste0("Rifampin-\nsusceptible\n(N = ",round(sum(SAfc$Freq[SAfc$RIF==0&SAfc$MOXI==1])*100000),")"),
                                   paste0("Rifampin-\nresistant\n(N = ", round(sum(SAfc$Freq[SAfc$RIF==1&SAfc$MOXI==1])*100000),")"),
                                   paste0("All\n(N = ",round(sum(SAfc$Freq[SAfc$MOXI==1])*100000),")\n ")), 
     lwd.ticks = F, line = 1, col = "white")

# legend("topleft", "A", bty="n", cex = 1) 

b <- boxplot(list("No FQ DST" = dstu1RS$noxxdr, "FQ DST in RIF-R" = dstu1RS$stepxxdr, "FQ DST for all" = dstu1RS$fullxxdr,
                  "No FQ DST" = dstu1RR$noxxdr, "FQ DST in RIF-R" = dstu1RR$stepxxdr, "FQ DST for all" = dstu1RR$fullxxdr,
                  "No FQ DST" = dstu1all$noxxdr, "FQ DST in RIF-R" = dstu1all$stepxxdr, "FQ DST for all" = dstu1all$fullxxdr),
             at=c(1,2,3,6,7,8, 11,12,13), 
             col=cols, notch=FALSE, outline=F, range=0,pars=list(lwd=1),
             # ylab = "Unsuccessful treatments of MOXI-R TB",
             ylim=c(0,max(unlist(dstu1all))*1.05),
             xaxt ='n',
             main="Southeast Asia")
text(x=c(1,2,3, 6,7,8, 11,12,13)-0.5 + c(-0.2,-0.1,0.1), y=b$stats[5,]-max(b$stats[5,])/100+30 + c(20,2.5,-15), pos=4, 
     labels = paste0(round(100*b$stats[3,]/rep(c(sum(SEAc$Freq[SEAc$MOXI==1&SEAc$RIF==0]), sum(SEAc$Freq[SEAc$MOXI==1&SEAc$RIF==1]), sum(SEAc$Freq[SEAc$MOXI==1])), each=3)/1e5, 1),"%",
                     # c("", rep(" unsuccessful",2),rep("",6)),
                     ""), 
     cex=0.9, font=2, 
     col=labcols, xpd=NA)

axis(1, at = c(2,7,12), labels = c(paste0("Rifampin-\nsusceptible\n(N = ",round(sum(SEAc$Freq[SAfc$RIF==0&SEAc$MOXI==1])*100000),")"),
                                   paste0("Rifampin-\nresistant\n(N = ", round(sum(SEAc$Freq[SAfc$RIF==1&SEAc$MOXI==1])*100000),")"),
                                   paste0("All\n(N = ",round(sum(SEAc$Freq[SEAc$MOXI==1])*100000),")\n ")), 
     lwd.ticks = F, line = 1, col = "white")
# legend("topleft", "B", bty="n", cex = 1) 

mtext("Subset of moxifloxacin-resistant TB",side=1, outer=T, line=-5.5, font=2)
mtext("       Absolute number of unsuccessful treatments\n       of moxifloxacin-resistant TB",side=2, outer=T, line=0, font=2)

par(mar=c(1,1,1,1))
plot.new()
legend('top', 
       inset=0.3, 
       legend = c("None", "RR-TB", "All TB"), 
       pch=22, 
       col=c(1,1,1), pt.bg=cols, title = expression(bold("With fluoroquinolone DST  performed for:")), xpd=NA, ncol=3, pt.cex = 3,
      border="black")

dev.off()
  


pdf("dstFig2Stime.pdf", width=11, height=8)
layout(array(c(1,3,2,3), dim=c(2,2)), heights = c(6,1))
par(mar=c(4,3,1,1), oma=rep(c(0.5,2,0.5,0.5)), lwd=0.3, cex=1)
cols <- brewer.pal(4, name = "Dark2")
par(mar=c(4,3,1,1), oma=rep(c(0.5,2,0.5,0.5)), lwd=0.3, cex=1)
b4 <- boxplot(list("No FQ DST" = dsttimeRSSAf$noxxdr, "FQ DST in RIF-R" = dsttimeRSSAf$stepxxdr, "FQ DST for all" = dsttimeRSSAf$fullxxdr,
                   "No FQ DST" = dsttimeRRSAf$noxxdr, "FQ DST in RIF-R" = dsttimeRRSAf$stepxxdr, "FQ DST for all" = dsttimeRRSAf$fullxxdr,
                   "No FQ DST" = dsttimeallSAf$noxxdr, "FQ DST in RIF-R" = dsttimeallSAf$stepxxdr, "FQ DST for all" = dsttimeallSAf$fullxxdr),
              at=c(1,2,3,6,7,8, 11,12,13), 
              col=cols[1:3], notch=FALSE,
              # ylab = "Person-months of active MOXI-R TB",
              ylim=c(0,max(unlist(dsttimeall))),
              main="South Africa",
              xaxt ='n')
axis(1, at = c(2,7,12), labels = c("Rifampin-\nsusceptible", "Rifampin-\nresistant", "All\n "), 
     lwd.ticks = F)

b2 <- boxplot(list("No FQ DST" = dsttimeRS$noxxdr, "FQ DST in RIF-R" = dsttimeRS$stepxxdr, "FQ DST for all" = dsttimeRS$fullxxdr,
                   "No FQ DST" = dsttimeRR$noxxdr, "FQ DST in RIF-R" = dsttimeRR$stepxxdr, "FQ DST for all" = dsttimeRR$fullxxdr,
                   "No FQ DST" = dsttimeall$noxxdr, "FQ DST in RIF-R" = dsttimeall$stepxxdr, "FQ DST for all" = dsttimeall$fullxxdr),
              at=c(1,2,3,6,7,8, 11,12,13), 
              col=cols[1:3], notch=FALSE,
              # ylab = "Person-months of active MOXI-R TB",
              ylim=c(0,max(unlist(dsttimeall))),
              main="Southeast Asia",
              xaxt ='n')
axis(1, at = c(2,7,12), labels = c("Rifampin-\nsusceptible", "Rifampin-\nresistant", "All\n "), 
     lwd.ticks = F)
# legend("topleft", "B", bty="n", cex = 1) 

mtext("Subset of moxifloxacin-resistant TB",side=1, outer=T, line=-7.5, font=2)

mtext("Person-months of active moxifloxacin-resistant TB,\namong a cohort of 100,000 individuals with TB",side=2, outer=T, line=0, font=2)



par(mar=c(1,1,1,1))
plot.new()
legend("center", legend = c("None", "RR-RB", "All TB"), 
       pch=15, col=cols[1:3], title = expression(bold("With fluoroquinolone DST performed for:")), xpd=NA, ncol=3, pt.cex = 3)
dev.off()


# more text results:

# contribution  of RS to MOXI-R bpamz poor outcomes
quantile(dstu1RSSAf$noxxdr/dstu1allSAf$noxxdr, c(0.25,0.5,0.75))
quantile(dstu1RS$noxxdr/dstu1all$noxxdr, c(0.25,0.5,0.75))



#################
 

# dsts performed:
# assuming dst required diagnosis but not treatment initiation, need to find columns with a diagnosis (TBstate==2) and compare that column to the following column to determine wehtehr a DST was done
resSEA <- loutcomeboot(individualoutcomefunction = cohortres, simoutput=dst, c=SEAc)
resSAf <- loutcomeboot(individualoutcomefunction = cohortres, simoutput=dst, c=SAfc)
resRMZ <- loutcomeboot(individualoutcomefunction = cohortres, simoutput=dst, c=highrmzc)
reslowMZ <- loutcomeboot(individualoutcomefunction = cohortres, simoutput=dst, c=lowmzc)
save(resSEA, resSAf, resRMZ, reslowMZ, file=paste0("res.",date,".Rdata"))

tbdSEA <- loutcomeboot(individualoutcomefunction = tbdetected, simoutput=dst, c=SEAc)
tbdSAf <- loutcomeboot(individualoutcomefunction = tbdetected, simoutput=dst, c=SAfc)
tbdRMZ <- loutcomeboot(individualoutcomefunction = tbdetected, simoutput=dst, c=highrmzc)
tbdlowMZ <- loutcomeboot(individualoutcomefunction = tbdetected, simoutput=dst, c=lowmzc)
save(tbdSEA, tbdSAf, tbdRMZ, tbdlowMZ, file=paste0("tbd.",date,".Rdata"))


fqdSEA <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=SEAc)
fqdSAf <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=SAfc)
fqdRMZ <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=highrmzc)
fqdlowMZ <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=lowmzc)
fqdRM <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=highrmc)
fqdRZ <- loutcomeboot(individualoutcomefunction = fqrdetected, simoutput=dst, c=highrzc)
save(fqdSEA, fqdSAf, fqdRMZ, fqdlowMZ, fqdMZ,fqdRZ,file=paste0("fqd.",date,".Rdata"))

dstdstSEA <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=SEAc)
save(dstdstSEA, file=paste0("dstdstSEA.",date,".Rdata"))

dstdstSAf <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=SAfc)
save(dstdstSAf, file=paste0("dstdstSAf.",date,".Rdata"))

dstdstRMZ <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=highrmzc)
save(dstdstRMZ, file=paste0("dstdstRMZ.",date,".Rdata"))

dstdstlowMZ <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=lowmzc)
save(dstdstlowMZ, file=paste0("dstdstlowMZ.",date,".Rdata"))

dstdstRZ <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=highrzc)
save(dstdstRZ, file=paste0("dstdstRZ.",date,".Rdata"))

dstdstRM <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=highrmc)
save(dstdstRM, file=paste0("dstdstRM.",date,".Rdata"))


## RIF-R patients, DSTs, successfull treatmnets:

#SAf detection and utilization (for supplemental figure and for text)
lapply( tbdSAf, function(x) quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75))) # TB detected
lapply( dstdstSAf, function(x) quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75))) # RIF DSTs
quantile(apply(dstdstSAf$stepxxdr[,2,],2,sum),  c(0.25,0.5,0.75)  )  # RR Detected  (i.e. FQ DSTs in stepwise case)
lapply(dstdstSAf, function(x) quantile((apply(x[,2,],2,sum)),  c(0.25,0.5,0.75)))    # FQ DSTs
lapply(fqdSAf, function(x) quantile((apply(x[,1,],2,sum)),  c(0.25,0.5,0.75), 3))    # MFX-R detected (uses of backup regimen)

lapply(lsuccess1SAf, function(y) quantile((apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)), # MFX-R treated successfully, RS
                                          c(0.25,0.5,0.75)))
lapply(lsuccess1SAf, function(y) quantile((apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),  # MFX-R treated successfully, RR
                                          c(0.25,0.5,0.75)))


#SEA:
lapply(lsuccess1SEA, function(y) quantile((apply(y[,2,]*(y[,3,]==0)*(y[,4,]==1),2,sum)), # MFX-R treated successfully, RS
                                          c(0.25,0.5,0.75)))
lapply(lsuccess1SEA, function(y) quantile((apply(y[,2,]*(y[,3,]==1)*(y[,4,]==1),2,sum)),  # MFX-R treated successfully, RR
                                          c(0.25,0.5,0.75)))




# initial cycle ratios, step v none, SAf
quantile(apply(dstdstSAf$stepxxdr[,2,],2,sum), c(0.25,0.5,0.75))# tests
quantile(apply(dstdstSAf$stepxxdr[,2,],2,sum)/apply(fqdSAf$stepxxdr[,1,],2,sum), c(0.25,0.5,0.75))# tests per detection
1/quantile(1/(apply(dstdstSAf$stepxxdr[,2,],2,sum)/(dstu1RRSAf$noxxdr - dstu1RRSAf$stepxxdr)), c(0.25,0.5,0.75)); # tests per cure
1/quantile(1/(apply(dstdstSAf$stepxxdr[,4,],2,sum)/(dsttimeRRSAf$noxxdr - dsttimeRRSAf$stepxxdr)), c(0.25,0.5,0.75));  # tests per month mfx-r reduction

# initial cycle ratios, step v none,  SEA
quantile(apply(dstdstSEA$stepxxdr[,2,],2,sum), c(0.25,0.5,0.75)) # tests
quantile(apply(dstdstSEA$stepxxdr[,2,],2,sum)/apply(fqdSEA$stepxxdr[,1,],2,sum), c(0.25,0.5,0.75)) # tests per detection
1/quantile(1/(apply(dstdstSEA$stepxxdr[,2,],2,sum)/(dstu1RR$noxxdr - dstu1RR$stepxxdr)), c(0.25,0.5,0.75)); 
1/quantile(1/(apply(dstdstSEA$stepxxdr[,4,],2,sum)/(dsttimeRR$noxxdr - dsttimeRR$stepxxdr)), c(0.25,0.5,0.75));  # tests per month mfx-r reduction


# initial cycle ratios, full vs step, SAf
quantile(apply(dstdstSAf$fullxxdr[,2,],2,sum) - apply(dstdstSAf$stepxxdr[,2,],2,sum), c(0.25,0.5,0.75)) # tests
quantile((apply(dstdstSAf$fullxxdr[,2,],2,sum) - apply(dstdstSAf$stepxxdr[,2,],2,sum))/(apply(fqdSAf$fullxxdr[,1,],2,sum) - apply(fqdSAf$stepxxdr[,1,],2,sum)), c(0.25,0.5,0.75)) # tests per detection
1/quantile(1/((apply(dstdstSAf$fullxxdr[,2,],2,sum) - (apply(dstdstSAf$stepxxdr[,2,],2,sum)))/(dstu1allSAf$stepxxdr - dstu1allSAf$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstSAf$fullxxdr[,2,],2,sum) - (apply(dstdstSAf$stepxxdr[,2,],2,sum)))/(dsttimeallSAf$stepxxdr - dsttimeallSAf$fullxxdr)), c(0.25,0.5,0.75))   # tests per month mfx-r reduction

# initial cycle ratios, full vs step, SEA
quantile(apply(dstdstSEA$fullxxdr[,2,],2,sum) - apply(dstdstSEA$stepxxdr[,2,],2,sum), c(0.25,0.5,0.75)) # tests
quantile((apply(dstdstSEA$fullxxdr[,2,],2,sum) - apply(dstdstSEA$stepxxdr[,2,],2,sum))/(apply(fqdSEA$fullxxdr[,1,],2,sum) - apply(fqdSEA$stepxxdr[,1,],2,sum)), c(0.25,0.5,0.75)) # tests per detection
quantile((apply(dstdstSEA$fullxxdr[,2,],2,sum) - (apply(dstdstSEA$stepxxdr[,2,],2,sum)))/(dstu1all$stepxxdr - dstu1all$fullxxdr), c(0.25,0.5,0.75))
quantile((apply(dstdstSEA$fullxxdr[,2,],2,sum) - (apply(dstdstSEA$stepxxdr[,2,],2,sum)))/(dsttimeall$stepxxdr - dsttimeall$fullxxdr), c(0.25,0.5,0.75))   # tests per month mfx-r reduction

# full vs none
1/quantile(1/((apply(dstdstSAf$fullxxdr[,2,],2,sum) )/(dstu1allSAf$noxxdr - dstu1allSAf$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstSEA$fullxxdr[,2,],2,sum) )/(dstu1all$noxxdr - dstu1all$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstRMZ$fullxxdr[,2,],2,sum) )/(dstu1all_RMZ$noxxdr - dstu1all_RMZ$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstlowMZ$fullxxdr[,2,],2,sum) )/(dstu1all_lowMZ$noxxdr - dstu1all_lowMZ$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstRZ$fullxxdr[,2,],2,sum) )/(dstu1all_RZ$noxxdr - dstu1all_RZ$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(dstdstRM$fullxxdr[,2,],2,sum) )/(dstu1all_RM$noxxdr - dstu1all_RM$fullxxdr)), c(0.25,0.5,0.75))

# sensiviity analysis
quantile((apply(dstdstRMZ$stepxxdr[,2,],2,sum)/(dstu1RR_RMZ$noxxdr - dstu1RR_RMZ$stepxxdr)), c(0.25,0.5,0.75)); # tests per cure
quantile((apply(dstdstlowMZ$stepxxdr[,2,],2,sum)/(dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)), c(0.25,0.5,0.75)); # tests per cure
mean(apply(dstdstlowMZ$stepxxdr[,2,],2,sum)/(dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)>0)
mean(apply(dstdstlowMZ$stepxxdr[,2,],2,sum)>0)
mean((dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)>0)
median(apply(dstdstlowMZ$stepxxdr[,2,],2,sum))/median(dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)
temp <- (apply(dstdstlowMZ$stepxxdr[,2,],2,sum))/(dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)
temp[temp<0] <- 100000
quantile(temp, c(0.25,0.5,0.75))

quantile((apply(dstdstRMZ$fullxxdr[,2,],2,sum) - (apply(dstdstRMZ$stepxxdr[,2,],2,sum)))/(dstu1all_RMZ$stepxxdr - dstu1all_RMZ$fullxxdr), c(0.25,0.5,0.75))
quantile((apply(dstdstlowMZ$fullxxdr[,2,],2,sum) - (apply(dstdstlowMZ$stepxxdr[,2,],2,sum)))/(dstu1all_lowMZ$stepxxdr - dstu1all_lowMZ$fullxxdr), c(0.25,0.5,0.75))



# for table 2:
load("dstdstSAf.20190711.Rdata")
load("dstdstSEA.20190711.Rdata")
load("dstdstRMZ.20190711.Rdata")
load("dstdstlowMZ.20190711.Rdata")
load("dstdstRZ.20190711.Rdata")
load("dstdstRM.20190711.Rdata")
load("lsuccess.dst.20190711.Rdata")
load("ludst.20190711.Rdata")
load("ludstSAf.20190711.Rdata")
load("ludstRMZ.20190711.Rdata")
load("ludstlowMZ.20190711.Rdata")
load("ludstMZ.20190711.Rdata")
load("ludstRZ.20190711.Rdata")
load(paste0("fqd.",date,".Rdata"))
dstu1RSSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RRSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1allSAf <- lapply(ludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR <- lapply(ludst, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all <- lapply(ludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RMZ <- lapply(ludstRMZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_lowMZ <- lapply(ludstlowMZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RM <- lapply(ludstRM, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all_RZ <- lapply(ludstRZ, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))



#SAf
table <- cbind(
  round(unlist(lapply(dstdstSAf, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdSAf, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1SAf, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                c(0.25,0.5,0.75)), 2))),
  rev(round(1/quantile(1/(apply(dstdstSAf$stepxxdr[,2,],2,sum)/(dstu1RRSAf$noxxdr - dstu1RRSAf$stepxxdr)), c(0.25,0.5,0.75)),0)),
  rev(round(1/quantile(1/((apply(dstdstSAf$fullxxdr[,2,],2,sum) - apply(dstdstSAf$stepxxdr[,2,],2,sum))/(dstu1allSAf$stepxxdr - dstu1allSAf$fullxxdr)), c(0.25,0.5,0.75)),-1)) 
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2SAf <- c(table2[1:6], 
               paste0(table2[7:9]," [",table2[10:12],"]"), 
               paste0(table2[13:15]," [",table2[16:18],"]"),
               table2[c(19, 22)])


# need to repeat for #SEA and others, and then combine:
table <- cbind(
  round(unlist(lapply(dstdstSEA, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdSEA, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1SEA, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1SEA, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1SEA, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1SEA, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                             c(0.25,0.5,0.75)), 2))),
  rev(round(1/quantile(1/(apply(dstdstSEA$stepxxdr[,2,],2,sum)/(dstu1RR$noxxdr - dstu1RR$stepxxdr)), c(0.25,0.5,0.75)),0)),
  rev(round(1/quantile(1/((apply(dstdstSEA$fullxxdr[,2,],2,sum) - apply(dstdstSEA$stepxxdr[,2,],2,sum))/(dstu1all$stepxxdr - dstu1all$fullxxdr)), c(0.25,0.5,0.75)),-1)) 
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2SEA <- c(table2[1:6], 
  paste0(table2[7:9]," [",table2[10:12],"]"), 
  paste0(table2[13:15]," [",table2[16:18],"]"),
  table2[c(19, 22)])


table <- cbind(
  round(unlist(lapply(dstdstRMZ, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdRMZ, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1RMZ, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1RMZ, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1RMZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1RMZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                             c(0.25,0.5,0.75)), 2))),
  rev(round(1/quantile(1/(apply(dstdstRMZ$stepxxdr[,2,],2,sum)/(dstu1RR_RMZ$noxxdr - dstu1RR_RMZ$stepxxdr)), c(0.25,0.5,0.75)),1)),
  rev(round(1/quantile(1/((apply(dstdstRMZ$fullxxdr[,2,],2,sum) - apply(dstdstRMZ$stepxxdr[,2,],2,sum))/(dstu1all_RMZ$stepxxdr - dstu1all_RMZ$fullxxdr)), c(0.25,0.5,0.75)),-1) )
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2RMZ <- c(table2[1:6], 
               paste0(table2[7:9]," [",table2[10:12],"]"), 
               paste0(table2[13:15]," [",table2[16:18],"]"),
               table2[c(19, 22)])

lownum <- apply(dstdstlowMZ$stepxxdr[,2,],2,sum)
lowdenom <- (dstu1RR_lowMZ$noxxdr - dstu1RR_lowMZ$stepxxdr)
lowdenom[lowdenom<0] <- 0.0001
round(quantile(lownum/lowdenom, c(0.25,0.5,0.75)),-1)

table <- cbind(
  round(unlist(lapply(dstdstlowMZ, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdlowMZ, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1lowMZ, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1lowMZ, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1lowMZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1lowMZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                             c(0.25,0.5,0.75)), 2))),
  round(quantile(lownum/lowdenom, c(0.25,0.5,0.75)),-1),
  rev(round(1/quantile(1/((apply(dstdstlowMZ$fullxxdr[,2,],2,sum) - apply(dstdstlowMZ$stepxxdr[,2,],2,sum))/(dstu1all_lowMZ$stepxxdr - dstu1all_lowMZ$fullxxdr)), c(0.25,0.5,0.75)),-1) )
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2lowMZ <- c(table2[1:6], 
               paste0(table2[7:9]," [",table2[10:12],"]"), 
               paste0(table2[13:15]," [",table2[16:18],"]"),
               table2[c(19, 22)])



table <- cbind(
  round(unlist(lapply(dstdstRM, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdRM, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1RM, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1RM, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1RM, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1RM, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                             c(0.25,0.5,0.75)), 2))),
  rev(round(1/quantile(1/(apply(dstdstRM$stepxxdr[,2,],2,sum)/(dstu1RR_RM$noxxdr - dstu1RR_RM$stepxxdr)), c(0.25,0.5,0.75)),1)),
  rev(round(1/quantile(1/((apply(dstdstRM$fullxxdr[,2,],2,sum) - apply(dstdstRM$stepxxdr[,2,],2,sum))/(dstu1all_RM$stepxxdr - dstu1all_RM$fullxxdr)), c(0.25,0.5,0.75)),-1) )
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2RM <- c(table2[1:6], 
               paste0(table2[7:9]," [",table2[10:12],"]"), 
               paste0(table2[13:15]," [",table2[16:18],"]"),
               table2[c(19, 22)])


table <- cbind(
  round(unlist(lapply(dstdstRZ, function(x) quantile((apply(x[,2,],2,sum)), c(0.25,0.5,0.75)))), c(0,0,0,-1,-1,-1,-2,-2,-2,0,0,0)),
  unlist(lapply(fqdRZ, function(x) round(quantile((apply(x[,1,],2,sum)), c(0.25,0.5,0.75)), 0))),
  unlist(lapply(lsuccess1RZ, function(y) round(quantile((apply(y[,2,],2,sum)), c(0.25,0.5,0.75)), -2))),
  100*unlist(lapply(lsuccess1RZ, function(y) round(quantile((apply(y[,2,],2,sum)/apply(y[,1,],2,sum)), c(0.25,0.5,0.75)), 3))) , 
  unlist(lapply(lsuccess1RZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)),c(0.25,0.5,0.75))))) , 
  100*unlist(lapply(lsuccess1RZ, function(y) round(quantile((apply(y[,2,]*(y[,4,]==1),2,sum)/apply(y[,1,]*(y[,4,]==1),2,sum)),
                                                             c(0.25,0.5,0.75)), 2))),
  rev(round(1/quantile(1/(apply(dstdstRZ$stepxxdr[,2,],2,sum)/(dstu1RR_RZ$noxxdr - dstu1RR_RZ$stepxxdr)), c(0.25,0.5,0.75)),1)),
  rev(round(1/quantile(1/((apply(dstdstRZ$fullxxdr[,2,],2,sum) - apply(dstdstRZ$stepxxdr[,2,],2,sum))/(dstu1all_RZ$stepxxdr - dstu1all_RZ$fullxxdr)), c(0.25,0.5,0.75)),-1) )
)
table[c(2,3,5,6,8,9),c(4,6)] <- paste0(table[c(2,3,5,6,8,9),c(4,6)] ,"%")
table2 <- paste0(table[c(2,5,8),]," (", table[c(1,4,7),],"-",table[c(3,6,9),],")")
table2RZ <- c(table2[1:6], 
               paste0(table2[7:9]," [",table2[10:12],"]"), 
               paste0(table2[13:15]," [",table2[16:18],"]"),
               table2[c(19, 22)])


(DSTtable2 <- cbind(table2SAf, table2SEA, table2RMZ, table2lowMZ, table2RM, table2RZ))
write.csv(DSTtable2, file="DSTtable2.csv")



# increase in cure %:
# moxi (full versus step, step versus none, full versus none)
quantile (apply(lsuccess1SAf$fullxxdr[,2,]*(lsuccess1SAf$fullxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$fullxxdr[,1,]*(lsuccess1SAf$fullxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SAf$stepxxdr[,2,]*(lsuccess1SAf$stepxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$stepxxdr[,1,]*(lsuccess1SAf$stepxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
quantile (apply(lsuccess1SAf$stepxxdr[,2,]*(lsuccess1SAf$stepxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$stepxxdr[,1,]*(lsuccess1SAf$stepxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SAf$noxxdr[,2,]*(lsuccess1SAf$noxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$noxxdr[,1,]*(lsuccess1SAf$noxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
quantile (apply(lsuccess1SAf$fullxxdr[,2,]*(lsuccess1SAf$fullxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$fullxxdr[,1,]*(lsuccess1SAf$fullxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SAf$noxxdr[,2,]*(lsuccess1SAf$noxxdr[,4,]==1),2,sum)/apply(lsuccess1SAf$noxxdr[,1,]*(lsuccess1SAf$noxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# rr (step versus none)
quantile (apply(lsuccess1SAf$stepxxdr[,2,]*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$stepxxdr[,1,]*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum) - 
            apply(lsuccess1SAf$noxxdr[,2,]*(lsuccess1SAf$noxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$noxxdr[,1,]*(lsuccess1SAf$noxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# rr, mfx-r (step versus none)
quantile (apply(lsuccess1SAf$stepxxdr[,2,]*(lsuccess1SAf$stepxxdr[,4,]==1)*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$stepxxdr[,1,]*(lsuccess1SAf$stepxxdr[,4,]==1)*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum) - 
            apply(lsuccess1SAf$noxxdr[,2,]*(lsuccess1SAf$noxxdr[,4,]==1)*(lsuccess1SAf$noxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$noxxdr[,1,]*(lsuccess1SAf$noxxdr[,4,]==1)*(lsuccess1SAf$noxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# and all (full versus none)
quantile (apply(lsuccess1SAf$fullxxdr[,2,],2,sum)/apply(lsuccess1SAf$fullxxdr[,1,],2,sum) - 
            apply(lsuccess1SAf$noxxdr[,2,],2,sum)/apply(lsuccess1SAf$noxxdr[,1,],2,sum) , c(0.25,0.5,0.75), 3) 

# proportion cured:
quantile(apply(lsuccess1SAf$stepxxdr[,2,]*(lsuccess1SAf$stepxxdr[,4,]==1)*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$stepxxdr[,1,]*(lsuccess1SAf$stepxxdr[,4,]==1)*(lsuccess1SAf$stepxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3)
quantile(apply(lsuccess1SAf$noxxdr[,2,]*(lsuccess1SAf$noxxdr[,4,]==1)*(lsuccess1SAf$noxxdr[,3,]==1),2,sum)/apply(lsuccess1SAf$noxxdr[,1,]*(lsuccess1SAf$noxxdr[,4,]==1)*(lsuccess1SAf$noxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3) 


# moxi (full versus step, step versus none, full versus none)
quantile (apply(lsuccess1SEA$fullxxdr[,2,]*(lsuccess1SEA$fullxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$fullxxdr[,1,]*(lsuccess1SEA$fullxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SEA$stepxxdr[,2,]*(lsuccess1SEA$stepxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$stepxxdr[,1,]*(lsuccess1SEA$stepxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
quantile (apply(lsuccess1SEA$stepxxdr[,2,]*(lsuccess1SEA$stepxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$stepxxdr[,1,]*(lsuccess1SEA$stepxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SEA$noxxdr[,2,]*(lsuccess1SEA$noxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$noxxdr[,1,]*(lsuccess1SEA$noxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
quantile (apply(lsuccess1SEA$fullxxdr[,2,]*(lsuccess1SEA$fullxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$fullxxdr[,1,]*(lsuccess1SEA$fullxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1SEA$noxxdr[,2,]*(lsuccess1SEA$noxxdr[,4,]==1),2,sum)/apply(lsuccess1SEA$noxxdr[,1,]*(lsuccess1SEA$noxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# rr (step versus none)
quantile (apply(lsuccess1SEA$stepxxdr[,2,]*(lsuccess1SEA$stepxxdr[,3,]==1),2,sum)/apply(lsuccess1SEA$stepxxdr[,1,]*(lsuccess1SEA$stepxxdr[,3,]==1),2,sum) - 
            apply(lsuccess1SEA$noxxdr[,2,]*(lsuccess1SEA$noxxdr[,3,]==1),2,sum)/apply(lsuccess1SEA$noxxdr[,1,]*(lsuccess1SEA$noxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# rr, mfx-r (step versus none)
quantile (apply(lsuccess1SEA$stepxxdr[,2,]*(lsuccess1SEA$stepxxdr[,4,]==1)*(lsuccess1SEA$stepxxdr[,3,]==1),2,sum)/apply(lsuccess1SEA$stepxxdr[,1,]*(lsuccess1SEA$stepxxdr[,4,]==1)*(lsuccess1SEA$stepxxdr[,3,]==1),2,sum) - 
            apply(lsuccess1SEA$noxxdr[,2,]*(lsuccess1SEA$noxxdr[,4,]==1)*(lsuccess1SEA$noxxdr[,3,]==1),2,sum)/apply(lsuccess1SEA$noxxdr[,1,]*(lsuccess1SEA$noxxdr[,4,]==1)*(lsuccess1SEA$noxxdr[,3,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
# and all (full versus none)
quantile (apply(lsuccess1SEA$fullxxdr[,2,],2,sum)/apply(lsuccess1SEA$fullxxdr[,1,],2,sum) - 
            apply(lsuccess1SEA$noxxdr[,2,],2,sum)/apply(lsuccess1SEA$noxxdr[,1,],2,sum) , c(0.25,0.5,0.75), 3) 

quantile (apply(lsuccess1RMZ$fullxxdr[,2,]*(lsuccess1RMZ$noxxdr[,4,]==1),2,sum)/apply(lsuccess1RMZ$fullxxdr[,1,]*(lsuccess1RMZ$noxxdr[,4,]==1),2,sum) - 
            apply(lsuccess1RMZ$noxxdr[,2,]*(lsuccess1RMZ$noxxdr[,4,]==1),2,sum)/apply(lsuccess1RMZ$noxxdr[,1,]*(lsuccess1RMZ$noxxdr[,4,]==1),2,sum) , c(0.25,0.5,0.75), 3) 
quantile (apply(lsuccess1RMZ$fullxxdr[,2,],2,sum)/apply(lsuccess1RMZ$fullxxdr[,1,],2,sum) - 
            apply(lsuccess1RMZ$noxxdr[,2,],2,sum)/apply(lsuccess1RMZ$noxxdr[,1,],2,sum) , c(0.25,0.5,0.75), 3) 


# scenario analysis, 6 month rr duration
load("dst6m.20190711.Rdata")
ludst6m <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst6m,c = SEAc)
ludst6mSAf <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst6m,c = SAfc)
dstu1RS6mSEA <- lapply(ludst6m, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR6mSEA <- lapply(ludst6m, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all6mSEA <- lapply(ludst6m, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstu1RS6mSAf <- lapply(ludst6mSAf, function(y) apply(y[,1,]*(y[,3,]==0)*(y[,4,]==1),2,sum))
dstu1RR6mSAf <- lapply(ludst6mSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1all6mSAf <- lapply(ludst6mSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
dstdst6mSEA <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst6m, c=SEAc)
dstdst6mSAf <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst6m, c=SAfc)

1/quantile(1/(apply(dstdst6mSAf$stepxxdr[,2,],2,sum)/(dstu1RR6mSAf$noxxdr - dstu1RR6mSAf$stepxxdr)), c(0.25,0.5,0.75)); 
1/quantile(1/((apply(dstdst6mSAf$fullxxdr[,2,],2,sum) - (apply(dstdst6mSAf$stepxxdr[,2,],2,sum)))/(dstu1all6mSAf$stepxxdr - dstu1all6mSAf$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/(apply(dstdst6mSAf$fullxxdr[,2,],2,sum)/(dstu1all6mSAf$noxxdr - dstu1all6mSAf$fullxxdr)), c(0.25,0.5,0.75)); 

1/quantile(1/(apply(dstdst6mSEA$stepxxdr[,2,],2,sum)/(dstu1RR6mSEA$noxxdr - dstu1RR6mSEA$stepxxdr)), c(0.25,0.5,0.75)); 
1/quantile(1/((apply(dstdst6mSEA$fullxxdr[,2,],2,sum) - (apply(dstdst6mSEA$stepxxdr[,2,],2,sum)))/(dstu1all6mSEA$stepxxdr - dstu1all6mSEA$fullxxdr)), c(0.25,0.5,0.75))




# ###########################333 
# # sensitivity analysis, DST among RR -- but detailed results in supplement
# 
# #1. Impact of moxi res on regimen robustness
load(paste0("robustdst.",date,".Rdata"))
dstdstrobust <- loutcomeboot(individualoutcomefunction = dsts, simoutput=robustdst, c=SAfc)
ludstrobust <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= robustdst, c = SAfc)
dstu1RRrobust <- lapply(ludstrobust, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
dstu1robust <- lapply(ludstrobust, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
1/quantile(1/(apply(dstdstrobust$stepxxdr[,2,],2,sum)/(dstu1RRrobust$noxxdr - dstu1RRrobust$stepxxdr)), c(0.25,0.5,0.75));
1/quantile(1/((apply(dstdstrobust$fullxxdr[,2,],2,sum) - (apply(dstdstrobust$stepxxdr[,2,],2,sum)))/(dstu1robust$stepxxdr - dstu1robust$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/(apply(dstdstrobust$fullxxdr[,2,],2,sum)/(dstu1robust$noxxdr - dstu1robust$fullxxdr)), c(0.25,0.5,0.75)); 
# 
# 
# #2. Increase PZA or M res but not both
# dstdstm <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=highrmc)
# ludstm <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highrmc)
# dstu1RRm <- lapply(ludstm, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
# 1/quantile(1/(apply(dstdstm$stepxxdr[,2,],2,sum)/(dstu1RRm$noxxdr - dstu1RRm$stepxxdr)), c(0.25,0.5,0.75)); 
# 
# dstdstz <- loutcomeboot(individualoutcomefunction = dsts, simoutput=dst, c=highrzc)
# ludstz <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = highrzc)
# dstu1RRz <- lapply(ludstz, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
# 1/quantile(1/(apply(dstdstz$stepxxdr[,2,],2,sum)/(dstu1RRz$noxxdr - dstu1RRz$stepxxdr)), c(0.25,0.5,0.75)); 
# 
# #[3? Strength of correlation of Moxi and PZA resistance]
# 




###################################### function for infectious time with drug resistance: 
###############

# Second manuscript resistance outcomes and final adr figure:

lrestreatedSAf <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = dst, c=SAfc)
lrestreatedSEA <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = dst, c=SEAc)
lrestreatedhighmz <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = dst, c=highmzc)
lrestreatedhighrmz <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = dst, c=highrmzc)

save(lrestreatedSAf, lrestreatedSEA, lrestreatedhighmz, lrestreatedhighrmz, file=paste0("lrestreated_dst.",date,".Rdata"))

load(paste0("highresdst.",date,".Rdata"))
lhighrestreatedSAf <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = highresdst, c=SAfc)
lhighrestreatedSEA <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = highresdst, c=SEAc)
lhighrestreatedhighmz <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = highresdst, c=highmzc)
lhighrestreatedhighrmz <- loutcomeboot(individualoutcomefunction = infectioustime, includeinitial=FALSE, simoutput = highresdst, c=highrmzc)

save(lhighrestreatedSAf, lhighrestreatedSEA, lhighrestreatedhighmz, lhighrestreatedhighrmz, file=paste0("lhighrestreated_dst.",date,".Rdata"))
rm(highresdst)

# quick check: plots of bdq res
boxplot(lapply(lrestreatedSAf, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lrestreatedSEA, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lrestreatedhighmz, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lrestreatedhighrmz, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lhighrestreatedSAf, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lhighrestreatedSEA, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lhighrestreatedhighmz, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))
boxplot(lapply(lhighrestreatedhighrmz, function(x) apply(x[,6,],2,mean)/apply(x[,7,],2,mean)))

# for 2nd manu figure, want treated time (abs, not proportion),  all 3 scenarios, 4 settings, and 3 drugs (no PA, and maybe also eliminate rif)

a1 <- data.frame(lapply(lrestreatedSAf, function(y) apply(y[,1,],2,sum)))
a2 <- data.frame(lapply(lrestreatedSAf, function(y) apply(y[,2,],2,sum)))
a3 <- data.frame(lapply(lrestreatedSAf, function(y) apply(y[,4,],2,sum)))
a1$drug <- "Rifampin"
a2$drug <- "Moxifloxacin"
a3$drug <- "Bedaquiline"
a <- bind_rows(a1,a2,a3)
am <- melt(a, id.vars=c("drug"))
am$params <- "Low ADR"
am$setting <- "South Africa"

c1 <- data.frame(lapply(lrestreatedSEA, function(y) apply(y[,1,],2,sum)))
c2 <- data.frame(lapply(lrestreatedSEA, function(y) apply(y[,2,],2,sum)))
c3 <- data.frame(lapply(lrestreatedSEA, function(y) apply(y[,4,],2,sum)))
c1$drug <- "Rifampin"
c2$drug <- "Moxifloxacin"
c3$drug <- "Bedaquiline"
c <- bind_rows(c1,c2,c3)
cm <- melt(c, id.vars=c("drug"))
cm$params <- "Low ADR"
cm$setting <- "Southeast Asia"

# e1 <- data.frame(lapply(lrestreatedhighmz, function(y) apply(y[,1,],2,sum)))
# e2 <- data.frame(lapply(lrestreatedhighmz, function(y) apply(y[,2,],2,sum)))
# e3 <- data.frame(lapply(lrestreatedhighmz, function(y) apply(y[,4,],2,sum)))
# e1$drug <- "Rifampin"
# e2$drug <- "Moxifloxacin"
# e3$drug <- "Bedaquiline"
# e <- bind_rows(e1,e2,e3)
# em <- melt(e, id.vars=c("drug"))
# em$params <- "Low ADR"

g1 <- data.frame(lapply(lrestreatedhighrmz, function(y) apply(y[,1,],2,sum)))
g2 <- data.frame(lapply(lrestreatedhighrmz, function(y) apply(y[,2,],2,sum)))
g3 <- data.frame(lapply(lrestreatedhighrmz, function(y) apply(y[,4,],2,sum)))
g1$drug <- "Rifampin"
g2$drug <- "Moxifloxacin"
g3$drug <- "Bedaquiline"
g <- bind_rows(g1,g2,g3)
gm <- melt(g, id.vars=c("drug"))
gm$params <- "Low ADR"
gm$setting <- "High-resistance setting (~Belarus)"

# b1 <- data.frame(lapply(lhighrestreatedSAf, function(y) apply(y[,1,],2,sum)))
# b2 <- data.frame(lapply(lhighrestreatedSAf, function(y) apply(y[,2,],2,sum)))
# b3 <- data.frame(lapply(lhighrestreatedSAf, function(y) apply(y[,4,],2,sum)))
# b1$drug <- "Rifampin"
# b2$drug <- "Moxifloxacin"
# b3$drug <- "Bedaquiline"
# b <- bind_rows(b1,b2,b3)
# bm <- melt(b, id.vars=c("drug"))
# bm$params <- "High ADR"
# require(dplyr)
# require(data.table)
# mm <- bind_rows(am,bm)
# colnames(mm) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
# mm$Parameters <- factor(mm$Parameters, levels=c("Low ADR", "High ADR"))
# levels(mm$Parameters) <- c("Low", "High")
# mm$setting <- "South Africa"
# 
# d1 <- data.frame(lapply(lhighrestreatedSEA, function(y) apply(y[,1,],2,sum)))
# d2 <- data.frame(lapply(lhighrestreatedSEA, function(y) apply(y[,2,],2,sum)))
# d3 <- data.frame(lapply(lhighrestreatedSEA, function(y) apply(y[,4,],2,sum)))
# d1$drug <- "Rifampin"
# d2$drug <- "Moxifloxacin"
# d3$drug <- "Bedaquiline"
# d <- bind_rows(d1,d2,d3)
# dm <- melt(d, id.vars=c("drug"))
# dm$params <- "High ADR"
# require(dplyr)
# require(data.table)
# mm2 <- bind_rows(cm,dm)
# colnames(mm2) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
# mm2$Parameters <- factor(mm2$Parameters, levels=c("Low ADR", "High ADR"))
# levels(mm2$Parameters) <- c("Low", "High")
# mm2$setting <- "Southeast Asia"
# 
# f1 <- data.frame(lapply(lhighrestreatedhighmz, function(y) apply(y[,1,],2,sum)))
# f2 <- data.frame(lapply(lhighrestreatedhighmz, function(y) apply(y[,2,],2,sum)))
# f3 <- data.frame(lapply(lhighrestreatedhighmz, function(y) apply(y[,4,],2,sum)))
# f1$drug <- "Rifampin"
# f2$drug <- "Moxifloxacin"
# f3$drug <- "Bedaquiline"
# f <- bind_rows(f1,f2,f3)
# fm <- melt(f, id.vars=c("drug"))
# fm$params <- "High ADR"
# require(dplyr)
# require(data.table)
# mm3 <- bind_rows(em,fm)
# colnames(mm3) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
# mm3$Parameters <- factor(mm3$Parameters, levels=c("Low ADR", "High ADR"))
# levels(mm3$Parameters) <- c("Low", "High")
# mm3$setting <- "Higher MOXI-R and PZA-R Prevalence"
# 
# h1 <- data.frame(lapply(lhighrestreatedhighrmz, function(y) apply(y[,1,],2,sum)))
# h2 <- data.frame(lapply(lhighrestreatedhighrmz, function(y) apply(y[,2,],2,sum)))
# h3 <- data.frame(lapply(lhighrestreatedhighrmz, function(y) apply(y[,4,],2,sum)))
# h1$drug <- "Rifampin"
# h2$drug <- "Moxifloxacin"
# h3$drug <- "Bedaquiline"
# h <- bind_rows(h1,h2,h3)
# hm <- melt(h, id.vars=c("drug"))
# hm$params <- "High ADR"
# require(dplyr)
# require(data.table)
# mm4 <- bind_rows(gm,hm)
# colnames(mm4) <- c("Drug", "Scenario", "Infectious_Time", "Parameters")
# mm4$Parameters <- factor(mm4$Parameters, levels=c("Low ADR", "High ADR"))
# levels(mm4$Parameters) <- c("Low", "High")
# mm4$setting <- "Higher RIF-R, MFX-R and PZA-R Prevalence"
# 

# mm5 <- rbind(mm, mm2, mm3, mm4)
mm5 <- rbind(am, cm, gm)
mm5$setting <- as.factor(mm5$setting); #levels(mm5$setting) <- c("Higher MFX-R and\nPZA-R Prevalence", "Higher MFX-R, PZA-R,\nand RR prevalence", "South Africa", "Southeast Asia")
mm5$setting <- factor(mm5$setting, levels = levels(mm5$setting)[c(2,3,1)])

mm5$Drug <- factor(mm5$drug)
mm5$Parameters <- factor(mm5$params)
mm5 <- subset(mm5, variable != "baseline")
mm5$Scenario <- factor(mm5$variable, levels=rev(levels(mm5$variable)))
mm5$Infectious_Time <- mm5$value
levels(mm5$Scenario) <- rev(c("None",
                              "Targeted",
                              "Universal", "x"))

save(am, bm, cm, dm, em, fm, gm, hm, mm, mm2, mm3, mm4, mm5, file="dstFig3data.Rdata")
save(am, cm, gm, mm5, file="dstFig3data.Rdata")

require(ggplot2)
require(gridExtra)
require(wesanderson)
require(RColorBrewer)
pdf("dstFig3.pdf", width=10, height=8)
colors <- brewer.pal(5, "Blues")[3:5]
par(mar=rep(2,4))
theme = theme_set(theme_bw())
glow <-
  ggplot(data = transform(subset(mm5, Parameters=="Low ADR" & setting %in% c("South Africa", "Southeast Asia") & 
                                   Drug %in% c("Moxifloxacin","Bedaquiline")), Drug=factor(Drug, levels=c("Moxifloxacin","Bedaquiline"))),
         aes(x=Scenario, y=Infectious_Time/1000)) + 
  geom_boxplot(aes(fill=Scenario, colour=Scenario), outlier.colour = NULL) +
  facet_grid(Drug ~ setting) +
  scale_fill_manual(values=colors) + scale_colour_manual(values=colors)  +
  # scale_y_continuous(trans = 'log2', limits = c(1,90)) +
  scale_y_continuous(limits = c(0,10)) +
  coord_flip() + 
  theme(legend.position = "none") + 
  # stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max, aes(fill=Scenario), size = 0.2, position = position_dodge(width=0.8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Use of MOXI DST") + 
  ylab("Post-treatment person-months with active TB resistant to the specified drug, per 100 TB cases (log scale)") +
  ggtitle("A. Assuming lower risks of acquiring resistance") +
  theme(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))+ 
  theme(plot.title = element_text(hjust=0)) +
  theme(plot.margin=unit(c(0.5,0.5,0.75,0.5), "cm"))
ghigh <-
  ggplot(data = transform(subset(mm5, Parameters=="High" & 
                                   Drug %in% c("Moxifloxacin","Bedaquiline")), Drug=factor(Drug, levels=c("Moxifloxacin","Bedaquiline"))),
         aes(x=Scenario, y=Infectious_Time/1000)) + 
  geom_boxplot(aes(fill=Scenario, colour=Scenario), outlier.colour = NULL) +
  facet_grid(Drug ~ setting) +
  scale_fill_manual(values=colors) + scale_colour_manual(values=colors)  +
  scale_y_continuous(trans = 'log2', limits=c(1,90)) +
  coord_flip() + 
  theme(legend.position = "none") + 
  # stat_summary(fun.y = median, fun.ymin = min, fun.ymax = max, aes(fill=Scenario), size = 0.2, position = position_dodge(width=0.8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Use of MOXI DST") + 
  ylab("Post-treatment person-months with active TB resistant to the specified drug, per 100 TB cases (log scale)") +
  ggtitle("B. Assuming higher risks of acquiring resistance") +
  theme(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11)) + 
  theme(plot.title = element_text(hjust=0)) + 
  theme(plot.margin=unit(c(.75,0.5,0.5,0.5), "cm"))
grid.arrange(
  glow, ghigh,
  ncol=1
)
dev.off()


# FOR DST ADR text:
##FQ DST reduced potential MFX-R transmission:
nonetime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="High", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="High", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Moxifloxacin", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

#bdq
nonetime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Bedaquiline", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Bedaquiline", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Bedaquiline", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Bedaquiline", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting==levels(mm5$setting)[3], Drug=="Bedaquiline", Parameters=="Low", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting==levels(mm5$setting)[3], Drug=="Bedaquiline", Parameters=="Low", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Bedaquiline", Parameters=="High", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="South Africa", Drug=="Bedaquiline", Parameters=="High", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Bedaquiline", Parameters=="High", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting=="Southeast Asia", Drug=="Bedaquiline", Parameters=="High", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)

nonetime <- unlist((mm5 %>% filter(setting==levels(mm5$setting)[3], Drug=="Bedaquiline", Parameters=="High", Scenario=="None") %>%  select(Infectious_Time) ))/1000
univtime <- unlist((mm5 %>% filter(setting==levels(mm5$setting)[3], Drug=="Bedaquiline", Parameters=="High", Scenario=="Universal") %>%  select(Infectious_Time) ))/1000
mean(nonetime-univtime); sd(nonetime-univtime)
mean((univtime-nonetime)/nonetime); sd((univtime-nonetime)/nonetime)




# sensitivity analysis cohort cmposition
sum(subset(highmzc, RIF==1)$Freq)/sum(highmzc$Freq)
sum(subset(highmzc, MOXI==1)$Freq)/sum(highmzc$Freq)
sum(subset(highmzc, RIF==0 & MOXI==1)$Freq)/sum(subset(highmzc, RIF==0)$Freq)
sum(subset(SEAc, RIF==1 & MOXI==1)$Freq)/sum(subset(SEAc, RIF==1)$Freq)
sum(subset(SAfc, RIF==1 & MOXI==1)$Freq)/sum(subset(SAfc, RIF==1)$Freq)
sum(subset(highmzc, RIF==1 & MOXI==1)$Freq)/sum(subset(highmzc, RIF==1)$Freq)
sum(subset(highmzc, RIF==1 & MOXI==1)$Freq)/sum(subset(highmzc, RIF==1)$Freq)/
  (sum(subset(highmzc, RIF==0 & MOXI==1)$Freq)/sum(subset(highmzc, RIF==0)$Freq))

sum(subset(highrmzc, RIF==1)$Freq)/sum(highrmzc$Freq)
sum(subset(highrmzc, MOXI==1)$Freq)/sum(highrmzc$Freq)
sum(subset(highrmzc, RIF==1 & MOXI==1)$Freq)/sum(subset(highrmzc, RIF==1)$Freq)
sum(subset(highrmzc, RIF==1 & PZA==1)$Freq)/sum(subset(highrmzc, RIF==1)$Freq)

# for contribution to poor outcomes, go down to ludst
# sensitivity analysis: outcomes with higher m/z res, still strong correlation with RR:
tSAfdsthighmz <- tableSreg(impact = combined, c = highmzc, save = F)
save(tSAfdsthighmz, file=paste0("tSAfdsthighmz.",date,".Rdata"))
tSAfdsthighmz[1:16,] <- round(tSAfdsthighmz[1:16,]*100,1)
tSAfdsthighmz[17:24,] <- round(tSAfdsthighmz[17:24,],2)
tSAf2dsthighmz <- cbind(paste0(tSAfdsthighmz[,1], " \u00b1 ",tSAfdsthighmz[,2]), paste0(tSAfdsthighmz[,3], " \u00b1 ",tSAfdsthighmz[,4]), paste0(tSAfdsthighmz[,5], " \u00b1 ",tSAfdsthighmz[,6]), paste0(tSAfdsthighmz[,7], " \u00b1 ",tSAfdsthighmz[,8]), paste0(tSAfdsthighmz[,9], " \u00b1 ",tSAfdsthighmz[,10]), paste0(tSAfdsthighmz[,11], " \u00b1 ",tSAfdsthighmz[,12]))
tSAf2dsthighmz[1:16,] <- paste0(tSAf2dsthighmz[1:16,],"%")
Encoding(tSAf2dsthighmz)<-"UTF-8"
tSAf2dst <- tSAf2dsthighmz[,c(4,5,6,1,2)]
colnames(tSAf2dsthighmz) <- c("RIF-R-based BPaMZ duration for all",
                              " MOXI DST for RIF-R TB", "MOXI DST for all TB", "4 months BPaMZ for all", "6 months BPaMZ for all")
write.table(tSAf2dsthighmz, file=paste0("tSAf2dsthighmz.",date,".csv"), sep = ",")
save(tSAf2dsthighmz, file=paste0("tSAf2dsthighmz.",date,".Rdata"))
tSAf2dsthighmz

tSAfdsthighrmz <- tableSreg(impact = combined, c = highrmzc, save = F)
save(tSAfdsthighrmz, file=paste0("tSAfdsthighrmz.",date,".Rdata"))
tSAfdsthighrmz[1:16,] <- round(tSAfdsthighrmz[1:16,]*100,1)
tSAfdsthighrmz[17:24,] <- round(tSAfdsthighrmz[17:24,],2)
tSAf2dsthighrmz <- cbind(paste0(tSAfdsthighrmz[,1], " \u00b1 ",tSAfdsthighrmz[,2]), paste0(tSAfdsthighrmz[,3], " \u00b1 ",tSAfdsthighrmz[,4]), paste0(tSAfdsthighrmz[,5], " \u00b1 ",tSAfdsthighrmz[,6]), paste0(tSAfdsthighrmz[,7], " \u00b1 ",tSAfdsthighrmz[,8]), paste0(tSAfdsthighrmz[,9], " \u00b1 ",tSAfdsthighrmz[,10]), paste0(tSAfdsthighrmz[,11], " \u00b1 ",tSAfdsthighrmz[,12]))
tSAf2dsthighrmz[1:16,] <- paste0(tSAf2dsthighrmz[1:16,],"%")
Encoding(tSAf2dsthighrmz)<-"UTF-8"
tSAf2dst <- tSAf2dsthighrmz[,c(4,5,6,1,2)]
colnames(tSAf2dsthighrmz) <- c("RIF-R-based BPaMZ duration for all",
                               " MOXI DST for RIF-R TB", "MOXI DST for all TB", "4 months BPaMZ for all", "6 months BPaMZ for all")
write.table(tSAf2dsthighrmz, file=paste0("tSAf2dsthighrmz.",date,".csv"), sep = ",")
save(tSAf2dsthighrmz, file=paste0("tSAf2dsthighrmz.",date,".Rdata"))
tSAf2dsthighrmz


# hypermutator
# combine the pan-S subset of safc/seac from dst, with the any-drug-r subset of safc/seac from hypermutator.

# SEA first:

noresludst <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = SEAc, 
             include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0),
             desiredsize = round(1e5*sum(SEAc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0])))
# number of poor outcomes in moxi-R:
noresdstu1RR <- lapply(noresludst, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
noresdstu1all <- lapply(noresludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
# number of DSTs:
noresdstdstSEA <- loutcomeboot(individualoutcomefunction = dsts, simoutput = dst, c = SEAc,
                               include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0),
                               desiredsize = round(1e5*sum(SEAc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0])))

load(paste0("hypermutator.",date,".Rdata"))
hypermutludst <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= hypermutator, c = SEAc,
                              include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0),
                              desiredsize = round(1e5*sum(SEAc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0])))
hypermutdstu1RR <- lapply(hypermutludst, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
hypermutdstu1all <- lapply(hypermutludst, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
hypermutdstdst <- loutcomeboot(individualoutcomefunction = dsts, simoutput=hypermutator, c = SEAc,
                               include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0),
                               desiredsize = round(1e5*sum(SEAc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0])))

1/quantile(1/((apply(noresdstdstSEA$stepxxdr[,2,],2,sum) + apply(hypermutdstdst$stepxxdr[,2,],2,sum))/
             (noresdstu1RR$noxxdr + hypermutdstu1RR$noxxdr - noresdstu1RR$stepxxdr - hypermutdstu1RR$stepxxdr)), c(0.25,0.5,0.75));
1/quantile(1/((apply(noresdstdstSEA$fullxxdr[,2,],2,sum) + apply(hypermutdstdst$fullxxdr[,2,],2,sum))/
                (noresdstu1all$stepxxdr + hypermutdstu1all$stepxxdr - noresdstu1all$fullxxdr - hypermutdstu1all$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(noresdstdstSEA$fullxxdr[,2,],2,sum) + apply(hypermutdstdst$fullxxdr[,2,],2,sum))/
                (noresdstu1all$noxxdr + hypermutdstu1all$noxxdr - noresdstu1all$fullxxdr - hypermutdstu1all$fullxxdr)), c(0.25,0.5,0.75))


# Then SAf:

noresludstSAf <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= dst, c = SAfc, 
                           include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0),
                           desiredsize = round(1e5*sum(SAfc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0])))
# number of poor outcomes in moxi-R:
noresdstu1RRSAf <- lapply(noresludstSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
noresdstu1allSAf <- lapply(noresludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
# number of DSTs:
noresdstdstSAf <- loutcomeboot(individualoutcomefunction = dsts, simoutput = dst, c = SEAc,
                               include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0),
                               desiredsize = round(1e5*sum(SAfc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA==0])))

hypermutludstSAf <- loutcomeboot(individualoutcomefunction = unsuccessfulrx, simoutput= hypermutator, c = SEAc,
                              include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0),
                              desiredsize = round(1e5*sum(SAfc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0])))
hypermutdstu1RRSAf <- lapply(hypermutludstSAf, function(y) apply(y[,1,]*(y[,3,]==1)*(y[,4,]==1),2,sum))
hypermutdstu1allSAf <- lapply(hypermutludstSAf, function(y) apply(y[,1,]*(y[,4,]==1),2,sum))
hypermutdstdstSAf <- loutcomeboot(individualoutcomefunction = dsts, simoutput=hypermutator, c = SEAc,
                               include=(SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0),
                               desiredsize = round(1e5*sum(SAfc$Freq[SEAc$RIF+SEAc$INH+SEAc$PZA+SEAc$MOXI+SEAc$BDQ+SEAc$PA>0])))

1/quantile(1/((apply(noresdstdstSAf$stepxxdr[,2,],2,sum) + apply(hypermutdstdstSAf$stepxxdr[,2,],2,sum))/
                (noresdstu1RRSAf$noxxdr + hypermutdstu1RRSAf$noxxdr - noresdstu1RRSAf$stepxxdr - hypermutdstu1RRSAf$stepxxdr)), c(0.25,0.5,0.75));
1/quantile(1/((apply(noresdstdstSAf$fullxxdr[,2,],2,sum) + apply(hypermutdstdstSAf$fullxxdr[,2,],2,sum))/
                (noresdstu1allSAf$stepxxdr + hypermutdstu1allSAf$stepxxdr - noresdstu1allSAf$fullxxdr - hypermutdstu1allSAf$fullxxdr)), c(0.25,0.5,0.75))
1/quantile(1/((apply(noresdstdstSAf$fullxxdr[,2,],2,sum) + apply(hypermutdstdstSAf$fullxxdr[,2,],2,sum))/
                (noresdstu1allSAf$noxxdr + hypermutdstu1allSAf$noxxdr - noresdstu1allSAf$fullxxdr - hypermutdstu1allSAf$fullxxdr)), c(0.25,0.5,0.75))
