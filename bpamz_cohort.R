
# setwd("C:/Users/ekendal2/OneDrive - Johns Hopkins University/Research/universal regimen")

require(arm)
require(dplyr)

allparams <- read.csv("allparams.csv", header=T, stringsAsFactors = F)
if(setting=="SEA") params <- as.numeric(allparams[,2])
if(setting=="SAf") params <- as.numeric(allparams[,5])
names(params) <- allparams[,1]
params["Tbdxtime_recurrence"] <- params["Tbdxtime_recurrenceratio"]*params["Tbdxtime"]

params["Regimendelay"] <- 0
params["Regimenloss"] <- 0

increaseodds <- function(recurrenceprobs, OR)
{
  initialodds <- (1-recurrenceprobs)/recurrenceprobs
  newodds <- OR*initialodds
  newrecurrenceprobs <- 1/(1+newodds)
  return(newrecurrenceprobs)
}

drugs <- c("RIF", "INH", "PZA", "MOXI", "partialmoxi","BDQ", "PA")
regimens <- c("hrze", "mdr", "bpamz4", "bpamz6", "bpal", "backup", "backup_mdr")
# activeregimens <- c("HR(ZE)","R(ZE)","(ZE)","MDR, FQ-S","MDR, FQ-R","BPaMZ","BPaM","BPaZ","BMZ","PaMZ",
#                     "BPa","BZ","BM","PaZ","PaM","MZ","B","Pa","M","Z","BPaL","BL","PaL","L","none")
activeregimens <- c("HR(ZE)","R(ZE)","H(ZE)","(ZE)","MDR, FQ-S","MDR, FQ-low","MDR, FQ-R","BPaMZ","BPaM","BPaZ","BPamZ", "BMZ","PaMZ",
                    "BPa","BPam", "BZ","BmZ", "BM","PaZ","PamZ", "PaM","MZ","B","Bm", "Pa","Pam","M","Z","mZ","BPaL","BL","PaL","L","m","none")
patientvars <- c("RxHist", "NovelHist","HIV", "SmearStatus", drugs)

regimendurations <- c(6, 18, 4, 6, 6, 6, 18)
monthsmodeled <- c(1:6, 9, 12, 18)
# monthsfollowed <- 36

eventtypes <- list(TBonset=1, TBdiagnosis=2, treatmentstart=3, treatmentend=4, death=5, pretreatmentltfu=6, carryforward=7)
# diagnosis will incorporate DST and regimen selection), treatment start will incorporate ADR, 
# For each event, track RxHist, NovelHist, current susceptibilities, TB state.
statetypes <- list(undiagnosed=1, diagnosed=2, treating=3, treating_adr=4, failed=5, pendingrelapse=6, cured=7, deceased=8, relapsed=9)
regimentypes <- as.list(1:length(regimens)); names(regimentypes) <- regimens; regimentypes$none <- 0


# Set up cohort of patients
# to look only at a particular patient type, use diferent params, 
# or to keep other ratios correct (e.g. new:retreatment whensubsetting to RR), make a huge cohort and then filter.
make.cohort <- function(params, patientvars, N=10000)
{
  cohort <- array(NA, dim=c(N,length(patientvars))); colnames(cohort)=patientvars
  cohort <- data.frame(cohort)
  
  cohort$RxHist <- sample(c(0,1), size = nrow(cohort), replace = T, 
                          prob=c(1-params['Retreatment_in_All'], params['Retreatment_in_All']))
  cohort$NovelHist <- 0
  cohort$HIV <- sample(c(0,1), size = nrow(cohort), replace = T, 
                       prob=c(1-params['HIV_in_all'], params['HIV_in_all']))
  cohort$SmearStatus[cohort$HIV==0] <- sample(c(0, 1), size = sum(cohort$HIV==0), replace = T, 
                                                 prob=c(1-params['Smearpos_in_HIVneg'], params['Smearpos_in_HIVneg']))
  cohort$SmearStatus[cohort$HIV==1] <- sample(c(0, 1), size = sum(cohort$HIV==1), replace = T, 
                                              prob=c(1-params['Smearpos_in_HIVpos'], params['Smearpos_in_HIVpos']))
  cohort$RIF[cohort$RxHist==0] <- sample(c(0, 1), size = sum(cohort$RxHist==0), replace = T, 
                                         prob=c(1-params['RIF-R_in_New'], params['RIF-R_in_New']))
  cohort$RIF[cohort$RxHist==1] <- sample(c(0, 1), size = sum(cohort$RxHist==1), replace = T, 
                                         prob=c(1-params['RIF-R_in_Retreatment'], params['RIF-R_in_Retreatment']))
  cohort$INH[cohort$RxHist==0&cohort$RIF==0] <- sample(c(0, 1), size = sum(cohort$RxHist==0&cohort$RIF==0), replace = T, 
                                                       prob=c(1-params['INH-R_in_NewRIF-S'], params['INH-R_in_NewRIF-S']))
  cohort$INH[cohort$RxHist==1&cohort$RIF==0] <- sample(c(0, 1), size = sum(cohort$RxHist==1&cohort$RIF==0), replace = T, 
                                                       prob=c(1-params['INH-R_in_RetreatmentRIF-S'], params['INH-R_in_RetreatmentRIF-S']))
  cohort$INH[cohort$RIF==1] <- sample(c(0, 1), size = sum(cohort$RIF==1), replace = T, 
                                      prob=c(1-params['INH-R_in_RIF-R'], params['INH-R_in_RIF-R']))

# if modeling PZA-moxi correlation, will correlate only with MOXI (ignore partialmoxi). 
# Within each RR stratum, I'll have a proportion m of Moxi R vs S, and an overall proportion z of PZA-R I need to assign among the rr
# If I assign proportion a (oro the RR or RS) to be both z and m resistant, then OR = a(1-z-m+a)/(m-a)(z-a), and solving a quadratic for a:
  imputecombo <- function(OR, z, m)
  { C1 <- (OR-1); C2 <- (z+m-OR*z-OR*m-1); C3 <- OR*m*z
    a2 = (-C2 - sqrt(C2^2 - 4*C1*C3))/(2*C1)
    if (0 <= a2 & a2 <= 1)
      return(c(a2, z-a2, m-a2, 1-z-m+a2)) # z and m, then z not m, then m not z, then neither
    else (return(NA))
  }
  
  cohort$PZA <- cohort$MOXI <- NA
  cohort[cohort$RIF==1,c("PZA", "MOXI")] <- array(c(1,1,0,0,1,0,1,0), dim=c(4,2))[sample(1:4, size=sum(cohort$RIF==1), replace=T,
                                                                            prob=imputecombo(OR=params["OR-PZA-if_MOXI"], z=params['PZA-R_in_RIF-R'], m=params['MOXI-R-any_in_RIF-R']))
                                                                            ,]
  cohort[cohort$RIF==0,c("PZA", "MOXI")] <- array(c(1,1,0,0,1,0,1,0), dim=c(4,2))[sample(1:4, size=sum(cohort$RIF==0), replace=T,
                                                                                         prob=imputecombo(OR=params["OR-PZA-if_MOXI"], z=params['PZA-R_in_RIF-S'], m=params['MOXI-R-any_in_RIF-S']))
                                                                                  ,]

    cohort$partialmoxi[cohort$MOXI==1&cohort$RIF==1] <- sample(c(0, 1), size = sum(cohort$MOXI==1&cohort$RIF==1), replace = T, 
                                                             prob=c(params['MOXI-R-highlevel_in_RIF-R']/params['MOXI-R-any_in_RIF-R'], 
                                                                    1 - params['MOXI-R-highlevel_in_RIF-R']/params['MOXI-R-any_in_RIF-R']))
  cohort$partialmoxi[cohort$MOXI==1&cohort$RIF==0] <- sample(c(0, 1), size = sum(cohort$MOXI==1&cohort$RIF==0), replace = T, 
                                                             prob=c(params['MOXI-R-highlevel_in_RIF-S']/params['MOXI-R-any_in_RIF-S'], 
                                                                    1 - params['MOXI-R-highlevel_in_RIF-S']/params['MOXI-R-any_in_RIF-S']))
  cohort$partialmoxi[cohort$MOXI==0] <- 0
  
  # # here, MOXI indicates any moxi resistance, and partialmoxi indicates only low-level/partial resistance (and partial activity)
  # cohort$MOXI[cohort$RIF==1] <- sample(c(0, 1), size = sum(cohort$RIF==1), replace = T, 
  #                                      prob=c(1-params['MOXI-R-any_in_RIF-R'], params['MOXI-R-any_in_RIF-R']))
  # cohort$MOXI[cohort$RIF==0] <- sample(c(0, 1), size = sum(cohort$RIF==0), replace = T, 
  #                                      prob=c(1-params['MOXI-R-any_in_RIF-S'], params['MOXI-R-any_in_RIF-S']))
  # cohort$PZA[cohort$RIF==1& cohort$MOXI==1] <- sample(c(0, 1), size = sum(cohort$RIF==1), replace = T, 
  #                                     prob=c(1-params['PZA-R_in_RIF-R'], params['PZA-R_in_RIF-R']))
  # cohort$PZA[cohort$RIF==0] <- sample(c(0, 1), size = sum(cohort$RIF==0), replace = T, 
  #                                     prob=c(1-params['PZA-R_in_RIF-S'], params['PZA-R_in_RIF-S']))

  cohort$BDQ <- sample(c(0, 1), size = nrow(cohort), replace = T,
                       prob=c(1-params['BDQ-R_in_All'], params['BDQ-R_in_All']))
  cohort$PA <- sample(c(0, 1), size = nrow(cohort), replace = T, 
                      prob=c(1-params['PA-R_in_All'], params['PA-R_in_All']))
  
  return(cohort)
}



set.scenario <- function(scenario="0")
{
  # No new regimen (0)
  # New regimen considered for known RR only
  #   Without Xpert XDR (1a)
  #   With Xpert XDR if RR (and alternative regimen if FQ-R) (1b)
  # New regimen considered for all TB diagnoses
  #   Without up-front DST (2)
  #     Using shorter duration for everyone (2a)
  #     Using longer duration for everyone (2b)
  #   With up-front RIF DST only (3)
  #     Would use longer duration only if RIF-R detected
  #   With up-front step-wise RIF DST, then FQ DST if RIF-R (4)
  #     Would use longer duration if RIF-R or if no DST, FQ-S; alternative regimen if RIF-R, FQ-R
  #   With up-front DST for RIF and FQ for everyone (5)
  #     Would use longer duration if RIF-R or if no DST, FQ-S; alternative regimen if FQ-R (RIF-R or RIF-S) 
  ## need to define regimen_s, regimen_r, regimen_m, xpertxdr
  intervention <- list()
  intervention <- within(intervention,
                         {
                           # treatment focused analyses
                           if (scenario=="0") { regimen_s <- "hrze"; regimen_r <- "mdr"; xpertxdr <- "none"; allxpert<-FALSE }
                           if (scenario=="1a") { regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-FALSE }  
                           if (scenario=="1x") {regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-TRUE}
                           if (scenario=="3") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-FALSE }  
                            #possibly also including no Xpert at all approaches:
                           if (scenario=="2a") { regimen_s <- "bpamz4"; regimen_r <- "bpamz4"; xpertxdr <- "none"; allxpert<-FALSE  }  
                           if (scenario=="2b") { regimen_s <- "bpamz6"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-FALSE  }  

                           # DST focused analyses, universal regimens:
                           if (scenario=="3x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "none" ; allxpert<-TRUE }
                           if (scenario=="4x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "stepwise"; allxpert<-TRUE  }  
                           if (scenario=="5x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "simultaneous"; allxpert<-TRUE  }  
                           
                           # # and could consider DST question for MDR only (but too much complexity re baseline DST andconventional XDR regimens):
                           # if (scenario=="1a") { regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none" }  
                           # if (scenario=="1b") { regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "stepwise" }  
                           # 
                         })
  return(intervention)
}

prob <- function(p, N, stochasticmode=TRUE)
{
  if(stochasticmode) return(rbinom(N, 1, p)) else return(p)
}

Novelavailable <- function(patient, params, stochasticmode=TRUE)
{
  N <- nrow(patient)
  
  Considernovel <- numeric(N)
  
  Considernovel[1:N] <- prob(params["novelcoverage"], N, stochasticmode)
  
  return(Considernovel)
}

Xpertdone <- function(scenario, patient, params, stochasticmode=TRUE)
{
  if (missing(scenario)) scenario <- "0"
  intervention <- set.scenario(scenario)
  if(intervention$allxpert) params["Xpert_current_new"] <- params["Xpert_current_rerx"] <- params["Xpert_max"]
  
  N <- nrow(patient)
  
  # probability of Xpert DST:
  
  Xpert <- (patient$RxHist==1) * prob( params["Xpert_current_rerx"], N, stochasticmode)  + 
            (patient$RxHist==0) * prob( params["Xpert_current_new"], N, stochasticmode) 
  return(Xpert)
}

# get TB and RIFR results of Xpert (will be zero if Xpert not done, or if Xpert negative) -- so Xpert negative is Xpert - Xpertresult. 
Xpertresult <- function(Xpert, patient, params, stochasticmode=TRUE)
{
  N <- nrow(patient)
  
  Result <- array(NA, dim=c(N, 2))
  
  # Xpert shows TB: 
  Result[,1] <- Xpert*
    ((patient$SmearStatus==1) * prob( params["Sens_TB_smearpos"], N, stochasticmode)  + 
    (patient$SmearStatus==0) * prob( params["Sens_TB_smearneg"], N, stochasticmode) )
  Result[,2] <- Result[,1]*prob(p = (patient$RIF=="1")*params["Sens_RR"], N, stochasticmode)  
  
  return(Result)
}


allDSTresults <- function(scenario, Result, patient, params, stochasticmode=TRUE) #asuming they are considered for novel intervention, and assuming Xpert is done and detects TB
{
  if (missing(scenario)) scenario <- "0"
  intervention <- set.scenario(scenario)
  
  N <- nrow(patient)
  
  FullDST <- array(NA, dim=c(N, 4)); colnames(FullDST) = c("none", "RR", "FQR", "both")
  
  # probability of detecting FQR if Xpert XDR is done (for now, we won't have treatment selection depend on type of FQ resistance, but only treatment outcome):
  fqdetect <- (patient$partialmoxi==1)*prob(params["Sens_FQ_low"], N, stochasticmode) + (patient$MOXI==1&patient$partialmoxi==0)*prob(params["Sens_FQ_high"], N, stochasticmode)
  
  FullDST[,4] <- Result[,2]*(intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")*fqdetect
  
  FullDST[,2] <- Result[,2]*( (intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")*(1 - fqdetect) + 
                                (1-(intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")) )
  
  FullDST[,3] <- (1-Result[,2])*(intervention$xpertxdr=="simultaneous")*fqdetect
  
  FullDST[,1] <- 1-rowSums(FullDST[,2:4])
  
  return(FullDST)  
}

Regimenselect <- function(scenario, Result, Considernovel, FullDST, patient, stochasticmode) # take DST results + scenario + patient history, and choose regimen
  ##** need to modify array handling to allow just a single patient
{
  if (missing(scenario)) scenario <- "0"
  intervention <- set.scenario(scenario)
  
  regs <- array(0, dim=c(nrow(patient), length(regimens))); colnames(regs) <- regimens

  # Consider novel
  # Detect FQ and RR
  regs[,'bpal'] <- regs[,'bpal'] + Considernovel*FullDST[,"both"]*(1-patient$NovelHist)
  regs[,'backup_mdr'] <- regs[,'backup_mdr'] + Considernovel*FullDST[,"both"]*patient$NovelHist
  # Detect FQ only
  regs[,'bpal'] <- regs[,'bpal'] + Considernovel*FullDST[,"FQR"]*(1-patient$NovelHist) 
  regs[,'backup'] <- regs[,'backup'] + Considernovel*FullDST[,"FQR"]*(patient$NovelHist)
  # Detect RR only
  regs[,intervention$regimen_r] <- regs[,intervention$regimen_r] + Considernovel*FullDST[,"RR"]
  # Detect neither
  regs[,intervention$regimen_s] <-regs[,intervention$regimen_s] + Considernovel*FullDST[,"none"]
  
  # Don't consider novel
  # Detect RR
  regs[,'mdr'] <- regs[,'mdr'] + (1-Considernovel)*Result[,2]
  # Don't detect RR
  regs[,'hrze'] <- regs[,'hrze'] + (1-Considernovel)*(1-Result[,2])
  
  if(stochasticmode) return(colnames(regs)[unlist(apply(regs, 1, function(x) which(x==1)))]) else  return(regs)
}

# Treatment outcome functions:
wallis <- function(cxconv, months, coefs, params)
{
  if(missing(coefs)) {coefs <- params[c("coef1","coef2","coef3")]}
  return(exp(coefs["coef1"]-coefs["coef2"]*log(months)+coefs["coef3"]*log((1-cxconv)/(cxconv)))/(1+exp(coefs["coef1"]-coefs["coef2"]*log(months)+coefs["coef3"]*log((1-cxconv)/(cxconv))))) 
}
reversewallis <- function(recurrence, months, coefs, params) # gets 2 month cxconv from recurrence after n months of treatment
{
  if(missing(coefs)) {coefs <- params[c("coef1","coef2","coef3")]}
  return(1-invlogit((logit(recurrence)-coefs["coef1"]+coefs["coef2"]*log(months))/coefs["coef3"]))
}

# # predict (i.e. randomly assign based on an adherence model) their adherence; if the regimen duration is <= this, they'll complete it
# monthsadherent <- function(params, N=1)
# {  return(sample(monthsmodeled, size = N, prob = params['Monthlyloss']*(1-params['Monthlyloss'])^(monthsmodeled-1), replace=T)) 
# }
# will be included elsewhere

make.active.regimen.matrix <- function() 
{
  # Identify the "active regimen" for each combination of patient and regimen
  resistance <- expand.grid( rep( list( 0:1), length(drugs))) 
  colnames(resistance) <- drugs
  resistance <- resistance %>% filter(MOXI==1 | partialmoxi==0)
  active <- cbind(resistance, array(NA, dim=c(nrow(resistance), length(regimens))))
  colnames(active) <- c(drugs, regimens)
  active[,"hrze"] <- ifelse(active[,"RIF"]==1, ifelse(active[,"INH"]==1, "(ZE)", "H(ZE)"), ifelse(active[,"INH"]==1, "R(ZE)", "HR(ZE)"))
  active[,"mdr"] <- ifelse(active[,"MOXI"]==1, ifelse(active[,"partialmoxi"],"MDR, FQ-low", "MDR, FQ-R"), "MDR, FQ-S")
  active[,"bpamz4"] <- active[,"bpamz6"] <- paste0(ifelse(active[,"BDQ"]==1, "", "B"),ifelse(active[,"PA"]==1, "", "Pa"), 
                                                   ifelse(active[,"MOXI"]==1, ifelse(active[,"partialmoxi"]==1, "m", ""), "M"), 
                                                   ifelse(active[,"PZA"]==1, "", "Z") )
  active[,"bpal"] <- paste0(ifelse(active[,"BDQ"]==1, "", "B"),ifelse(active[,"PA"]==1, "", "Pa"), "L")
  active[,"backup"] <- active[,"hrze"]
  active[,"backup_mdr"] <- active[,"mdr"]
  
  
  active["resistance"] <- data.matrix(active[,1:length(drugs)])%*%10^((length(drugs)-1):0)
  return(active)
}

activeregimen <- function(matrix, patient, Regimen)
{
  if(missing(matrix)) matrix <- make.active.regimen.matrix()
  
  patient$resistance <- data.matrix(patient[,drugs])%*%10^((length(drugs)-1):0)
  
  rows <- apply(patient$resistance, 1, function(x) which(matrix$resistance==x))
  Regimenarray <- array(Regimen, dim=c(length(Regimen),1))
  cols <- apply(Regimenarray, 1, function(x) which(colnames(matrix)==x))
  
  return(matrix[cbind(rows, cols)])
}


make.recurrence.matrix <- function()
{
  # Model treatment outcomes (probablities of cure, failure, relapse, and acquired resistance, for a given phenotype and a given months completed)
  recurrencematrix <- array(NA, dim=c(length(activeregimens),length(monthsmodeled))); dimnames(recurrencematrix) <- list(activeregimens,monthsmodeled)

    recurrencematrix["HR(ZE)",] <- wallis(cxconv = reversewallis(params["HRZE_pooled_relapse"]/(params["INH_multiplier"]*params["pooled_INH_fraction"] + (1-params["pooled_INH_fraction"])),6, params=params), months = monthsmodeled, params=params)
  # HR(ZE)_relapse * (INH_multiplier*pooled_INH_fraction + (1-pooled_INH_fraction)) = HRZE_pooled_relapse
  recurrencematrix["R(ZE)",] <- params["INH_multiplier"]*recurrencematrix["HR(ZE)",] 
  recurrencematrix["BPaMZ",] <- wallis(reversewallis(params["HRZE_pooled_relapse"],4,params=params), monthsmodeled, params=params)
  recurrencematrix["BPaM",] <- wallis(params["BPaM_cxconv"], monthsmodeled,params=params)
  recurrencematrix["BPaZ",] <- wallis(params["BPaZ_cxconv"], monthsmodeled,params=params)
  recurrencematrix["BMZ",] <- recurrencematrix["PaMZ",] <- recurrencematrix["BPaZ",]
  recurrencematrix["MDR, FQ-S",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-S"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/3, params = params)
  recurrencematrix["MDR, FQ-R",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-R"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/3, params = params)
  recurrencematrix["H(ZE)",] <- wallis(reversewallis(params["INHmono_failrelapse"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/1.5, params = params)
  recurrencematrix["BPa",] <- recurrencematrix["BZ",] <- recurrencematrix["PaZ",] <- recurrencematrix["BM",] <- recurrencematrix["PaM",] <- recurrencematrix["MZ",] <- recurrencematrix["R(ZE)",] 
  
  recurrencematrix["(ZE)",] <- recurrencematrix["B",] <- recurrencematrix["Pa",] <- recurrencematrix["M",] <- recurrencematrix["Z",] <- 
    c(wallis(cxconv = reversewallis(recurrence = params["Highrecurrence"],months = 6,params = params),months = monthsmodeled,params = params)[1:6], rep(params["Highrecurrence"], length(monthsmodeled)-6))
  
  recurrencematrix["BPaL",] <- recurrencematrix["HR(ZE)",]
  recurrencematrix["BL",] <- recurrencematrix["PaL",] <- recurrencematrix["BPa",]
  recurrencematrix["m",] <- recurrencematrix["L",] <- recurrencematrix["none",] <- recurrencematrix["B",]
  
  recurrencematrix[recurrencematrix>1] <- 1
  
  recurrencematrix["BPamZ",] <- increaseodds(recurrencematrix["BPaZ",], params["partialmoxiOR"])
  recurrencematrix["MDR, FQ-low",] <- increaseodds(recurrencematrix["MDR, FQ-R",], params["partialmoxiOR"])
  recurrencematrix["BPam",] <- increaseodds(recurrencematrix["BPa",], params["partialmoxiOR"])
  recurrencematrix["BmZ",] <- increaseodds(recurrencematrix["BZ",], params["partialmoxiOR"])
  recurrencematrix["PamZ",] <- increaseodds(recurrencematrix["PaZ",], params["partialmoxiOR"])
  recurrencematrix["Bm",] <- increaseodds(recurrencematrix["B",], params["partialmoxiOR"])
  recurrencematrix["Pam",] <- increaseodds(recurrencematrix["Pa",], params["partialmoxiOR"])
  recurrencematrix["mZ",] <- increaseodds(recurrencematrix["Z",], params["partialmoxiOR"])
  
  recurrencematrix[recurrencematrix>params["Highrecurrence"]] <- params["Highrecurrence"]
  
  return(recurrencematrix)
} 

make.adr.matrix <- function()
{
  adrmatrix <- array(NA, dim=c(length(activeregimens), 5)); dimnames(adrmatrix) <- list(activeregimens, c("HR", "RR", "FQR", "BR", "PaR"))
  adrmatrix[,"RR"] <- c(params["adr_r"]*c(1, params["adrfactor_other"]), rep(0, length(activeregimens)-2))
  adrmatrix["MDR, FQ-S","FQR"] <- adrmatrix["MDR, FQ-low","FQR"] <- params["adr_mdr"]
  adrmatrix["BPaMZ","FQR"] <- adrmatrix["BPamZ","FQR"] <- params["adr_bpamz"]
  adrmatrix["BPaM","FQR"] <- adrmatrix["BPam","FQR"] <- params["adr_bpamz"]*params["adrfactor_z"]
  adrmatrix["BMZ","FQR"] <- adrmatrix["PaMZ","FQR"] <- adrmatrix["BmZ","FQR"] <- adrmatrix["PamZ","FQR"] <- params["adr_bpamz"]*params["adrfactor_other"]
  adrmatrix["BM","FQR"] <- adrmatrix["PaM","FQR"] <- adrmatrix["Bm","FQR"] <- adrmatrix["Pam","FQR"] <- params["adr_bpamz"]*params["adrfactor_other"]*params["adrfactor_z"]
  adrmatrix["MZ","FQR"] <- adrmatrix["mZ","FQR"] <- params["adr_bpamz"]*params["adrfactor_twodrugs"]
  adrmatrix["M","FQR"] <- adrmatrix["m","FQR"] <- params["Highrecurrence"]*(1+params[ "Failures_per_recurrence"])/(1+params["Highrecurrence"]*params[ "Failures_per_recurrence"])
  adrmatrix[c("BPaMZ", "BPaM", "BPaZ", "BMZ", "BPa", "BZ", "BM", "BPaL", "BL"),"BR"] <- 
    params["adr_bpamz"]*c(1, params["adrfactor_z"], params["adrfactor_other"], params["adrfactor_other"], params["adrfactor_other"]*params["adrfactor_z"], 
                          params["adrfactor_twodrugs"], params["adrfactor_other"]*params["adrfactor_z"], 1, params["adrfactor_other"])
  adrmatrix[c("BPaMZ", "BPaM", "BPaZ", "PaMZ", "BPa", "PaZ", "PaM", "BPaL", "PaL"),"PaR"] <- 
    params["adr_bpamz"]*c(1, params["adrfactor_z"], params["adrfactor_other"], params["adrfactor_other"], params["adrfactor_other"]*params["adrfactor_z"], 
                          params["adrfactor_twodrugs"], params["adrfactor_other"]*params["adrfactor_z"], 1, params["adrfactor_other"])
  adrmatrix["H(ZE)","HR"] <- adrmatrix["B", "BR"] <- adrmatrix["Pa", "PaR"] <- params["Highrecurrence"]*(1+params[ "Failures_per_recurrence"])/(1+params["Highrecurrence"]*params[ "Failures_per_recurrence"])
  adrmatrix[c("BPamZ", "BPam", "BmZ", "Bm"),"BR"] <- 
    params["adrfactor_partialmoxi"]*adrmatrix[c("BPaZ", "BPa", "BZ", "B"),"BR"] 
  adrmatrix[c("BPamZ", "BPam", "PamZ", "Pam"),"PaR"] <- 
    params["adrfactor_partialmoxi"]*adrmatrix[c("BPaZ", "BPa", "PaZ", "Pa"),"PaR"] 

  return(adrmatrix)
}

###########################
# Steps for running model #

## FIRST EVENT ##
# TB diagnosis, not with tbdxtime from an exponential distribution
# ignoring here any who will spontaneously resolve before diagnosis or death, and assuming that those who fail or felapse won't spontaneously resolve (or won't be inncluded among documented relapse rates)
diagnosisevent <- function(last, params, N, eventtypes, statetypes, regimentypes)
{
  current <- last
  
  # ** could later change the diagnosis time distribution (and same for mortallity etc?) to something less skewed
  needdiagnosis <- ((last[,"TBstate"]==statetypes$undiagnosed)|(last[,"TBstate"]==statetypes$failed)| (last[,"TBstate"]==statetypes$relapsed))
  dxtime <- rexp(N, 1/(last[,"RxHist"]*params["Tbdxtime_recurrence"] + (1-last[,"RxHist"])*params["Tbdxtime"]))
  
  # check for mortality before initial TB diagnosis event:
  morttime <- rexp(N, params["Monthlymortality_TB"] + last[,"HIV"]*params["Monthlymortality_HIV"])
  deaths <- needdiagnosis&(morttime < dxtime  )
  current[deaths, "eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]; 
  current[deaths, "eventtype"] <- eventtypes$death; current[deaths, "TBstate"] <- statetypes$deceased; # leave others as diagnosis and diagnosed
  
  current[needdiagnosis&(!deaths), "eventtime"] <- last[needdiagnosis&(!deaths), "eventtime"] +  dxtime[needdiagnosis&(!deaths)]; 
  current[needdiagnosis&(!deaths),"eventtype"] <- eventtypes$TBdiagnosis; current[needdiagnosis&(!deaths),"TBstate"] <- statetypes$diagnosed
  
  carryforward <- !needdiagnosis
  current[carryforward,"eventtype"] <- eventtypes$carryforward
  
  return(current)
}

## SECOND EVENT ##
## Treatment initiation attempt (will account for any DST, regimen selection, associated delays, possibility of pretreatment LTFU, and possibility of ADR)
# assumes we already saved current to course
treatmentinitiationattempt <- function(last, params, N, eventtypes, statetypes, regimentypes, scenario)
{
  intervention <- set.scenario(scenario)
  
  current <- last
  currentcohort <- data.frame(current[, c("RxHist", "NovelHist", "HIV","SmearStatus","RIF", "INH", "PZA", "MOXI","partialmoxi", "BDQ", "PA")])
  
  # perform any DST and select a regimen
  Considernovel <- Novelavailable(currentcohort, params, stochasticmode)
  Xpert <- Xpertdone(scenario, currentcohort, params, stochasticmode) # is xpert done? depends on scenario and patient
  Result <- Xpertresult(Xpert, currentcohort, params, stochasticmode) # Xpert MTB/RIF result. Dependent on Xpertdone (And therefore indirectly on scenario)
  FullDST <- allDSTresults(scenario, Result, currentcohort, params, stochasticmode) # Full results of DST, depending on Xpert MTB/RIF result (indirectly on scenario and patient) and on Xpert XDR use scenario
  Regimen <- Regimenselect(scenario, Result, Considernovel, FullDST, currentcohort, stochasticmode) # regimen selected for each patient, depending on the DST results, the novel availability, the patient's prior novel regimen exposure, and the scenario's preferred regimens
  
  if(!stochasticmode) stop("Can't use the matrix regimen output in the lines that follow, still need to code this version")
  
  # event will now be treatment start, but we'll wait to add that only for those who start treatment
  # determine time of this event
  # DST, and selection of a nonstandard regimen, may incur delays or pretreatment losses to follow up
  needtreatment <- (last[,"TBstate"] == statetypes$diagnosed)
  
  starttime <- Xpert*params["DSTdelay"] 
  # delays for non-first-line regimens:
  if (intervention$regimen_s=="hrze") starttime <- starttime + (Regimen!=intervention$regimen_s)*params["Regimendelay"]
  if (intervention$regimen_s %in% c("bpamz4", "bpamz6")) starttime <- starttime + (!(Regimen %in% c("bpamz4", "bpamz6")))*params["Regimendelay"]
  morttime <- rexp(N, params["Monthlymortality_TB"] + last[,"HIV"]*params["Monthlymortality_HIV"])
  
  deaths <-  needtreatment&(morttime < starttime)
  current[deaths, "eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]; 
  current[deaths, "eventtype"] <- eventtypes$death
  current[deaths, "TBstate"] <- statetypes$deceased; # leave others as assigned above (treated or lost)
  
  lost <- needtreatment&(!deaths)&(prob(params["Unavoidableloss"] + Xpert*params["DSTloss"] + (Regimen!=intervention$regimen_s)*params["Regimenloss"], N, stochasticmode))
  current[lost,"eventtype"] <- eventtypes$pretreatmentltfu
  current[lost,"TBstate"] <- statetypes$undiagnosed
  
  current[needtreatment&(!deaths), "eventtime"] <- last[needtreatment&(!deaths), "eventtime"] + starttime[needtreatment&(!deaths)];
  
  treated <- needtreatment&(!deaths)&(!lost)
  current[treated,"eventtype"] <- eventtypes$treatmentstart
  
  #those who aren't lost start treatment, possibly with ADR
  # assign regimen
  current[treated, "Currentregimen"] <- as.numeric(regimentypes[Regimen])[treated]
  active <- activeregimen(make.active.regimen.matrix(), currentcohort, Regimen)
  active[active==""] <- "none"
  ## for each patient, sample risks of ADR based on active regimen, and if any occur, change phenotype in cohort and set state as _adr
  adrrisks <- make.adr.matrix()[active,]
  adr <- array(rbinom(N*4, 1, adrrisks[]), dim=c(N,4)); adr[is.na(adr)] <- 0
  # update state accordingly (haven't removed deaths yet)
  current[treated, "TBstate"] <- ((rowSums(adr)>=1)*statetypes$treating_adr + (rowSums(adr)<1)*statetypes$treating)[treated]
  # update resistance accordingly, if ADR: change 0s (for both FQ columns) in current to 1 if adr is 1 (and !lost) 
  current[treated, c("RIF", "MOXI", "BDQ", "PA")] <- pmax(last[treated, c("RIF", "MOXI", "BDQ", "PA")], adr[treated,])
  # if moxi adr but already were moxi, then partialmoxi needs to change to zero
  # (but first make sure we're allowing moxi adr for baseline partialmoxi -- yes, because the "m"s have a risk of FQ-R ADR in adrmatrix.)
  current[(treated&(adr[,2]==1)&(last[,"partialmoxi"]==1)), "partialmoxi"] <- 0
  
    
  current[!needtreatment,"eventtype"] <- eventtypes$carryforward
  
  # save before moving to next event
  
  return(current)
}


## THIRD EVENT ###
## treatment end (due to end of regimen or loss to follow up)
## will carry forward those who aren't dead or getting treated

treatmentend <- function(last, params, N, eventtypes, statetypes, regimentypes)
{
  current <- last
  
  onRx <- (last[,"TBstate"] == statetypes$treating | last[,"TBstate"] == statetypes$treating_adr); 
  
  # decide when their regimen should end, when/if they'll be lost, and when/if they'll die
  endtime <- numeric(N); endtime[onRx] <- regimendurations[current[,"Currentregimen"][onRx]]; endtime[!onRx] <- NA
  ltfutime <- rexp(N, params["Monthlyloss"])
  morttime <- rexp(N, params["Monthlymortality_background"] + last[,"HIV"]*params["Monthlymortality_HIV"]) #**for now, assuming only background mortality; could make this dependend on current failure probability, but this is okay as long as the pre-treatment mortality is high enough to account for TB mortality during early treatment
  
  deaths <- onRx&(morttime < pmin(endtime, ltfutime))
  current[deaths, "eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]; 
  current[deaths, "eventtype"] <- eventtypes$death
  current[deaths, "TBstate"] <- statetypes$deceased
  current[deaths, "Currentregimen"] <- 0
  
  # all who ADR become by definition relapses or failures
  failrelapseadr <- (!deaths)&(last[,"TBstate"]==statetypes$treating_adr) # those who have acquired resistance during current treatment and don't die during current treatment
  fail_if_adr <- rbinom(N,1, prob = params['Failures_per_recurrence']/(1+params['Failures_per_recurrence'])) #"if they are in failrelapseadr, will they fail as opposed to relapse?"
  
  current[failrelapseadr, "eventtime"] <- last[failrelapseadr, "eventtime"] + pmin(endtime, ltfutime)[failrelapseadr]; 
  current[failrelapseadr, "eventtype"] <- eventtypes$treatmentend
  current[failrelapseadr,"TBstate"] <- fail_if_adr[failrelapseadr]*statetypes$failed + ((!fail_if_adr)[failrelapseadr]*statetypes$pendingrelapse)
  current[failrelapseadr,"RxHist"] <- 1
  current[failrelapseadr & (last[,"Currentregimen"] %in% which(regimens %in% c("bpamz4", "bpamz6", "bpal"))), "NovelHist"] <- 1
  current[failrelapseadr,"Currentregimen"] <- 0
  current[failrelapseadr,"SmearStatus"] <- prob(current[failrelapseadr,"HIV"]*params['Smearpos_in_HIVpos'] + 
                                                  (1-current[failrelapseadr,"HIV"])*params['Smearpos_in_HIVneg'],
                                                N=sum(failrelapseadr), T) 

  # the rest get an outcome assigned when they stop treatment at end of regimen or at ltfu time
  needoutcome <- onRx&(!failrelapseadr)&(!deaths)
  
  current[needoutcome, "eventtime"] <- last[needoutcome, "eventtime"] + pmin(endtime, ltfutime)[needoutcome]; 
  current[needoutcome, "eventtype"] <- eventtypes$treatmentend
  # decide failure or relapser TBstate:
  ## adjust remaining risk of relapse/failure based on overall risk of adr
  adrrisks <- make.adr.matrix()
  adrrisks[is.na(adrrisks)] <- 0
  totaladrrisks <- adrrisks[,1] + (1-adrrisks[,1])*(adrrisks[,2] + (1-adrrisks[,2])*(adrrisks[,3] + (1-adrrisks[,3])*adrrisks[,4]))
  # if recurrences are R of those at risk, and failures are p*R (p=failures_per_Rec), 
  # then overall failures or recurrences will be f + (1-f)*R = (p+1)*(1-f)*R --> f=pR/(1+pR) -->f+(1-f)R=R(1+p)/(1+pR)
  # so then really, the cap on poor outcomes is a little over highrecurrence, at 83% rather than 80, but I think that's okay.
  # and it will be higher still when ADR=1, so I'll reduce the max ADR to highrecurrence (p(1+R)/(1+pR) to allow for self-cure there as well. 
  failures_or_relapses_remaining <- apply(make.recurrence.matrix()*(1+params["Failures_per_recurrence"])/(1+make.recurrence.matrix()*params["Failures_per_recurrence"]), 2, "-", totaladrrisks)
  failures_or_relapses_remaining[failures_or_relapses_remaining<0] <- 0
  ## assign relapse/failure based on (rounded) time completed:
  Rxtime <- pmin(endtime, monthsmodeled[unlist(lapply(floor(ltfutime), function(x) which.min(abs(monthsmodeled-x))))])
  ## these three will have the length of needoutcome
  active <- activeregimen(make.active.regimen.matrix(), data.frame(current[needoutcome,,drop=FALSE]), regimens[current[needoutcome,"Currentregimen"]])
  active[active==""] <- "none"
  cures <- rbinom(n=sum(needoutcome), size=1, prob = 1 - failures_or_relapses_remaining[cbind(active,as.character(Rxtime[needoutcome]))])
  relapses <- (!cures)&rbinom(sum(needoutcome), 1, 1/(1+params["Failures_per_recurrence"]))
  ## assign:
  current[needoutcome, "TBstate"][cures==1] <- statetypes$cured
  current[needoutcome, "TBstate"][relapses==1] <- statetypes$pendingrelapse
  current[needoutcome, "TBstate"][(!cures)&(!relapses)] <- statetypes$failed
  current[needoutcome,"RxHist"] <- 1
  current[needoutcome & (last[,"Currentregimen"] %in% which(regimens %in% c("bpamz4", "bpamz6", "bpal"))), "NovelHist"] <- 1
  current[needoutcome,"Currentregimen"] <- 0
  
  current[!onRx,"eventtype"] <- eventtypes$carryforward
  
  return(current)
}

## EVENT 4 = RELAPSES (vs other mortality while pending relapse), and carry others forward if they're not dead ## 
relapseevent <-  function(last, params, N, eventtypes, statetypes, regimentypes)
{
  current <- last
  
  relapsers <- last[,"TBstate"]==statetypes$pendingrelapse
  
  relapsetime <- rexp(N, 1/params["Tbdxtime_recurrence"]) 
  morttime <- rexp(N, params["Monthlymortality_background"] + last[,"HIV"]*params["Monthlymortality_HIV"]) 
  
  deaths <- relapsers&(morttime < relapsetime)
  current[deaths, "eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]; 
  current[deaths, "eventtype"] <- eventtypes$death; current[deaths, "TBstate"] <- statetypes$deceased
  
  current[relapsers&(!deaths), "eventtime"] <- last[relapsers&(!deaths), "eventtime"] + relapsetime[relapsers&(!deaths)]; 
  current[relapsers&(!deaths), "eventtype"] <- eventtypes$TBonset
  current[relapsers&(!deaths), "TBstate"] <- statetypes$relapsed
  
  current[!relapsers,"eventtype"] <- eventtypes$carryforward
  
  return(current)
}

# add mortality during intervals without TB
#(already included above for those pending relapse, and for those with TB)
addnaturalmortality <- function(last, params, N, eventtypes, statetypes, maxtime=1000*12)
{
  current <- last
  
  cures <- last[,"TBstate"]==statetypes$cured
  
  morttime <- rexp(N, params["Monthlymortality_background"] + last[,"HIV"]*params["Monthlymortality_HIV"]) 
  deaths <- cures & (morttime<maxtime)
    
  current[deaths,"eventtype"] <- eventtypes$death
  current[deaths,"eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]
  current[deaths, "TBstate"] <- statetypes$deceased
  
  current[!deaths,"eventtype"] <- eventtypes$carryforward
  
  return(current)
}

# time.in.state <- function (states, course, characteristics=c(), cutofftime=100*12, carryforward=F) 
# {
#   elapsedtimes <- course[,"eventtime",2:((dim(course))[[3]])] - course[,"eventtime",1:(((dim(course))[[3]])-1)] 
#   elapsedtimes[course[,"eventtime",2:((dim(course))[[3]])]>cutofftime] <- cutofftime - (course[,"eventtime",1:(((dim(course))[[3]])-1)])[course[,"eventtime",2:((dim(course))[[3]])]>cutofftime]
#   elapsedtimes[elapsedtimes<0] <- 0
#   t <- rep(0, dim(course)[[1]])
#   for (state in states[1:length(states)]){
#     indices <-  course[,"TBstate",1:((dim(course))[[3]]-1)] == statetypes[[state]]
#     if(length(characteristics)==1) indices <- indices & course[,characteristics,1:((dim(course))[[3]]-1)] >0
#     if(length(characteristics)>1) indices <- indices & apply(course[,characteristics,1:((dim(course))[[3]]-1)] >0, 1, all)
#     t <- t + rowSums(elapsedtimes*indices)
#     if (carryforward) {
#       laststep <- course[,"TBstate",dim(course)[[3]]] == statetypes[[state]]
#       extratime <- cutofftime-course[laststep,"eventtime",dim(course)[[3]]]
#       extratime[extratime<0] <- 0
#       t[laststep] <- t[laststep] + extratime
#     }
#   }
#   return( summary(t))
# }


# function to model a set of patients through a TB episode
require(abind)
require(rlist)
require(prodlim)

modelcourse <- function(scenario="0", cohort, params, reps=1, steplimit=16, stochasticmode=TRUE) # will need to repeat same for each rep of a given patient, if probabilistic events
{
  repcohort <- do.call(rbind, replicate(reps, cohort, simplify=FALSE))
  intervention <- set.scenario(scenario)
  N=nrow(repcohort)
  
  # Structure this as a third dimension of cohort (or second dimension of patient) array, with one row per event date
  
  trackeditems <- c("eventtime", "eventtype", "TBstate", "Currentregimen", colnames(cohort))
  
  current <- array(dim=c(nrow(repcohort), length(trackeditems))); colnames(current) <- trackeditems
  course <- array(NA, dim=c(dim(current),1)); dimnames(course) <- list("patient"=1:nrow(repcohort), "trackeditem"=trackeditems, "event"=numeric())
  
  # TB onset for all patients
  course[,,1] <- current[] <- cbind(0, eventtypes$TBonset, statetypes$undiagnosed, regimentypes$none, data.matrix(repcohort))
  
  # for each separate patient, determine their course: 
  # (act on current event, save current to course when move to next event)
  course <- abind(course, diagnosisevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  course <- abind(course, treatmentinitiationattempt(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, scenario), along=3)
  course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)

  ## and back to event 5= DIAGNOSIS AGAIN ## 
  course <- abind(course, diagnosisevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  course <- abind(course, treatmentinitiationattempt(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, scenario), along=3)
  course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  
  #continue to loop up to sooner of steplimit/4 loops or no more (uncured or dead) or those still alive and uncured have run for 20 years:
  while(dim(course)[3]<steplimit & 
        mean(course[,"eventtype",dim(course)[[3]]] %in% c(statetypes$deceased, statetypes$cured))<1 & 
        mean(course[,"eventtime",dim(course)[[3]]][!course[,"eventtype",dim(course)[[3]]] %in% c(statetypes$deceased, statetypes$cured)] > 20*12)<1)
  {
    course <- abind(course, diagnosisevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
    course <- abind(course, treatmentinitiationattempt(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, scenario), along=3)
    course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
    course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  }
    
  # course <- abind(course, addnaturalmortality(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, maxtime = max(course[,"eventtime",])), along=3)
  course <- abind(course, addnaturalmortality(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes), along=3)
  
  newcourse <- array(course, dim=c(dim(course)[1]/reps,reps,dim(course)[2:3]))
  newcourse <- aperm(newcourse, c(1,3,4,2))
  dimnames(newcourse) <- list("patient"=c(), "characteristic"=trackeditems,"timestep"=c(),"rep"=c())
  
  return(newcourse)
}

## ** this model is assuming no correlation between moxi and pza resistance in RR, which is probably unrealistic 
# but which is optimistic re: BPaMZ impact and re: importance of Fq DST (bc FQ-R could also be predictive of Z resistance).
# Could add a parameter, correlation between moxi and pza resistance, and incorporate it into make.cohort(), but would be complex with low vs high resistance and adr calibration. 

# set the model type
stochasticmode <- T



## ** sensitivity analysis to add: correlate outcomes with Xpert status (Xpert negative shortens required treatment, Xpert positive extends it), \\\
# made patients Sm+ or Sm- (which will depend on HIV status), and will differentiate Xpert sensitivity (but not DST result yield once Xpert is positive) and, later, infectiousness, based on smear. Status will be reassigned at time of failure/relapse. 

## could consider correlation of resistance with other factors that promote relapse (but dont have data beyond what's already be captured in the data on resistance outcomes?)

# single patient model:
#  time of TB onset = t0
#  determine timing of TB diagnosis (can later make this and other dependent on smear status)
#  determine DST done and results if any (per algorithms, which may result in delay or loss)
#  determine regimen chosen
#  determine treatment start date (may be none, return to active pool if pretreatment loss)
#  determine acq res (will change susceptibilities, and mark as adr)
#  determine monthly probs of failure if adherent (may later tie this to infectiousness and/or mortality)
#  determine how many months will be completed if survives
#  determine outcome of cure/relapse/failure, with outcome date
#  update Rxhist and novelhist
#  determine date of relapse if occurs
#  determine date of next TB diagnosis if relapse or failure
#  determine whether they died before any of these events

# Notes:
# Because of uncertainties about infectiousness during treatment and how different regimens affect it (whether/how much it correlates with outcomes), our infectiousness calculations won't include time on treatment
# Failures will become immediately infectious when they stop treatment (if adding smear, P(Xpert+)etc, will start smear neg)
# ADR (acquired drug resistance) is a fraction of all treatments, so it will happen at the start of treatment,
# Those who ADR  will be destined to fail or relapse (divided into same ratio, and with same temporal dynamics after treatment ends) no matter how much treatment they get.
# When TB diagnosis occurs, any DST and regimen selection are modeled immediately, but treatment start may not be immediate
