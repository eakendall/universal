require(arm)
require(dplyr)

allparams <- read.csv("allparams.csv", header=T, stringsAsFactors = F)
params <- as.numeric(allparams[,2])
if(setting=="SAf") params[!is.na(allparams[,5]) & allparams[,5]!=""] <- as.numeric(allparams[,5][!is.na(allparams[,5]) & allparams[,5]!=""])
names(params) <- allparams[,1]
params["Tbdxtime_recurrence"] <- params["Tbdxtime_recurrenceratio"]*params["Tbdxtime"]


increaseodds <- function(recurrenceprobs, OR) # increase odds of success by OR
{
  initialodds <- (1-recurrenceprobs)/recurrenceprobs
  newodds <- OR*initialodds
  newrecurrenceprobs <- 1/(1+newodds)
  return(newrecurrenceprobs)
}

splitodds <- function(recurrenceprobsa, recurrenceprobsb)
{
  recurrenceprobsa[recurrenceprobsa>0.99] <- 0.99
  recurrenceprobsb[recurrenceprobsb>0.99] <- 0.99
  # average the logit of recurrence
  meanlogit <-( log((recurrenceprobsa)/(1-recurrenceprobsa)) + log((recurrenceprobsb)/(1-recurrenceprobsb)) )/2
  # and invert:
  return(exp(meanlogit)/(1+exp(meanlogit)))
}

drugs <- c("RIF", "INH", "PZA", "MOXI", "partialmoxi","BDQ", "PA")
adrs <- c("RIF", "INH", "MOXI", "BDQ", "PA")

regimens <- c("hrze", "mdr", "bpamz4", "bpamz4r", "bpamz6", "backup")
regimendurations <- c(6, 18, 4, 4, 6, 6)
# regimendurations <- c(6, 18, 4, 6, 6, 6, 6, 18)
activeregimens <- c("HR(ZE)","R(ZE)","H(ZE)","(ZE)","MDR, FQ-S","MDR, FQ-low","MDR, FQ-R","BPaMZ","BPaM","BPaZ","BPamZ", "BMZ","PaMZ",
                    "BPa","BPam", "BZ","BmZ", "BM","PaZ","PamZ", "PaM","MZ","B","Bm", "Pa","Pam","M","Z","mZ","m","none")
patientvars <- c("RxHist", "NovelHist","HIV", "SmearStatus", "DSTs", drugs)

monthsmodeled <- c(1:6, 9, 12, 18)

eventtypes <- list(TBonset=1, TBdiagnosis=2, treatmentstart=3, treatmentend=4, death=5, pretreatmentltfu=6, carryforward=7)
# diagnosis will incorporate DST and regimen selection), treatment start will incorporate ADR, 
# For each event, track RxHist, NovelHist, current susceptibilities, TB state.
statetypes <- list(undiagnosed=1, diagnosed=2, treating=3, treating_adr=4, failed=5, pendingrelapse=6, cured=7, deceased=8, relapsed=9, pretreatmentlost=10)
regimentypes <- as.list(1:length(regimens)); names(regimentypes) <- regimens; regimentypes$none <- 0


# Set up cohort of patients by probability (can choose N later)
cohort.probs <- function(params, patientvars)
{

  if(params['MOXI-R-highlevel_in_RIF-S'] > params['MOXI-R-any_in_RIF-S']) params['MOXI-R-highlevel_in_RIF-S'] <- 0.99* params['MOXI-R-any_in_RIF-S']
  if(params['MOXI-R-highlevel_in_RIF-R'] > params['MOXI-R-any_in_RIF-R']) params['MOXI-R-highlevel_in_RIF-R'] <- 0.99* params['MOXI-R-any_in_RIF-R']
    
  types <- expand.grid(c(0,1), #RxHist
                       c(0), # Novelhist
                       c(0,1), # HIV
                       c(0,1), # smearstatus
                       c(0), #DSTs
                       c(0,1), # RIF
                       c(0,1), # INH
                       c(0,1), # PZA
                       c(0,1), #MOXI
                       c(0,1), #Partialmoxi
                       c(0), #B
                       c(0)) # Pa
  colnames(types) <- patientvars
  types <- data.frame(types)
  
  probs <- numeric(nrow(types))
  
  probs <- ifelse(types$RxHist, params['Retreatment_in_All'], 1-params['Retreatment_in_All'])
  probs <- probs* ifelse(types$HIV, params['HIV_in_all'], 1-params['HIV_in_all'])
  probs[types$HIV==0] <- probs[types$HIV==0]*ifelse(types$SmearStatus[types$HIV==0], params['Smearpos_in_HIVneg'], 1-params['Smearpos_in_HIVneg'])
  probs[types$HIV==1] <- probs[types$HIV==1]*ifelse(types$SmearStatus[types$HIV==1], params['Smearpos_in_HIVpos'], 1-params['Smearpos_in_HIVpos'])
  
  probs[types$RxHist==0] <- probs[types$RxHist==0]*ifelse(types$RIF[types$RxHist==0], params['RIF-R_in_New'], 1-params['RIF-R_in_New'])
  probs[types$RxHist==1] <- probs[types$RxHist==1]*ifelse(types$RIF[types$RxHist==1], params['RIF-R_in_Retreatment'], 1-params['RIF-R_in_Retreatment'])

  probs[types$RxHist==0&types$RIF==0] <- probs[types$RxHist==0&types$RIF==0]*ifelse(types$INH[types$RxHist==0&types$RIF==0], params['INH-R_in_NewRIF-S'], 1-params['INH-R_in_NewRIF-S'])
  probs[types$RxHist==1&types$RIF==0] <- probs[types$RxHist==1&types$RIF==0]*ifelse(types$INH[types$RxHist==1&types$RIF==0], params['INH-R_in_RetreatmentRIF-S'], 1-params['INH-R_in_RetreatmentRIF-S'])
  probs[types$RIF==1] <- probs[types$RIF==1]*ifelse(types$INH[types$RIF==1], params['INH-R_in_RIF-R'], 1-params['INH-R_in_RIF-R'])
  
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
  
  zm1 <- imputecombo(OR=params["OR-PZA-if_MOXI"], z=params['PZA-R_in_RIF-R'], m=params['MOXI-R-any_in_RIF-R'])
  probs[types$RIF==1] <- probs[types$RIF==1]*
    ((zm1[1]*types$PZA*types$MOXI + zm1[2]*types$PZA*(1-types$MOXI) + zm1[3]*(1-types$PZA)*types$MOXI + zm1[4]*(1-types$PZA)*(1-types$MOXI))[types$RIF==1])
                                                                                  
  zm2 <- imputecombo(OR=params["OR-PZA-if_MOXI"], z=params['PZA-R_in_RIF-S'], m=params['MOXI-R-any_in_RIF-S'])
  probs[types$RIF==0] <- probs[types$RIF==0]*
    ((zm2[1]*types$PZA*types$MOXI + zm2[2]*types$PZA*(1-types$MOXI) + zm2[3]*(1-types$PZA)*types$MOXI + zm2[4]*(1-types$PZA)*(1-types$MOXI))[types$RIF==0])
  
  probs[types$MOXI==1&types$RIF==1] <- probs[types$MOXI==1&types$RIF==1]*ifelse(types$partialmoxi[types$MOXI==1&types$RIF==1],  1 - params['MOXI-R-highlevel_in_RIF-R']/params['MOXI-R-any_in_RIF-R'], 
                                                                    params['MOXI-R-highlevel_in_RIF-R']/params['MOXI-R-any_in_RIF-R'])

  probs[types$MOXI==1&types$RIF==0] <- probs[types$MOXI==1&types$RIF==0]*ifelse(types$partialmoxi[types$MOXI==1&types$RIF==0], 1 - params['MOXI-R-highlevel_in_RIF-S']/params['MOXI-R-any_in_RIF-S'], 
                                                                                params['MOXI-R-highlevel_in_RIF-S']/params['MOXI-R-any_in_RIF-S'])
  
  probs[types$partialmoxi==1&types$MOXI==0] <- 0
  # # here, MOXI indicates any moxi resistance, and partialmoxi indicates only low-level/partial resistance (and partial activity)

  typesandprobs <- cbind(types, probs)
  colnames(typesandprobs) <- c(colnames(types), "Freq")
  
  typesandprobs <- typesandprobs[do.call(order, as.data.frame(typesandprobs[,patientvars])),] 
  
  # Split the most common types into 10, 10x smaller rows, so they'll be run more often in simulations
  typesandprobs <- rbind(typesandprobs, do.call(rbind, replicate(subset(typesandprobs, RxHist==0&NovelHist==0&RIF==0&INH==0&PZA==0&MOXI==0&partialmoxi==0), n = 9, simplify = F)))
  typesandprobs[typesandprobs$RxHist==0&typesandprobs$NovelHist==0&typesandprobs$RIF==0&typesandprobs$INH==0&typesandprobs$PZA==0&typesandprobs$MOXI==0&typesandprobs$partialmoxi==0, "Freq"] <- 
    typesandprobs[typesandprobs$RxHist==0&typesandprobs$NovelHist==0&typesandprobs$RIF==0&typesandprobs$INH==0&typesandprobs$PZA==0&typesandprobs$MOXI==0&typesandprobs$partialmoxi==0, "Freq"]/10
  
  return(typesandprobs)
}



set.scenario <- function(scenario="0")
{
  # No new regimen (0)
  # with less xpert (0o)
  # with max xpert (0x)
  
  # New regimen considered for known RR only (1)
  #   With lowered Xpert (1o)
  #   With max Xpert (1x)
  #   
  
  # New regimen considered for all TB diagnoses
  #   Without up-front DST (2)
  #     Using shorter duration for everyone (2a)
  #     Using longer duration for everyone (2b)
  # other DST approaches, but with 4 or 6 months for all without known FQ_R: add a (4mo) or b (6mo)
  
  #   With up-front RIF DST only (3)
  #     Would use longer duration only if RIF-R detected
  #   With lowered Xpert (3o)
  #   With max Xpert (3x)
  
  #   With up-front step-wise RIF DST, then FQ DST if RIF-R (4)
  #     Would use longer duration if RIF-R or if no DST, FQ-S; alternative regimen if RIF-R, FQ-R
  #   With up-front DST for RIF and FQ for everyone (5)
  #     Would use longer duration if RIF-R or if no DST, FQ-S; alternative regimen if FQ-R (RIF-R or RIF-S) 
  ## need to define regimen_s, regimen_r, regimen_m, xpertxdr
  intervention <- list()
  intervention <- within(intervention,
                         {
                           # treatment focused analyses
                           if (scenario=="0") { regimen_s <- "hrze"; regimen_r <- "mdr"; xpertxdr <- "none"; allxpert<-"current" }
                           if (scenario=="0x") { regimen_s <- "hrze"; regimen_r <- "mdr"; xpertxdr <- "none"; allxpert<-"max" }
                           if (scenario=="0o") { regimen_s <- "hrze"; regimen_r <- "mdr"; xpertxdr <- "none"; allxpert<-"low" }
                           if (scenario=="1") { regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"current" }  
                           if (scenario=="1x") {regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"max"}
                           if (scenario=="1o") {regimen_s <- "hrze"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"low"}
                           if (scenario=="3") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"current" }  
                           if (scenario=="3x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "none" ; allxpert<-"max" }
                           if (scenario=="3xa") { regimen_s <- "bpamz4"; regimen_r <- "bpamz4r"; xpertxdr <- "none" ; allxpert<-"max" }
                           if (scenario=="3xb") { regimen_s <- "bpamz6"; regimen_r <- "bpamz6"; xpertxdr <- "none" ; allxpert<-"max" }
                           if (scenario=="3o") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"low" }  
                            #possibly also including no Xpert at all approaches:
                           if (scenario=="2a") { regimen_s <- "bpamz4"; regimen_r <- "bpamz4r"; xpertxdr <- "none"; allxpert<-"current"  }  
                           if (scenario=="2b") { regimen_s <- "bpamz6"; regimen_r <- "bpamz6"; xpertxdr <- "none"; allxpert<-"current"  }  

                           # DST focused analyses, universal regimens:
                           if (scenario=="4x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "stepwise"; allxpert<-"max"  }  
                           if (scenario=="5x") { regimen_s <- "bpamz4"; regimen_r <- "bpamz6"; xpertxdr <- "simultaneous"; allxpert<-"max"  }  
                           if (scenario=="4xa") { regimen_s <- "bpamz4"; regimen_r <- "bpamz4r"; xpertxdr <- "stepwise"; allxpert<-"max"  }  
                           if (scenario=="5xa") { regimen_s <- "bpamz4"; regimen_r <- "bpamz4r"; xpertxdr <- "simultaneous"; allxpert<-"max"  }  
                           if (scenario=="4xb") { regimen_s <- "bpamz6"; regimen_r <- "bpamz6"; xpertxdr <- "stepwise"; allxpert<-"max"  }  
                           if (scenario=="5xb") { regimen_s <- "bpamz6"; regimen_r <- "bpamz6"; xpertxdr <- "simultaneous"; allxpert<-"max"  }  
                           
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
  if(intervention$allxpert=="max") params["Xpert_current_new"] <- params["Xpert_current_rerx"] <- params["Xpert_max"]
  if(intervention$allxpert=="low") {params["Xpert_current_new"]  <- params["Xpert_low"]; params["Xpert_current_rerx"] <- params["Xpert_current_rerx"]/2 }
  
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


allDSTresults <- function(scenario, Result, patient, params, stochasticmode=TRUE) 
  #asuming they are considered for novel intervention, and assuming Xpert is done, but now (** changed 20190414 **) factoring TB detection (Result[1,]) into allDSTresults.
{
  if (missing(scenario)) scenario <- "0"
  intervention <- set.scenario(scenario)
  
  N <- nrow(patient)
  
  FullDST <- array(NA, dim=c(N, 4)); colnames(FullDST) = c("none", "RR", "FQR", "both")
  # make it 2 if it's high-level FQ-R, 1 otherwise (for FQR and for both, and for fqdetect)  
  
  fqdetect <- Result[,1]*((patient$partialmoxi==1)*prob(params["Sens_FQ_low"], N, stochasticmode) + 2*(patient$MOXI==1&patient$partialmoxi==0)*prob(params["Sens_FQ_high"], N, stochasticmode))
  # RR and FQR
  FullDST[,4] <- Result[,2]*(intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")*fqdetect
  ## RR only
  FullDST[,2] <- Result[,2]*( (intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")*as.numeric(fqdetect==0) + 
                                (1-(intervention$xpertxdr=="simultaneous"|intervention$xpertxdr=="stepwise")) )
  # FQ only
  FullDST[,3] <- (1-Result[,2])*(intervention$xpertxdr=="simultaneous")*fqdetect
  
# neither 
  FullDST[,1] <- as.numeric(rowSums(FullDST[,2:4])==0)
  
  return(FullDST)  
}


Regimenselect <- function(scenario, Result, Considernovel, FullDST, patient, stochasticmode) 
  # take DST results + scenario + patient history, and choose regimen. FullDST now includes TB detection or no.
  ##** need to modify array handling to allow just a single patient
{
  if (missing(scenario)) scenario <- "0"
  intervention <- set.scenario(scenario)
  
  regs <- array(0, dim=c(nrow(patient), length(regimens))); colnames(regs) <- regimens

  # Consider novel
  # Detect FQ (with or without RR)
  regs[,'backup'] <- regs[,'backup'] + Considernovel*as.numeric(FullDST[,"both"]>0)
  regs[,'backup'] <- regs[,'backup'] + Considernovel*as.numeric(FullDST[,"FQR"]>0)
  
  # # Detect FQ and RR
  # regs[,'bpal'] <- regs[,'bpal'] + Considernovel*as.numeric(FullDST[,"both"]==2)*(1-patient$NovelHist)
  # regs[,'bpammz6'] <- regs[,'bpammz6'] + Considernovel*as.numeric(FullDST[,"both"]==1)*(1-patient$NovelHist)
  # regs[,'backup_mdr'] <- regs[,'backup_mdr'] + Considernovel*as.numeric(FullDST[,"both"]>0)*patient$NovelHist
  # # Detect FQ only
  # regs[,'bpal'] <- regs[,'bpal'] + Considernovel*as.numeric(FullDST[,"FQR"]==2)*(1-patient$NovelHist) 
  # regs[,'bpammz6'] <- regs[,'bpammz6'] + Considernovel*as.numeric(FullDST[,"FQR"]==1)*(1-patient$NovelHist) 
  # regs[,'backup'] <- regs[,'backup'] + Considernovel*as.numeric(FullDST[,"FQR"]>0)*(patient$NovelHist)
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
  active[,"bpamz4"] <- active[,"bpamz4r"] <- active[,"bpamz6"] <- paste0(ifelse(active[,"BDQ"]==1, "", "B"),ifelse(active[,"PA"]==1, "", "Pa"), 
                                                   ifelse(active[,"MOXI"]==1, ifelse(active[,"partialmoxi"]==1, "m", ""), "M"), 
                                                   ifelse(active[,"PZA"]==1, "", "Z") )
  active[,"backup"] <- paste0(ifelse(active[,"BDQ"]==1, "", "B"), ifelse(active[,"PA"]==1, "", "Pa"), 
                              "M", 
                              ifelse(active[,"PZA"]==1, "", "Z") )

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

make.recurrence.matrix <- function(params, bpamzhigh = FALSE)
{
  # this just gives the fraction p that recur, of those (1-f, where f is failures) who are completing the specified amount of therapy without failure.
  # assume there are also f who fail, and that f = failures_per_relapse*(p*(1-f)) => f = p*(failures_pre_repalse)/(1 + p*failures_per_relapse)
  # for a total of f + p*(1-f) = f + f/failures_per_Relapse = p*(reailures_per_relapse + 1)/(1 + p*failures_per_relapse) --> gets used later in ADR, and total adr will get divided the same way. 
  # when our data combine failures and relapses, we estimate just the relapse fraction here
  recurrencematrix <- array(NA, dim=c(length(activeregimens),length(monthsmodeled))); dimnames(recurrencematrix) <- list(activeregimens,monthsmodeled)

  recurrencematrix["HR(ZE)",] <- wallis(cxconv = reversewallis(params["HRZE_pooled_relapse"]/(params["INH_multiplier"]*params["pooled_INH_fraction"] + (1-params["pooled_INH_fraction"])),6, params=params), months = monthsmodeled, params=params)
  # HR(ZE)_relapse * (INH_multiplier*pooled_INH_fraction + (1-pooled_INH_fraction)) = HRZE_pooled_relapse
  recurrencematrix["R(ZE)",] <- params["INH_multiplier"]*recurrencematrix["HR(ZE)",] 
  if(bpamzhigh) recurrencematrix["BPaMZ",] <- wallis(reversewallis(params["HRZE_pooled_relapse"],5,params=params), monthsmodeled, params=params) else
    recurrencematrix["BPaMZ",] <- wallis(reversewallis(params["HRZE_pooled_relapse"],4,params=params), monthsmodeled, params=params)
  recurrencematrix["BPaM",] <- wallis(params["BPaM_cxconv"], monthsmodeled,params=params)
  recurrencematrix["BPaZ",] <- wallis(params["BPaZ_cxconv"], monthsmodeled,params=params)
  recurrencematrix["BMZ",] <- recurrencematrix["PaMZ",] <- recurrencematrix["BPaZ",]
  if(round(params["MDR_failrelapse_FQ-S"],3)==.070){ # for sensitivity analysis with better mdr soc, a quick hack solution that won't handle changes in parameters well
    recurrencematrix["MDR, FQ-S",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-S"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/2, params = params)
    recurrencematrix["MDR, FQ-R",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-R"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/2, params = params)
  } else
    {
      recurrencematrix["MDR, FQ-S",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-S"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/3, params = params)
      recurrencematrix["MDR, FQ-R",] <- wallis(reversewallis(params["MDR_failrelapse_FQ-R"]*1/(1+params["Failures_per_recurrence"]),6,params=params), months = monthsmodeled/3, params = params)
      }
  recurrencematrix["H(ZE)",] <- wallis(reversewallis(params["INHmono_relapse"],6,params=params), months = monthsmodeled/1.5, params = params)
  recurrencematrix["BPa",] <- recurrencematrix["BZ",] <- recurrencematrix["PaZ",] <- recurrencematrix["BM",] <- recurrencematrix["PaM",] <- recurrencematrix["MZ",] <- recurrencematrix["R(ZE)",] 
  
  recurrencematrix["(ZE)",] <- recurrencematrix["B",] <- recurrencematrix["Pa",] <- recurrencematrix["M",] <- recurrencematrix["Z",] <- 
    c(wallis(cxconv = reversewallis(recurrence = params["Highrecurrence"],months = 6,params = params),months = monthsmodeled,params = params)[1:6], rep(params["Highrecurrence"], length(monthsmodeled)-6))
  
  # recurrencematrix["BPaL",] <- recurrencematrix["HR(ZE)",]
  # recurrencematrix["BL",] <- recurrencematrix["PaL",] <- recurrencematrix["BPa",]
  recurrencematrix["m",] <- recurrencematrix["none",] <- recurrencematrix["B",]
  # recurrencematrix["L",] <- recurrencematrix["B",]
  
  recurrencematrix[recurrencematrix>1] <- 1
  
  recurrencematrix["BPamZ",] <- increaseodds(recurrencematrix["BPaZ",], params["partialmoxiOR"])
  recurrencematrix["MDR, FQ-low",] <- increaseodds(recurrencematrix["MDR, FQ-R",], params["partialmoxiOR"])
  recurrencematrix["BPam",] <- increaseodds(recurrencematrix["BPa",], params["partialmoxiOR"])
  recurrencematrix["BmZ",] <- increaseodds(recurrencematrix["BZ",], params["partialmoxiOR"])
  recurrencematrix["PamZ",] <- increaseodds(recurrencematrix["PaZ",], params["partialmoxiOR"])
  recurrencematrix["Bm",] <- increaseodds(recurrencematrix["B",], params["partialmoxiOR"])
  recurrencematrix["Pam",] <- increaseodds(recurrencematrix["Pa",], params["partialmoxiOR"])
  recurrencematrix["mZ",] <- increaseodds(recurrencematrix["Z",], params["partialmoxiOR"])

    recurrencematrix[recurrencematrix>params["Highrecurrence"] | recurrencematrix=="NaN"] <- params["Highrecurrence"]
  
  return(recurrencematrix)
} 

make.adr.matrix <- function(params)
{
  adrmatrix <- array(NA, dim=c(length(activeregimens), 5)); dimnames(adrmatrix) <- list(activeregimens, paste0(adrs,"-R")) 
  adrmatrix[,"RIF-R"] <- c(params["adr_r"]*c(1, params["adrfactor_other"]), rep(0, length(activeregimens)-2))
  adrmatrix["MDR, FQ-S","MOXI-R"] <- adrmatrix["MDR, FQ-low","MOXI-R"] <- params["adr_mdr"]
  adrmatrix["BPaMZ","MOXI-R"] <- adrmatrix["BPamZ","MOXI-R"] <- params["adr_bpamz"]
  adrmatrix["BPaM","MOXI-R"] <- adrmatrix["BPam","MOXI-R"] <- params["adr_bpamz"]*params["adrfactor_z"]
  adrmatrix["BMZ","MOXI-R"] <- adrmatrix["PaMZ","MOXI-R"] <- adrmatrix["BmZ","MOXI-R"] <- adrmatrix["PamZ","MOXI-R"] <- params["adr_bpamz"]*params["adrfactor_other"]
  adrmatrix["BM","MOXI-R"] <- adrmatrix["PaM","MOXI-R"] <- adrmatrix["Bm","MOXI-R"] <- adrmatrix["Pam","MOXI-R"] <- params["adr_bpamz"]*params["adrfactor_other"]*params["adrfactor_z"]
  adrmatrix["MZ","MOXI-R"] <- adrmatrix["mZ","MOXI-R"] <- params["adr_bpamz"]*params["adrfactor_twodrugs"]
  adrmatrix["M","MOXI-R"] <- adrmatrix["m","MOXI-R"] <- params["Highrecurrence"]*(1+params[ "Failures_per_recurrence"])/(1+params["Highrecurrence"]*params[ "Failures_per_recurrence"])
  adrmatrix[c("BPaMZ", "BPaM", "BPaZ", "BMZ", "BPa", "BZ", "BM"),"BDQ-R"] <- 
    params["adr_bpamz"]*c(1, params["adrfactor_z"], params["adrfactor_other"], params["adrfactor_other"], params["adrfactor_other"]*params["adrfactor_z"], 
                          params["adrfactor_twodrugs"], params["adrfactor_other"]*params["adrfactor_z"])
  adrmatrix[c("BPaMZ", "BPaM", "BPaZ", "PaMZ", "BPa", "PaZ", "PaM"),"PA-R"] <- 
    params["adr_bpamz"]*c(1, params["adrfactor_z"], params["adrfactor_other"], params["adrfactor_other"], params["adrfactor_other"]*params["adrfactor_z"], 
                          params["adrfactor_twodrugs"], params["adrfactor_other"]*params["adrfactor_z"])
  adrmatrix["HR(ZE)","INH-R"] <- params["adr_r"] # added assumption
  adrmatrix["H(ZE)","INH-R"] <- params["adr_r"]*params["adrfactor_twodrugs"] # this shouldn't spell universal failure?
  adrmatrix["B", "BDQ-R"] <- adrmatrix["Pa", "PA-R"] <- params["Highrecurrence"]*(1+params[ "Failures_per_recurrence"])/(1+params["Highrecurrence"]*params[ "Failures_per_recurrence"])
  adrmatrix[c("BPamZ", "BPam", "BmZ", "Bm"),"BDQ-R"] <- 
    params["adrfactor_partialmoxi"]*adrmatrix[c("BPaZ", "BPa", "BZ", "B"),"BDQ-R"] 
  adrmatrix[c("BPamZ", "BPam", "PamZ", "Pam"),"PA-R"] <- 
    params["adrfactor_partialmoxi"]*adrmatrix[c("BPaZ", "BPa", "PaZ", "Pa"),"PA-R"] 

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
  needdiagnosis <- ((last[,"TBstate"]==statetypes$undiagnosed)|(last[,"TBstate"]==statetypes$failed)| (last[,"TBstate"]==statetypes$relapsed)| (last[,"TBstate"]==statetypes$pretreatmentlost))
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
  # document DSTs performed. Add 1 if xpert done, and 0.1 if xdr done (i.e. Xpert&ifelse(intervention$xpertxdr=="simultaneous",1,Result[2])).
  # note: we're documenting all the DSTs that would be done if treated, even for those not getting treatment. I'll need to interpret these by combining with either diagnosed (if iltfu gets DST) or treated or Currentregimen (if dst means starting regimen)
  Considernovel <- Novelavailable(currentcohort, params, stochasticmode)
  Xpert <- Xpertdone(scenario, currentcohort, params, stochasticmode) # is xpert done? depends on scenario and patient
  Result <- Xpertresult(Xpert, currentcohort, params, stochasticmode) # Xpert MTB/RIF result. Dependent on Xpertdone (And therefore indirectly on scenario)
  FullDST <- allDSTresults(scenario, Result, currentcohort, params, stochasticmode) # Full results of DST, depending on Xpert MTB/RIF result (indirectly on scenario and patient) and on Xpert XDR use scenario
  Regimen <- Regimenselect(scenario, Result, Considernovel, FullDST, currentcohort, stochasticmode) # regimen selected for each patient, depending on the DST results, the novel availability, the patient's prior novel regimen exposure, and the scenario's preferred regimens
  
  if(!stochasticmode) stop("Can't use the matrix regimen output in the lines that follow, still need to code this version")

  current[,"DSTs"] <- last[,"DSTs"] + Xpert + 0.1*Xpert*(as.numeric(intervention$xpertxdr=="simultaneous") + Result[,2]*as.numeric(intervention$xpertxdr=="stepwise"))
    
  # event will now be treatment start, but we'll wait to add that only for those who start treatment
  # determine time of this event
  # DST, and selection of a nonstandard regimen, may incur delays or pretreatment losses to follow up
  needtreatment <- (last[,"TBstate"] == statetypes$diagnosed)
  
  starttime <- Xpert*params["DSTdelay"] 
  # delays for non-first-line regimens:
  if (intervention$regimen_s=="hrze") starttime <- starttime + (Regimen!=intervention$regimen_s)*params["Regimendelay"]
  if (intervention$regimen_s %in% c("bpamz4", "bpamz6")) starttime <- starttime + (!(Regimen %in% c("bpamz4", "bpamz4r", "bpamz6")))*params["Regimendelay"]
  morttime <- rexp(N, params["Monthlymortality_TB"] + last[,"HIV"]*params["Monthlymortality_HIV"])
  
  deaths <-  needtreatment&(morttime < starttime)
  current[deaths, "eventtime"] <- last[deaths, "eventtime"] + morttime[deaths]; 
  current[deaths, "eventtype"] <- eventtypes$death
  current[deaths, "TBstate"] <- statetypes$deceased; # leave others as assigned above (treated or lost)
  
  lost <- needtreatment&(!deaths)&(prob(params["Unavoidableloss"] + Xpert*params["DSTloss"] + (Regimen!=intervention$regimen_s)*params["Regimenloss"], N, stochasticmode))
  current[lost,"eventtype"] <- eventtypes$pretreatmentltfu
  current[lost,"TBstate"] <- statetypes$pretreatmentlost
  
  current[needtreatment&(!deaths), "eventtime"] <- last[needtreatment&(!deaths), "eventtime"] + starttime[needtreatment&(!deaths)]; #assumes pretreatment losses happen at starttime
  
  treated <- needtreatment&(!deaths)&(!lost)
  current[treated,"eventtype"] <- eventtypes$treatmentstart
  
  #those who aren't lost start treatment, possibly with ADR
  # assign regimen
  current[treated, "Currentregimen"] <- as.numeric(regimentypes[Regimen])[treated]
  active <- activeregimen(make.active.regimen.matrix(), currentcohort, Regimen)
  active[active==""] <- "none"
  ## for each patient, sample risks of ADR based on active regimen, and if any occur, change phenotype in cohort and set state as _adr
  adrrisks <- make.adr.matrix(params)[active,]
  adr <- array(rbinom(N*length(adrs), 1, adrrisks[]), dim=c(N,length(adrs))); adr[is.na(adr)] <- 0 # this samples resistance risk for each drug
  # update state accordingly (haven't removed deaths yet)
  current[treated, "TBstate"] <- ((rowSums(adr)>=1)*statetypes$treating_adr + (rowSums(adr)<1)*statetypes$treating)[treated]
  # update resistance accordingly, if ADR: change 0s (for both FQ columns) in current to 1 if adr is 1 (and !lost) 
  current[treated, adrs] <- pmax(last[treated, adrs], adr[treated,])
  # if moxi adr but already were moxi, then partialmoxi now needs to change to zero
  # (but first make sure we're allowing moxi adr for baseline partialmoxi -- yes, because the "m"s have a risk of FQ-R ADR in adrmatrix.)
  current[(treated&(adr[,2]==1)&(last[,"partialmoxi"]==1)), "partialmoxi"] <- 0
  
    
  current[!needtreatment,"eventtype"] <- eventtypes$carryforward
  
  # save before moving to next event
  
  return(current)
}


## THIRD EVENT ###
## treatment end (due to end of regimen or loss to follow up)
## will carry forward those who aren't dead or getting treated

treatmentend <- function(last, params, N, eventtypes, statetypes, regimentypes, bpamzhigh=FALSE)
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
  current[failrelapseadr & (last[,"Currentregimen"] %in% which(regimens %in% c("bpamz4", "bpamz4r","bpamz6"))), "NovelHist"] <- 1
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
  adrrisks <- make.adr.matrix(params)
  adrrisks[is.na(adrrisks)] <- 0
  totaladrrisks <- adrrisks[,1] + (1-adrrisks[,1])*(adrrisks[,2] + (1-adrrisks[,2])*(adrrisks[,3] + (1-adrrisks[,3])*(adrrisks[,4] + (1-adrrisks[,4])*adrrisks[,5])))
  # if recurrences are R of those at risk, and failures are p*R (p=failures_per_Rec), 
  # then overall failures or recurrences will be f + (1-f)*R = (p+1)*(1-f)*R --> f=pR/(1+pR) -->f+(1-f)R=R(1+p)/(1+pR)
  # so then really, the cap on poor outcomes is a little over highrecurrence, at 83% rather than 80, but I think that's okay.
  # and it will be higher still when ADR=1, so I'll reduce the max ADR to highrecurrence (p(1+R)/(1+pR) to allow for self-cure there as well. 
  failures_or_relapses_remaining <- apply(make.recurrence.matrix(params, bpamzhigh)*(1+params["Failures_per_recurrence"])/(1+make.recurrence.matrix(params, bpamzhigh)*params["Failures_per_recurrence"]), 2, "-", totaladrrisks)
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
  current[needoutcome & (last[,"Currentregimen"] %in% which(regimens %in% c("bpamz4", "bpamz4r", "bpamz6"))), "NovelHist"] <- 1
  current[needoutcome,"Currentregimen"] <- 0
  
  current[!onRx,"eventtype"] <- eventtypes$carryforward
  
  return(current)
}

## EVENT 4 = RELAPSES (vs other mortality while pending relapse), and carry others forward if they're not dead ## 
relapseevent <-  function(last, params, N, eventtypes, statetypes, regimentypes)
{
  current <- last
  
  relapsers <- last[,"TBstate"]==statetypes$pendingrelapse
  
  relapsetime <- rexp(N, 1/params["Recurrence_interval"]) 
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

# add mortality after cure
#(already included above for those pending relapse, and for those with TB, so will just run this as the final step of the model)
addnaturalmortality <- function(last, params, N, eventtypes, statetypes, maxtime=100*12)
{
  current <- last
  
  cures <- last[,"TBstate"]==statetypes$cured
  
  morttime <- rexp(N, params["Monthlymortality_background"] + last[,"HIV"]*params["Monthlymortality_HIV"]) 
  morttime[morttime>maxtime] <- maxtime
  deaths <- cures
    
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

modelcourse <- function(scenario="0", cohort, params, reps=1, steplimit=16, stochasticmode=TRUE, bpamzhigh=FALSE) # will need to repeat same for each rep of a given patient, if probabilistic events
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
  course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, bpamzhigh), along=3)
  course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)

  ## and back to event 5= DIAGNOSIS AGAIN ## 
  course <- abind(course, diagnosisevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  course <- abind(course, treatmentinitiationattempt(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, scenario), along=3)
  course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, bpamzhigh), along=3)
  course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  
  #continue to loop up to sooner of steplimit/4 loops or no more (uncured or dead) or those still alive and uncured have run for 20 years:
  while(dim(course)[3]<steplimit & 
        mean(course[,"eventtype",dim(course)[[3]]] %in% c(statetypes$deceased, statetypes$cured))<1 & 
        mean(course[,"eventtime",dim(course)[[3]]][!course[,"eventtype",dim(course)[[3]]] %in% c(statetypes$deceased, statetypes$cured)] > 20*12)<1)
  {
    course <- abind(course, diagnosisevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
    course <- abind(course, treatmentinitiationattempt(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, scenario), along=3)
    course <- abind(course, treatmentend(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes, bpamzhigh), along=3)
    course <- abind(course, relapseevent(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes, regimentypes), along=3)
  }
    
  course <- abind(course, addnaturalmortality(course[,,dim(course)[[3]]], params, N, eventtypes, statetypes), along=3)
  
  newcourse <- array(course, dim=c(dim(course)[1]/reps,reps,dim(course)[2:3]))
  newcourse <- aperm(newcourse, c(1,3,4,2))
  dimnames(newcourse) <- list("patient"=c(), "characteristic"=trackeditems,"timestep"=c(),"rep"=c())
  
  return(newcourse)
}


# set the model type
stochasticmode <- T
