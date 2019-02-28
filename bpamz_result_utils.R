require(reshape2); require(ggplot2)

require(parallel)
no_cores <- detectCores()


bootsample <- function(c, include, desiredsize, course_outcomes) 
{
  bootindex <- sample((1:nrow(c))[include], desiredsize, replace=T, prob=c$Freq[include]/sum(c$Freq[include]))
  randindex <-  sample(1:dim(course_outcomes)[[3]], size = desiredsize, replace=T)
  return(sapply(X = 1:dim(course_outcomes)[[1]], function(x) course_outcomes[cbind(x, bootindex, randindex)]))
}

#...which then feeds into
outcomeboot <- function(individualoutcomefunction, course, c, include=rep(TRUE, dim(course)[[1]]), copies=50, desiredsize=1e5)
{
  course_outcomes <- apply(X=course, FUN=individualoutcomefunction, MARGIN = c(1,4))
  return(replicate( copies, bootsample(c, include, desiredsize, course_outcomes)))
}

#... or for a list of courses for different scenarios
loutcomeboot <- function(individualoutcomefunction, simoutput, c, 
                         include=rep(TRUE, nrow(c)), copies=50, desiredsize=1e5, course_outcomes=NA, 
                         oncluster=F, ...) 
{
  freqweights <- c$Freq # simoutput[[1]][,"Freq",1,1]
  if(oncluster) 
  {
    cores <- detectCores()
    clust <- makeCluster(cores)
    clusterExport(clust, c("course_outcomes", "individualoutcomefunction", "copies", "bootsample", "c", "include", "desiredsize"), envir=environment())
    l <- parLapply(clust, simoutput, function(y) {  
      if(is.na(course_outcomes)) course_outcomes <- apply(X=y, FUN=individualoutcomefunction, MARGIN = c(1,4), ...)
      return(replicate( copies, bootsample(c, include, desiredsize, course_outcomes)))
    } )
  }
  else 
  {
    l <- lapply(simoutput, function(y) {  
      if(is.na(course_outcomes)) course_outcomes <- apply(X=y, FUN=individualoutcomefunction, MARGIN = c(1,4), ...)
      return(replicate( copies, bootsample(c, include, desiredsize, course_outcomes)))
    } )
  }
  return(l)
}

step4status <- function(patiententry) {patiententry[c("TBstate", "eventtime"),4] } 
step8status <- function(patiententry) {patiententry[c("TBstate", "eventtime"),8] }
step4812cure <- function(patiententry) {c(patiententry[c("TBstate"),c(4, 8, 12)]=="7", patiententry[patientvars,1]) } 
step4cure <- function(patiententry) {c(patiententry[c("TBstate"),4]=="7", patiententry[patientvars,1]) } 

# redefining from cohort file, for single-patient course:
# but note that need to keep characteristics here, and call this function with an include of all pts at risk for the characteristic, not just those who have it at baseline (e.g. when talking about adr)
# will return time in states of interest, final tbstate, and initial characteristics
time.in.state <- function(patiententry, states=8, characteristics=c(), characteristicvals=1, 
                          cutofftime=100*12, carryforward=F)
{
  elapsedtimes <- patiententry["eventtime",2:ncol(patiententry)] - patiententry["eventtime",1:((ncol(patiententry))-1)] 
  elapsedtimes[patiententry["eventtime",2:ncol(patiententry)]>cutofftime] <- cutofftime - (patiententry["eventtime",1:(ncol(patiententry)-1)])[patiententry["eventtime",2:(ncol(patiententry))]>cutofftime]
  elapsedtimes[elapsedtimes<0] <- 0
  indices <-  patiententry["TBstate",1:(ncol(patiententry)-1)] %in% states
  if(length(characteristics)==1) indices <- (indices & patiententry[characteristics,1:(ncol(patiententry)-1)]==characteristicvals)
  if(length(characteristics)>1) indices <- indices & apply(patiententry[characteristics,1:(ncol(patiententry)-1)]==characteristicvals, 2, all)
  t <- sum(elapsedtimes*indices)
  if (carryforward)  # if carryforward, will maintain final state incl infectous states from end of data to cutofftime.
    if (patiententry["TBstate",ncol(patiententry)] %in% states)
      t <- t + max(0, cutofftime-patiententry["eventtime",ncol(patiententry)])
  return(c(t, patiententry["TBstate",ncol(patiententry)], patiententry[patientvars,1]))
}

still.in.state <- function(patiententry, states=8, t=36)
{
  # find last entry with time < time
  maxed <- (max(patiententry["eventtime",])>=t)
  # is it in states?
  return(c(
    patiententry["TBstate", ifelse(maxed, which.max(patiententry["eventtime",]>t)-1, length(patiententry["eventtime",]))] %in% states,
    maxed,
    patiententry[patientvars,1]))
}

time.and.courses.to.cure <- function(patiententry) # also include number of treatment courses
{  c(
  ifelse(patiententry["TBstate",ncol(patiententry)]==statetypes$cured, patiententry["eventtime",which.max(patiententry["TBstate",]==statetypes$cured)], NA),
  sum(patiententry["eventtype",]==eventtypes$treatmentstart)  
)
}



treatmentsuccess <- function(patiententry) #  treated then cured, and append rif and moxi status, and partialmoxi and pza
{ return(c(patiententry["eventtype",3]==eventtypes$treatmentstart, 
           patiententry["TBstate",4]==7,
           patiententry["RIF",1],
           patiententry["MOXI",1],
           patiententry["partialmoxi",1],
           patiententry["PZA",1]) ) }


# for sensitivty analysis:



shortmodelcourse <- function(scenario="0", cohort, params, reps=1, steplimit=4, stochasticmode=TRUE) # will need to repeat same for each rep of a given patient, if probabilistic events
{
  repcohort <- do.call(rbind, replicate(reps, cohort, simplify=FALSE))
  intervention <- set.scenario(scenario)
  N=nrow(repcohort)
  
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
  
  newcourse <- array(course, dim=c(dim(course)[1]/reps,reps,dim(course)[2:3]))
  newcourse <- aperm(newcourse, c(1,3,4,2))
  dimnames(newcourse) <- list("patient"=c(), "characteristic"=trackeditems,"timestep"=c(),"rep"=c())
  
  return(newcourse)
}


# 
# # can look at this equally-weighted "cohort":
# lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",16]==7)))) 
# lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",4]==1)))) 
# # and particular subsets
# lapply(impact, FUN = function(x) summary(apply(x[c$RIF==1,,,], 4, function(y) mean(y[,"TBstate",10]==7)))) #RR
# lapply(impact, FUN = function(x) summary(apply(x[c$MOXI==0,,,], 4, function(y) mean(y[,"TBstate",10]==7)))) 
# 
# #** Possible To do: add worse HRZE outcomes for PZA-R. Would need to differentiate the active regimens for (ZE) to +-Z(E)
# 
# 
# # example:  proportion diagnosed by step two:
# fun <- function(patiententry) { return(c(patiententry["TBstate",2]==2, NA))}
# l <- loutcomeboot(individualoutcomefunction = fun, simoutput = impact, c = c, include=include, copies=50, desiredsize=1e5)
# lapply(l, function(x) summary(apply(x[,1,], 2, mean)))
# #very consistent across copies, but slight differences between scenarios (not consistent between runs of original simulation), suggesting need for reps>1000
# 
# ####
# # plot states over time? this is by model step rather than time
# par(mfrow=c(2,2))
# lapply(impact, function(y) barplot(array(unlist(apply(y[,"TBstate",2:20,], 2, tabulate)), dim = c(8,19)), beside=F, col = rainbow(8), legend.text = names(statetypes)))
# 
