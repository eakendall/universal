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
outcomeboot <- function(individualoutcomefunction, course, c, include=rep(1, dim(course)[[1]]), copies=10, desiredsize=1e4)
{
  course_outcomes <- apply(X=course, FUN=individualoutcomefunction, MARGIN = c(1,4))
  return(replicate( copies, bootsample(c, include, desiredsize, course_outcomes)))
}
# or for a list of courses for different scenarios
loutcomeboot <- function(individualoutcomefunction, simoutput, c, 
                         include=rep(1, nrow(c)), copies=10, desiredsize=1e4, course_outcomes=NA, 
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
step4812cure <- function(patiententry) {patiententry[c("TBstate"),c(4, 8, 12)]=="7" } 

# redefining from cohort file, for single-patient course:
# but note that need to keep characteristics here, and call this function with an include of all pts at risk for the characteristic, not just those who have it at baseline (e.g. when talking about adr)
# will return time in states of interest, and final tbstate
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
  return(c(t, patiententry["TBstate",ncol(patiententry)]))
}

still.in.state <- function(patiententry, states=8, time=36)
{
  # look at last entry with time < time, is it in states?
  maxed <- (max(patiententry["eventtime",])>t)
  return(c(
    patiententry["TBstate", ifelse(maxed, which.max(patiententry["eventtime",]>t)-1, length(patiententry["eventtime",]))] %in% states,
    maxed))
}
time.and.courses.to.cure <- function(patiententry) # also include number of treatment courses
{  c(
  ifelse(patiententry["TBstate",ncol(patiententry)]==statetypes$cured, patiententry["eventtime",which.max(patiententry["TBstate",]==statetypes$cured)], NA),
  sum(patiententry["eventtype",]==eventtypes$treatmentstart)  
)
}


require(RColorBrewer)
colors <- c(palette(brewer.pal(length(statetypes), name = "Set1")),"black")
# colors <- c(palette(brewer.pal(length(statetypes)-1, name = "Set2")), "black")
# colors <- c("black", palette(brewer.pal(length(statetypes)-1, name = "Accent")))



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
