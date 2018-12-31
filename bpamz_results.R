date <- "20181229"

load(paste0("cohort.",date,".Rdata"))
impact <- readRDS(file = paste0("impact.",date,".RDS"))
impact <- lapply(impact, function(x) x[])
dst <- readRDS(file = paste0("dst.",date,".RDS"))
delays <- readRDS(file = paste0("delays.",date,".RDS"))
                                    

# can look at this equally-weighted "cohort":
lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",20]==7)))) 
lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",4]==1)))) 
# and particular subsets
lapply(impact, FUN = function(x) summary(apply(x[c$RIF==1,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) #RR
lapply(impact, FUN = function(x) summary(apply(x[c$MOXI==0,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) 

#** Possible To do: add worse HRZE outcomes for PZA-R. Would need to differentiate the active regimens for (ZE) to +-Z(E)

# But really we want to bootstrap to construct the cohort whose outcomes will be modeled: 
#  for each element of the above lists, for each ith entry in 1st dimension (which corresponds to the rows of c), 
# sample dims 2 and 3 (across dim 4) based on ith frequncy in c$Freq, and abind the result into a course (3 dims), 
# and do this copies times to get the 4th dim (variability between simulations of the same cohort).

include <- c$RIF==1 # define a sub-cohort of interest.  
copies <- 10 # how many bootstrap copies of the whole cohort to create (for looking at variation between realizations)?
desiredsize <- 1e5 # how big to make each bootstrapped copy?

# # sample  by cohort weight from the patient type number (dimension 1 in modelcourse output)
# # and for each, sample randomly from the reps of that patient type (i.e. across dimension 4)
# # make copies of these bound across 4th dimension to get multiple realizastions of the cohort:
# boot <- function(course)
#   {
#   a <- abind(lapply(1:copies, FUN=function(x) {
#   bootindices <- sample((1:nrow(c))[include], desiredsize, replace=T, prob=c$Freq[include]/sum(c$Freq[include]))
#   course[bootindices,,,sample(1:reps,1)]}), 
#   along=4)
#   return(a)
# }
# final <- (boot(impact$baseline)) # this may require lots of memory
# or, to save memory, calculate the outcome(S) of interest for each patient type (~100 types, versus >10k in a cohort)
# individualoutcomefunction takes one patient entry matrix (characteristic x timestep) from course (e.g. course[1,,,1])
# and generates a *vector* of outcomes (ideally named)
# bootsample will sample outcomes (one vector for each patient type and rep) based on the weights of the types and randomly for each rep
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
                         include=rep(1, nrow(c)), copies=10, desiredsize=1e4, course_outcomes=NA, ...) 
{
  freqweights <- c$Freq # simoutput[[1]][,"Freq",1,1]
  l <- lapply(simoutput, function(y) {  
    
    if(is.na(course_outcomes)) course_outcomes <- apply(X=y, FUN=individualoutcomefunction, MARGIN = c(1,4), ...)
    return(replicate( copies, bootsample(c, include, desiredsize, course_outcomes)))
  } )
  return(l)
}

# # example of an individualoutcomefunction:
# fun <- function(patiententry)
#   {
#     return(quantile(patiententry["eventtime",],c(.05,.5,.95)))
# }
# patiententry <- impact$baseline[1,,,1]
# fun(patiententry)
# l <- loutcomeboot(fun, impact, c, include=include, copies=10, desiredsize=1e5)
# str(l)
# #or,
# lastalive <- function(patiententry)
# {
#   return(max((data.frame(t(patiententry)) %>% filter(TBstate<8))["eventtime"]))
# }
# lastalive(patiententry)
# 
# # trying a fun with 2-dim array output:
# twodim <- function(patiententry)
# {
#   return(patiententry[c("eventtime","TBstate"),])
# }
# twodim(patiententry)
# o <- outcomeboot(twodim, course=impact$baseline, include=include, copies=10, desiredsize=1e5)  
# str(o)
# # answer: the apply in outcomeboot corced the output of individualoutcomefunction to a vector.


# example:  proportion diagnosed by step two:
fun <- function(patiententry) { return(c(patiententry["TBstate",2]==2, NA))}
l <- loutcomeboot(individualoutcomefunction = fun, simoutput = impact, c = c, include=include, copies=50, desiredsize=1e5)
lapply(l, function(x) summary(apply(x[,1,], 2, mean)))
#very consistent across copies, but slight differences between scenarios (not consistent between runs of original simulation), suggesting need for reps>1000

####
# plot states over time? this is by model step rather than time
par(mfrow=c(2,2))
lapply(impact, function(y) barplot(array(unlist(apply(y[,"TBstate",2:20,], 2, tabulate)), dim = c(8,19)), beside=F, col = rainbow(8), legend.text = names(statetypes)))


require(RColorBrewer)
colors <- palette(brewer.pal(length(statetypes), name = "Set2"))
# require(colorspace); colors <- palette(rainbow_hcl(n = length(statetypes))); colors <- choose_palette()

# snapshots by time interval: 
# plot average states after first treatment attempt, for the equal-wieghted cohort:
par(mfrow=c(1,1), oma=c(0,0,1,1))
barplot(array(unlist(lapply(impact, FUN = function(x) tabulate(x[,"TBstate",8,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Distribution of states after second round of treatment", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)

# do the same but for the "include"+reweighted subset cohort of interest:
# converting to weighted cohort:
# first get the outcome of interest, here distribution of states and times(Added this) at time step 8:
step4status <- function(patiententry) {patiententry[c("TBstate", "eventtime"),4] } 
step8status <- function(patiententry) {patiententry[c("TBstate", "eventtime"),8] }
# and get that across the cohort:
l <- loutcomeboot(function(patiententry) {patiententry[c("TBstate", "eventtime"),4] }, impact, c, include=(c$RIF==1), copies=100, desiredsize = 1e4)
barplot(array(unlist(lapply(l, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Distribution of states after first round of treatment, RR-TB cases", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)
l <- loutcomeboot(step8status, impact, c, include=(c$RIF==1), copies=100, desiredsize = 1e4)
par(mfrow=c(1,1), oma=c(0,0,1,1))
barplot(array(unlist(lapply(l, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Distribution of states after second round of treatment, RR-TB cases", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)
# benefit for novelrrx over novelrr is more than the increase in Xpert coverage, because the new pts whose xpert increases most had longest TBmortality exposure
# additional benefit of pantB here is those whose TB or RR is missed by Xpert? ~10% for TB and another 5% rif


l <- loutcomeboot(function(patiententry) {patiententry["TBstate",] }, impact, c, include=(c$RIF==0), copies=10, desiredsize = 1e4)
lapply(l, function(x) apply(x[,2,], 2, tabulate))
# I've checked code thoroughly, am confident differences here for timestep 2 are due to chance but indicate need to increase reps in initial simulatio

# can do the same for a different subset:
l <- loutcomeboot(step4status, impact, c, include=(c$MOXI==1), copies=5, desiredsize = 1e4)
#and then tabulate across the cohort and plot as before:
par(mfrow=c(1,1), oma=c(0,0,1,1))
barplot(array(unlist(lapply(l, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Outcomes among FQ-r patients:\n", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)
#for Moxi R, Rif S patients, novel regimen looks same or slightly worse

# Differences compared to baseline: 
# in proportion cured after 1st round:
par(mfrow=c(1,1), oma=c(0,0,1,1), mar=c(3,3,3,1))
step4812cure <- function(patiententry) {patiententry[c("TBstate"),c(4, 8, 12)]=="7" } 
l <- loutcomeboot(step4812cure, impact, c, include=(c$RIF==1), copies=50, desiredsize = 1e5)
boxplot(lapply(l, function(y) apply(y[,1,]- apply(l$baseline[,1,],2, mean), 2, mean)),
        main="Incremental proportion of RIF-R cohort cured after first treatment,\ncompared to baseline average", col=colors[2:6])
boxplot(lapply(l, function(y) apply(y[,1,], 2, mean)), ylim=c(0,1),
        main="Proportion of RIF-R cohort cured after first treatment", col=colors[2:6])



# showing with and without regimen delays:
# ** should look at the same population I'm considering in the main aalysis. Probably all patietns or RIF==1.
d <- lapply(loutcomeboot(step4812cure, delays, c, include=(rep(1, length=nrow(c))), copies=50, desiredsize = 1e5), function(y) apply(y[,1,], 2, mean))
l <- lapply(loutcomeboot(step4812cure, impact, c, include=(rep(1, length=nrow(c))), copies=50, desiredsize = 1e5), function(y) apply(y[,1,], 2, mean))
str(d)
ddf <- data.frame(d); ldf <- data.frame(l)
ddf$params <- "Delay for alternate regimens"
ldf$params <- "No delays"
require(dplyr)
require(data.table)
m <- bind_rows(ldf, ddf)
mm <- melt(m, id.vars="params")
colnames(mm) <- c("Parameters", "Scenario", "Cure1prob")
mm$Parameters <- factor(mm$Parameters, levels=c("No delays", "Delay for alternate regimens"))
require(ggplot2)
ggplot(data = mm, aes(x=Scenario, y=Cure1prob)) +
  geom_boxplot(aes(fill=Parameters)) +
  ggtitle("Proportion of cohort cured after first treatment, effect of regimen delays") + theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Proportion cured") + ylim(0,1)


# as a heat map for status over time (don't group by state)
par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
times <- 0:30
whichpts<-c$RIF==1 #** can change here, should change title to match
mapply(FUN = function(a, b) { y <- a[whichpts,,,1]; 
image( 
  do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
    ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))])),
  col=colors, xaxt='n', yaxt='n', main=paste0(length(whichpts)," patient types, ",b));
axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
},
impact, names(impact))
legend("bottomright", legend = rev(names(statetypes)), fill=rev(colors), cex=0.9)
# or better, could make this 100 reps of a common or important patient type:
w <- which.max(c$RIF==1 & c$MOXI==1)
mapply(FUN = function(a, b) { y <- a[w,,,1:100]; # * set number to plot here
# mapply(FUN = function(a, b) { y <- a[w,,,]; 
image( 
  do.call('rbind', lapply(times, function(t) y["TBstate",,][cbind(apply(y["eventtime",,], 2, function(x) 
    ifelse(max(x)>t, which.max(x>t) -1, length(x)) ), 1:dim(y)[[3]])])),
  col=colors, xaxt='n', yaxt='n', main=paste0(reps," duplicated patients, ",b));
axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
},
impact, names(impact))
legend("bottomright", legend = (names(statetypes)), fill=(colors), cex=0.8)
mtext(text = "months", side = 1, outer=T, line=-1)
mtext(text = "Each horizontal line represents one realization of the patient's course under the specified scenario", side = 2, outer=T, line=-1)

# for one run, sorted:
par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
times <- seq(0,30,by=1)
y <- loutcomeboot(function(x) x[c("eventtime","TBstate"),], impact, c, include=(c$RIF==1), copies=1, desiredsize = 1e5)$baseline
z <- abind(y[,seq(1,45,by=2),1], y[,seq(2,46,by=2),1], along = 3, new.names = list("patient"=1:1e5, "step"=1:23, "a"=c("time","state")))
barplot(xlab = "months elapsed", ylab="patients in cohort", border = NA, 
        t(array(do.call('rbind', lapply(times, function(t) 
          tabulate(unlist( z[, ,"state"][
            cbind(1:dim(z)[[2]], apply(z[,,"time"], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
        )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))), 
        col=colors, beside=F)
legend("topright", legend = rev(names(statetypes)), fill=rev(colors), cex=0.8)
mtext("Status of RIF-R patients over time, in baseline scenario (no novel regimen)", side=3)

# # for one run, sorted but in heatmap mode:
# msort <- function(m)
# {return(m[do.call(order, c(decreasing = FALSE, data.frame(m[,1:ncol(m)]))),])}
# 
# par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
# times <- seq(0,30,by=0.1)
# y <- loutcomeboot(function(x) x[c("eventtime","TBstate"),], impact, c, include=(c$RIF==1), copies=1, desiredsize = 1e4)$baseline
# z <- abind(y[,seq(1,45,by=2),1], y[,seq(2,46,by=2),1], along = 3, new.names = list("patient"=1:1e4, "step"=1:23, "a"=c("time","state")))
# m <- do.call('rbind', lapply(times, function(t) z[, ,"state"][cbind(1:dim(z)[[2]], apply(z[,,"time"], 1, function(x)
#   ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))]))
# image(t(msort(t(m))), col=colors, xaxt='n', yaxt='n')
# legend("bottomleft", legend = rev(names(statetypes)), fill=rev(colors), cex=0.8)
# axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
# mtext(text = "months", side = 1, line=2)
# mtext(text = "cohort", side = 2, line=2)

# for multiple scenarios, bar plots, :
par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,2,1))
times <- seq(0,30,by=1)
yl <- loutcomeboot(function(x) x[c("eventtime","TBstate"),], impact, c, include=(c$RIF==1), copies=1, desiredsize = 1e4)
zl <- lapply(yl, function(y) abind(y[,seq(1,45,by=2),1], y[,seq(2,46,by=2),1], along = 3, new.names = list("patient"=1:1e4, "step"=1:23, "a"=c("time","state"))))
mapply( function(z, b)
{ barplot(xlab = "months elapsed", ylab="patients in cohort", border = NA, 
          t(array(do.call('rbind', lapply(times, function(t) 
            tabulate(unlist( z[, ,"state"][
              cbind(1:dim(z)[[1]], apply(z[,,"time"], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
          )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))), 
          col=colors, beside=F, main=b) }, zl, names(zl))
legend("bottomleft", legend = rev(names(statetypes)), fill=rev(colors), cex=0.7)
mtext("Months elapsed", 1, outer=T, cex=0.7)
mtext("Patients in cohort", 2, outer=T, cex=0.7)
mtext("Status of RIF-R cohort over time", side=3, font=2, outer=T)

# can also look at above differences for other subsets e.g. +- moxi

# could change smear to 0 when treated (when on any treatment) if want it to reflect reduced transmission during that time. 
# As currently coded, need to combine TB status and (original, for the current or pendinng episode) smear status to judge current infectiousness. 


# redefining from cohort file, for single-patient course:
# but note that need to keep characteristics here, and call this function with an include of all pts at risk for the characteristic, not just those who have it at baseline (e.g. when talking about adr)
# will return time in states of interest, and final tbstate
time.in.state <- function(patiententry, states=8, characteristics=c(), characteristicvals=1, cutofftime=100*12, carryforward=F)
{
  elapsedtimes <- patiententry["eventtime",2:ncol(patiententry)] - patiententry["eventtime",1:((ncol(patiententry))-1)] 
  elapsedtimes[patiententry["eventtime",2:ncol(patiententry)]>cutofftime] <- cutofftime - (patiententry["eventtime",1:(ncol(patiententry)-1)])[patiententry["eventtime",2:(ncol(patiententry))]>cutofftime]
  elapsedtimes[elapsedtimes<0] <- 0
  indices <-  patiententry["TBstate",1:(ncol(patiententry)-1)] %in% states
  if(length(characteristics)==1) indices <- (indices & patiententry[characteristics,1:(ncol(patiententry)-1)]==characteristicvals)
  if(length(characteristics)>1) indices <- indices & apply(patiententry[characteristics,1:(ncol(patiententry)-1)]==characteristicvals, 2, all)
  t <- sum(elapsedtimes*indices)
  if (carryforward) 
    if (patiententry["TBstate",ncol(patiententry)] %in% states)
      t <- t + max(0, cutofftime-patiententry["eventtime",ncol(patiententry)])
  return(c(t, patiententry["TBstate",ncol(patiententry)]))
}

# tally up time alive: # need to cut all scenarios off at the same point, say 3 years - no, need a longer time window, see below:
followyears <- 20
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=1:7, cutofftime=12*followyears, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,summary))
# why is survival staying same or going down with novelrr for RIF==1 and cutofftime=36?? (same for MXOI==0) need to recheck with larger iniital sim **
## one possibility: I model low mortality during treatment, so longer treatment courses allow less opportunity for TB death. Esp with a cutoff at 3 years, there's no opportunity for death after unsuccessful conventional MDR-TB treament, but some do die of TB in the period after novel RR TB treatement (even though ultimately their risk is lower)
# but I'm seeing this to some extent even for very long cutofftimes -- novelrr doesn't reduce mortality. maybe bc most get hrze, so most mortality is before mdr treament starts, and of those who get to mdr treatment, there's still the efect that longer regimen defers relapse/failure mortality
# for MOXI==1, the lack of benefit of pantb is expected. 
#(there was also a problem with remaining diagnosed at the last step, which I've hopefuly fixed 12/30)
# for RIF==0 at baseline, there's a slight mortality increase with novel RR in 20181218 -- just chance? **

# incremental months of life per patient: 
lapply(l, function(x) apply(x[,1,] - l$baseline[,1,],2,mean))
# plot years of life gained within first [3] years:
par(mfrow=c(1,1), mar=c(4,4,3,1), oma=c(1,1,1,1))
boxplot( lapply(l, function(x) apply(x[,1,] - l$baseline[,1,],2,mean)), main=paste0("Incremental months of life gained per patient,\nover ",followyears," years after TB onset"))

# time on treatment (of some kind)
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(statetypes$treating, statetypes$treating_adr), cutofftime=12*followyears, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,summary))
boxplot(lapply(l, function(x) x[,1,1]), main="Total time on treatment per RIF-R patient", col=colors, ylim=c(0,20), ylab="months")
boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main="Average time on treatment for RIF-R patients", col=colors, ylim=c(0,20), ylab="months")

# time deceased i.e. YLL (compare to time alive, should see inverse)
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=8, cutofftime=12*followyears, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,mean))


# cures, as opposed to death or still on treatment at xx months:
t <- 36; 
still.in.state <- function(patiententry, states=8, time=36)
{
  # look at last entry with time < time, is it in states?
  maxed <- (max(patiententry["eventtime",])>t)
  return(c(
  patiententry["TBstate", ifelse(maxed, which.max(patiententry["eventtime",]>t)-1, length(patiententry["eventtime",]))] %in% states,
  maxed))
}
l <- loutcomeboot(individualoutcomefunction =still.in.state, states=7, time=t, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main=paste0("Proportion with cure (vs ongoing TB or death) ",t/12," years after TB onset"), col=colors, ylim=c(0,1))


# deaths within 5 years of TB onset:
t <- 120
l <- loutcomeboot(individualoutcomefunction =still.in.state, states=8, time=t, 
                  simoutput=impact, c=c, include=(c$RIF==1&c$RxHist==0), copies = 3, desiredsize = 1e5)
boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main=paste0("Mortality within ",t/12," years of TB onset, proportion of cohort"), col=colors, ylim=c(0,1))
## ** need to figure out why no mortality benefit for novelRR (here, RxHist and RIF do worse, and for new RIF there's no difference)
## are outcomes wrong? make.recurrence.matrix()[c("MDR, FQ-S", "MDR, FQ-R", "(ZE)","BPaMZ", "BPaM", "BPaZ"),c("6","18")] looks okay
## we saw above that we don't reduce YLL either.
## don't expect to reduce time to RR diagnosis (which in many cases will be long, with one or more (ZE) treatments first and taking 2 yrs or more),
## and once diagnosed as RR, expect TB mortality risk among ~10% after 18 months, vs among ~3% after 6 months --> 3% + 10%*3% over ~24 months,
## so by 24 months we should be seeing novelrr mortality benefit.

chisq.test(
  
  rbind(c(
table(impact$baseline[(1:50)[c$RIF==1],"TBstate",7,])
),c(
table(impact$novelrr[(1:50)[c$RIF==1],"TBstate",7,])
))

)
unlist(statetypes)

# time to cure, if cured (where cure is assumed to happen when treatment stops):
par(mar=c(3,3,3,1))
time.and.courses.to.cure <- function(patiententry) # also include number of treatment courses
{  c(
  ifelse(patiententry["TBstate",ncol(patiententry)]==statetypes$cured, patiententry["eventtime",which.max(patiententry["TBstate",]==statetypes$cured)], NA),
  sum(patiententry["eventtype",]==eventtypes$treatmentstart)  
  )
}
l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
                   simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,summary))
boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRR-TB patients only")
l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
                   simoutput=impact, c=c, include=(c$RIF==0), copies = 3, desiredsize = 1e5)
boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRS-TB patients only")


# total infectious time:
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=1:6, cutofftime=12*5, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e4)
lapply(l, function(x) apply(x[,1,],2,summary))
boxplot(lapply(l, function(x) x[,1,1]), main="Infectious time, RR-TB cases")



# estimated reduction in infectious time and force of infection (not accounting for infectiousness during treatment, changes in a person's infectiousness over time, etc): 
reduction <- lapply(outcomes[2:6], function(x)
  (apply(outcomes$baseline, 4, function(y) (time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"])) - 
     apply(x, 4, function(y) (time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"]))) /
    apply(outcomes$baseline, 4, function(y) time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12)["Mean"]))
b <- boxplot(reduction, main="Proportion of total infectious time averted, compared to baseline", col=colors[2:6])
text(1:5,unlist(lapply(reduction,median))+0.008, paste0(round(unlist(lapply(reduction,median))*100,1),"%"))

## Tracking acquired resistance
# number of acquired resistance events, overall
par(mfrow=c(3,2)); lapply(lapply(outcomes, function(x) apply(x, 4, function(y) sum(y[,"eventtype",]==eventtypes$treatmentstart&y[,"TBstate",]==statetypes$treating_adr))), 
                          hist)
## overall a little udnerunder 1% incidence in cohort
## happening more for baseline that for other scenarios, because of high fq-r acquisition with current rr regimen?

# number of acquired FQ events:
## increased moxi R
lapply(outcomes, function(x) apply(x, 4, function(y) sum(y[,"eventtype",2:dim(y)[[3]]]==eventtypes$treatmentstart & y[,"MOXI",2:dim(y)[[3]]]>y[,"MOXI",1:(dim(y)[[3]]-1)])))
# any moxi R, from none
par(mfrow=c(1,1))
boxplot(lapply(outcomes, function(x) apply(x, 4, function(y) sum(y[,"eventtype",2:dim(y)[[3]]]==eventtypes$treatmentstart & y[,"MOXI",2:dim(y)[[3]]]>0 & y[,"MOXI",1:(dim(y)[[3]]-1)]==0))),
        main="Instances of newly acquired moxifloxacin resistance,\nper 1000 patient cohort", col=colors)

# number of acquired RR events:
boxplot(lapply(outcomes, function(x) apply(x, 4, function(y) sum(y[,"eventtype",2:dim(y)[[3]]]==eventtypes$treatmentstart & y[,"RIF",2:dim(y)[[3]]]>y[,"RIF",1:(dim(y)[[3]]-1)]))), col=colors)

# number of acquired B/PA events:
boxplot(lapply(outcomes, function(x) apply(x, 4, function(y) sum(y[,"eventtype",2:dim(y)[[3]]]==eventtypes$treatmentstart & (y[,"BDQ",2:dim(y)[[3]]]>y[,"BDQ",1:(dim(y)[[3]]-1)] | y[,"PA",2:dim(y)[[3]]]>y[,"PA",1:(dim(y)[[3]]-1)])))), 
        main="Instances of acquired bedaquiline or pretonamid resistance,\nper 1000 patient cohort", col=colors)


# infectious time with drug resistance: 
par(mfrow=c(3,2), oma=c(1,2,3,1))
layout(array(c(1,1,3,2,2,4), dim=c(3,2)))
boxplot(lapply(lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("RIF"), cutofftime = 12*4))), function(z) z["Mean",]), main="Rifampin", col=colors, ylim=c(0,1))
boxplot(lapply(lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("MOXI"), cutofftime = 12*4))), function(z) z["Mean",]), main="Moxifloxacin", col=colors, ylim=c(0,1))
boxplot(lapply(lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("BDQ"), cutofftime = 12*4))), function(z) z["Mean",]), main="Bedaquiline", col=colors, ylim=c(0,0.2))
boxplot(lapply(lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("PA"), cutofftime = 12*4))), function(z) z["Mean",]), main="Pretomanid", col=colors, ylim=c(0,0.2))
mtext("Infectious time with drug resistance, among initial patient cohort over 3 years", side=3, font=2, outer=T)
mtext("Average infectious months with resistant strain, per patient in total cohort", side=2, cex=0.9, outer=T)
# lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("MOXI"), cutofftime = 12*4)))
# lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("BDQ"), cutofftime = 12*4)))
# lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("PA"), cutofftime = 12*4)))
# 


# *** adapt this for low/high ADR:
a1 <- data.frame(lapply(lapply(outcomes_nodelays, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("RIF"), cutofftime = 12*4))), function(z) z["Mean",]))
a2 <- data.frame(lapply(lapply(outcomes_nodelays, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("MOXI"), cutofftime = 12*4))), function(z) z["Mean",]))
a3 <- data.frame(lapply(lapply(outcomes_nodelays, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("BDQ"), cutofftime = 12*4))), function(z) z["Mean",]))
a4 <- data.frame(lapply(lapply(outcomes_nodelays, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("PA"), cutofftime = 12*4))), function(z) z["Mean",]))
a5 <- data.frame(lapply(lapply(outcomes_nodelays, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, cutofftime = 12*4))), function(z) z["Mean",]))
a1$drug <- "Rifampin"
a2$drug <- "Moxifloxacin"
a3$drug <- "Bedaquiline"
a4$drug <- "Pretomanid"
a5$drug <- "Any TB"
a5m <- melt(a5, id.vars=c("drug"))        
a <- bind_rows(a1,a2,a3,a4)
am <- melt(a, id.vars=c("drug"))
am$fraction <- am$value/rep(a5m$value,4)
am$params <- "Low ADR"
b1 <- data.frame(lapply(lapply(outcomes_highADR, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("RIF"), cutofftime = 12*4))), function(z) z["Mean",]))
b2 <- data.frame(lapply(lapply(outcomes_highADR, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("MOXI"), cutofftime = 12*4))), function(z) z["Mean",]))
b3 <- data.frame(lapply(lapply(outcomes_highADR, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("BDQ"), cutofftime = 12*4))), function(z) z["Mean",]))
b4 <- data.frame(lapply(lapply(outcomes_highADR, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, characteristics = c("PA"), cutofftime = 12*4))), function(z) z["Mean",]))
b5 <- data.frame(lapply(lapply(outcomes_highADR, function(x) apply(x, 4, function(y) time.in.state(states=c("undiagnosed", "diagnosed", "failed"), course = y, cutofftime = 12*4))), function(z) z["Mean",]))
b1$drug <- "Rifampin"
b2$drug <- "Moxifloxacin"
b3$drug <- "Bedaquiline"
b4$drug <- "Pretomanid"
b5$drug <- "Any TB"
b5m <- melt(b5, id.vars=c("drug"))        
b <- bind_rows(b1,b2,b3,b4)
bm <- melt(b, id.vars=c("drug"))
bm$fraction <- bm$value/rep(b5m$value)
bm$params <- "High ADR"
require(dplyr)
require(data.table)
mm <- bind_rows(am,bm)
# mm <- melt(m, id.vars=c("params","drug"))        
colnames(mm) <- c("Drug", "Scenario", "Infectious_Time", "Normalized", "Parameters")
require(ggplot2)
ggplot(data = mm, aes(x=Scenario, y=Infectious_Time)) + 
  geom_boxplot(aes(fill=Parameters)) + facet_wrap( ~ Drug, ncol=2) + 
  ggtitle("Drug resistance transmission potential") + theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Average months infectious with TB resistant to the specified drug")
ggplot(data = mm, aes(x=Scenario, y=Normalized)) + 
  geom_boxplot(aes(fill=Parameters)) + facet_wrap( ~ Drug, ncol=2) + 
  ggtitle("Drug resistance transmission potential") + theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Fraction resistant to the specified drug, of all TB transmitted by patient cohort ")





outcomes_nodelays <- outcomes
## then can add back treatment delays etc:
# params <- as.numeric(allparams[,2]); names(params) <- allparams[,1]
params["DSTdelay"] <- 0; params["DSTloss"] <- 0; params["Regimendelay"] <- 0.5; params["Regimenloss"] <- 0.1
# penalizing separate regimens, but not the performance of DST
outcomes <- list()
outcomes$baseline <- modelcourse(scenario = "0", cohort, params,reps = reps)
outcomes$novelrr <- modelcourse(scenario = "1a", cohort, params,reps = reps)
outcomes$novelrr_xxdr <- modelcourse(scenario = "1b", cohort, params,reps = reps)
outcomes$bpamz4 <- modelcourse(scenario = "2a", cohort, params,reps = reps)
outcomes$bpamz6 <- modelcourse(scenario = "2b", cohort, params,reps = reps)
outcomes$fullnovel <- modelcourse(scenario = "5", cohort, params,reps = reps)
outcomes_regimendelays <- outcomes

# rerun outcomes from above

params <- as.numeric(allparams[,2]); names(params) <- allparams[,1]
params["DSTdelay"] <- 0.5; params["DSTloss"] <- 0.1; params["Regimenloss"] <- params["Regimendelay"] <- 0
# penalizing separate regimens, but not the performance of DST
outcomes <- list()
outcomes$baseline <- modelcourse(scenario = "0", cohort, params,reps = reps)
outcomes$novelrr <- modelcourse(scenario = "1a", cohort, params,reps = reps)
outcomes$novelrr_xxdr <- modelcourse(scenario = "1b", cohort, params,reps = reps)
outcomes$bpamz4 <- modelcourse(scenario = "2a", cohort, params,reps = reps)
outcomes$bpamz6 <- modelcourse(scenario = "2b", cohort, params,reps = reps)
outcomes$fullnovel <- modelcourse(scenario = "5", cohort, params,reps = reps)
outcomes_DSTdelays <- outcomes


# get outcomes with no delays/ltfu and with high adr_bpamz==0.01, then add to acq res boxplot
params["DSTdelay"] <- 0; params["DSTloss"] <- 0; params["Regimenloss"] <- params["Regimendelay"] <- 0
params["adr_bpamz"]<-0.01
# penalizing separate regimens, but not the performance of DST
outcomes <- list()
outcomes$baseline <- modelcourse(scenario = "0", cohort, params,reps = reps)
outcomes$novelrr <- modelcourse(scenario = "1a", cohort, params,reps = reps)
outcomes$novelrr_xxdr <- modelcourse(scenario = "1b", cohort, params,reps = reps)
outcomes$bpamz4 <- modelcourse(scenario = "2a", cohort, params,reps = reps)
outcomes$bpamz6 <- modelcourse(scenario = "2b", cohort, params,reps = reps)
outcomes$fullnovel <- modelcourse(scenario = "5", cohort, params,reps = reps)
outcomes_highADR <- outcomes
