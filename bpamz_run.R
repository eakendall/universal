setwd("C:/Users/ekendal2/OneDrive - Johns Hopkins University/Research/universal regimen/universal")
source("bpamz_cohort_20181218.R")

date <- "20181218"

cohortsize <- 5e7 # 1e7 is barely enough to capture all but the rarest combinations at least once. Really using this to get freqs.
# make the cohort
cohort <- make.cohort(params=params, patientvars = patientvars, N=cohortsize)
c <- as.data.frame(table(cohort),stringsAsFactors = F)
c[rev(order(c$Freq)),]
c[] <- lapply(c, function(x) as.numeric(x))


distincts_with_freqs <- c[rev(order(c$Freq)),] # rather than modeling based on frequency in the population, model each patient type the same number of times and then bootstrap (if we increase the number of different types, may want to model the common types more than the rarer types)
c <- distincts_with_freqs %>% filter(MOXI==1|partialmoxi==0) # an imposible combo as coded
save(c, file="cohort.Rdata")
load("cohort.Rdata")
tail(c %>% filter(INH==1|RIF==0)) # Rif monos essentially don't exist in India, so checking that other combos aren't missing
nrow(c %>% filter(Freq>0)) # number of patient types to be run, ~140 for India
c <- c %>% filter(Freq>0)
reps <- 100 #1e3 # how many reps of each patient type to run? (will then sample Freq of each to recreate original cohort or subset thereof)
# 1e3 already gives a >1GB list for 3 scenarios, so may want to filter the cohort (e.g. to RR of FQ-R only) before the next step.  

impact <- list()
impact$baseline <- modelcourse(scenario = "0", c, params, reps = reps)
impact$novelrr <- modelcourse(scenario = "1a", c, params, reps = reps)
impact$novelpantb <- modelcourse(scenario = "3", c, params, reps = reps)
saveRDS(object = impact, file = paste0("impact.",date,".RDS"))

dst <- list()
dst$noxxdr <- modelcourse(scenario = "3", c, params, reps = reps)
dst$stepxxdr <- modelcourse(scenario = "4", c, params, reps = reps)
dst$fullxxdr <- modelcourse(scenario = "5", c, params, reps = reps)
saveRDS(object = dst, file = paste0("dst.",date,".RDS"))


# can look at this equally-weighted "cohort":
lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",20]==7)))) 
lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",20]==8)))) 
lapply(impact, FUN = function(x) summary(apply(x, 4, function(y) mean(y[,"TBstate",12]==1)))) 
# and particular subsets
lapply(impact, FUN = function(x) summary(apply(x[c$RIF==1,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) #RR
lapply(impact, FUN = function(x) summary(apply(x[c$MOXI==0,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) 
lapply(impact, FUN = function(x) summary(apply(x[c$RIF==1&c$MOXI==0,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) 
lapply(impact, FUN = function(x) summary(apply(x[c$RIF==1&c$MOXI==1,,,], 4, function(y) mean(y[,"TBstate",10]==8)))) 

#** To do: add worse HRZE outcomes for PZA-R. Would need to differentiate the active regimens for (ZE) to +-Z(E)

# But really we want to bootstrap to construct the cohort whose outcomes will be modeled: 
#  for each element of the above lists, for each ith entry in 1st dimension (which corresponds to the rows of c), 
# sample dims 2 and 3 (across dim 4) based on ith frequncy in c$Freq, and abind the result into a course (3 dims), 
# and do this copies times to get the 4th dim (variability between simulations of the same cohort).

include <- c$RIF==1 # define a sub-cohort of interest.  
copies <- 10 # how many bootstrap copies of the whole cohort to create?
desiredsize <- 1e5 # how big to make each bootstrapped copy?

# sample  by cohort weight from the patient type number (dimension 1 in modelcourse output)
# and for each, sample randomly from the reps of that patient type (i.e. across dimension 4)
# make copies of these bound across 4th dimension to get multiple realizastions of the cohort:
boot <- function(course)
  {
  a <- abind(lapply(1:copies, FUN=function(x) {
  bootindices <- sample((1:nrow(c))[include], desiredsize, replace=T, prob=c$Freq[include]/sum(c$Freq[include]))
  course[bootindices,,,sample(1:reps,1)]}), 
  along=4)
  return(a)
}
final <- (boot(impact$baseline)) # this may require lots of memory
str(final)
str(impact$baseline) # these have the same dimension types (through diff length for 1st and 4th)

# or, to save memory, calculate the outcome(S) of interest for each patient type (~100 types, versus >10k in a cohort)
# individualoutcomefunction takes one patient entry matrix (characteristic x timestep) from course (e.g. course[1,,,1])
  # and generates a *vector* of outcomes (ideally named)
outcomeboot <- function(individualoutcomefunction, course, include=rep(1, dim(course)[[1]]), copies=10, desiredsize=1e4)
{
    course_outcomes <- apply(X=course, FUN=individualoutcomefunction, MARGIN = c(1,4))
    a <- abind(lapply(1:copies, FUN=function(x) {
    bootindices <- sample((1:nrow(c))[include], desiredsize, replace=T, prob=c$Freq[include]/sum(c$Freq[include]))
    course_outcomes[,bootindices,sample(1:reps,1)]}), 
    along=3)
  return(a)
}

# example of an individualoutcomefunction:
fun <- function(patiententry)
  {
    return(quantile(patiententry["eventtime",],c(.05,.5,.95)))
}
patiententry <- impact$baseline[1,,,1]
fun(patiententry)
#or,
lastalive <- function(patiententry)
{
  return(max((data.frame(t(patiententry)) %>% filter(TBstate<8))["eventtime"]))
}
lastalive(patiententry)

# and bootstrap from the results of that outcome funcntion
# for one course:
o <- outcomeboot(fun, course=impact$baseline, include=include, copies=10, desiredsize=1e5)  
str(o) # here, 1st dimension is new reflecting the outcome(x) of interest, with 2nd and 3rd reflecting variability. 2nd is what was the 1st in include or boot, and 3rd is what was the 4th. (i.e. the 2nd and 3rd dims of outcomeboot are just patients, differing across second dim and replicated across 3rd.)
# or for a list of courses for different scenarios
loutcomeboot <- function(individualoutcomefunction, simoutput, include=rep(1, dim(l[[1]])[[1]]), copies=10, desiredsize=1e4) 
  {
    freqweights <- impact$baseline[,"Freq",1,1]
    l <- lapply(simoutput, function(y) {  
     course_outcomes <- apply(X=y, FUN=individualoutcomefunction, MARGIN = c(1,4))
    abind(lapply(1:copies, FUN=function(x) {
      bootindices <- sample((1:nrow(c))[include], desiredsize, replace=T, prob=freqweights[include]/sum(freqweights[include]))
      course_outcomes[,bootindices,sample(1:reps,1)]}), 
      along=3)} )
    return(l)
  }
  

l <- loutcomeboot(fun, impact, include=include, copies=10, desiredsize=1e5)
str(l)
# median of some component of the outcome vector, across the cohort:
lapply(l, function(x) apply(x["5%",,], 2, median))
lapply(l, function(x) summary(apply(x["5%",,], 2, median)))
lapply(l, function(x) apply(x["5%",,], 2, summary))
# can see that these are varying a lot between the 3 copies
lapply(l, function(x) apply(x["50%",,], 2, summary))
#median event times vary less than 5th precentiles between the bootstrapped cohorts, as expected.
lapply(l, function(x) summary(apply(x["50%",,], 2, median)))
#this shows the median event time going down (as expected because shorter time to treatment and higher and faster cure, through not all that interesting)

#  diagnosed by step two:
fun <- function(patiententry) { return(c(patiententry["TBstate",2]==2, NA))}
l <- lapply(impact, function(x) outcomeboot(individualoutcomefunction = fun, 
                                            course=x, include=include, copies=50, desiredsize=1e6)[1,,])  
lapply(l, function(x) summary(apply(x, 2, mean)))
# lots of stochasticity, but once I increase to 50 copies and size 1e6 the results are consistent between reps

####
# plot states over time? this is by model step rather than time
par(mfrow=c(1,1))
barplot(array(unlist(apply(impact$baseline[,"TBstate",2:20,], 2, tabulate)), dim = c(8,19)), beside=F, col = rainbow(8), legend.text = names(statetypes))
    
    
# snapshots by time interval: 
require(RColorBrewer)
colors <- palette(brewer.pal(length(statetypes), name = "Set2"))
# require(colorspace); colors <- palette(rainbow_hcl(n = length(statetypes))); colors <- choose_palette()
    
# plot average states after first treatment attempt, for the equal-wieghted cohort:
par(mfrow=c(1,1), oma=c(0,0,1,1))
barplot(array(unlist(lapply(impact, FUN = function(x) tabulate(x[,"TBstate",8,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Distribution of states after first round of treatment", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)

# do the same but for the "include"+reweighted subset cohort of interest:
par(mfrow=c(1,1), oma=c(0,0,1,1))
barplot(array(unlist(lapply(impact, FUN = function(x) tabulate(x[,"TBstate",8,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
        main="Distribution of states after first round of treatment", yaxt='n')
legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1)

    
# Differences compared to baseline: 
par(mfrow=c(1,1), oma=c(0,0,1,1), mar=c(3,3,3,1))
boxplot(lapply(outcomes, function(y) apply(y-outcomes$baseline, 4, function(x) mean(x[,"TBstate",4]==statetypes$cured))),
        main="Incremental proportion of cohort cured after first treatment,\ncompared to baseline", col=colors)
boxplot(lapply(outcomes[2:6], function(y) apply(y, 4, function(x) mean(x[,"TBstate",4]==statetypes$cured)-mean(outcomes$baseline[,"TBstate",4,]==statetypes$cured))),
        main="Incremental proportion of cohort cured after first treatment,\ncompared to baseline average", col=colors[2:6])
    
    # showing with and without regimen delays:
    boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), 
            boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1)) #invisible boxes
    boxplot(df[which(df$id=="Good"), -1], xaxt = "n", add = TRUE, boxfill="red", boxwex=0.25, 
            at = 1:ncol(df[,-1]) - 0.15) #shift these left by -0.15
    boxplot(df[which(df$id=="Bad"), -1], xaxt = "n", add = TRUE, boxfill="blue", boxwex=0.25,
            at = 1:ncol(df[,-1]) + 0.15) #shift these right by +0.15
    
    
    a <- data.frame(lapply(outcomes_nodelays[2:6], function(y) apply(y, 4, function(x) mean(x[,"TBstate",4]==statetypes$cured)-mean(outcomes$baseline[,"TBstate",4,]==statetypes$cured))))
    a$params <- "No delays"
    b <- data.frame(lapply(outcomes_regimendelays[2:6], function(y) apply(y, 4, function(x) mean(x[,"TBstate",4]==statetypes$cured)-mean(outcomes$baseline[,"TBstate",4,]==statetypes$cured))))
    b$params <- "Regimen delays"
    c <- data.frame(lapply(outcomes_DSTdelays[2:6], function(y) apply(y, 4, function(x) mean(x[,"TBstate",4]==statetypes$cured)-mean(outcomes$baseline[,"TBstate",4,]==statetypes$cured))))
    c$params <- "DST delays"
    require(dplyr)
    require(data.table)
    m <- bind_rows(a,b,c)
    mm <- melt(m, id.vars="params")        
    colnames(mm) <- c("Parameters", "Scenario", "Cureincrease")
    mm$Parameters <- factor(mm$Parameters, levels=c("No delays", "Regimen delays", "DST delays"))
    require(ggplot2)
    ggplot(data = mm, aes(x=Scenario, y=Cureincrease)) + 
      geom_boxplot(aes(fill=Parameters)) + 
      ggtitle("Increase in proportion of cohort cured after first treatment,\ncompared to baseline scenario") + theme(plot.title = element_text(hjust = 0.5)) + 
      ylab("Incremental proportion cured")
    
    
    # as a heat map (don't group by state)
    par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
    times <- 0:30
    whichpts<-1:100
    mapply(FUN = function(a, b) { y <- a[whichpts,,,1]; 
    image( 
      do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
        ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))])),
      col=colors, xaxt='n', yaxt='n', main=paste0(length(whichpts)," representative patients, ",b));
    axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
    },
    outcomes, names(outcomes))
    legend("bottomright", legend = (names(statetypes)), fill=(colors), cex=1)
    mtext(text = "months", side = 1, outer=T, line=-1)
    mtext(text = "Each horizontal line represents one patient's course", side = 2, outer=T, line=-1)
    
    par(mfrow=c(1,1))
    times <- 0:30
    whichpts<-1:100
    mapply(FUN = function(a, b) { y <- a[whichpts,,,1]; 
    image( 
      do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
        ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))])),
      col=colors, xaxt='n', yaxt='n', main=paste0(length(whichpts)," representative patients, ",b));
    axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
    },
    outcomes, names(outcomes))
    legend("bottomright", legend = (names(statetypes)), fill=(colors), cex=1)
    mtext(text = "months", side = 1, outer=T, line=-1)
    mtext(text = "Each horizontal line represents one patient's course", side = 2, outer=T, line=-1)
    
    # for one run, sorted:
    par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
    times <- seq(0,30,by=1)
    y <- outcomes$baseline[,,,2]
    barplot(xlab = "months elapsed", ylab="patients in cohort", border = NA, 
            t(array(do.call('rbind', lapply(times, function(t) 
              tabulate(unlist( y[, "TBstate",][
                cbind(1:nrow(cohort), apply(y[,"eventtime",], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
            )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))), 
            col=colors, beside=F)
    legend("topright", legend = rev(names(statetypes)), fill=rev(colors), cex=0.8)
    
    # # for one run, sorted but in heatmap mode:
    # msort <- function(m)
    # {return(m[do.call(order, c(decreasing = FALSE, data.frame(m[,1:ncol(m)]))),])}
    # 
    # par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
    # times <- seq(0,30,by=0.1)
    # y <- outcomes$baseline
    # m <- do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
    #   ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))]))
    # image(t(msort(t(m))), col=colors, xaxt='n', yaxt='n')
    # legend("bottomleft", legend = rev(names(statetypes)), fill=rev(colors), cex=0.8)
    # axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
    # mtext(text = "months", side = 1, line=2)
    # mtext(text = "cohort", side = 2, line=2)
    
    # for multiple scenarios, bar plots, :
    par(mfrow=c(2,3), mar=c(4,4,1,1), oma=c(1,1,2,1))
    times <- seq(0,24,by=1)
    mapply(function(a, b) { 
      y <- a[,,,1];
      barplot(xlab = "", ylab="", border = NA, 
              t(array(do.call('rbind', lapply(times, function(t) 
                tabulate(unlist( y[, "TBstate",][
                  cbind(1:nrow(cohort), apply(y[,"eventtime",], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
                # abind(1:nrow(cohort), apply(y[,"eventtime",,], c(1,3), function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
              )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))), 
              col=colors, beside=F, main=b)
    }, outcomes, as.list(names(outcomes)))
    legend("bottomleft", legend = rev(names(statetypes)), fill=rev(colors), cex=1)
    mtext("Months elapsed", 1, outer=T, cex=0.7)
    mtext("Patients in cohort", 2, outer=T, cex=0.7)
    mtext("Status of cohort over time", side=3, font=2, outer=T)
    
    # differences are too small to visualize clearly across this whole cohort. 
    # will want to look at incremental differences. 
    
    # {m <- do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
    #   ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))]))
    # image(t(msort(t(m))), col=colors, xaxt='n', yaxt='n')
    # }) 
    
    
    # #ultimate deaths
    # s <- do.call('rbind', lapply(outcomes, function(x) tabulate(x[,"TBstate",dim(x)[3],]))); colnames(s) <- names(statetypes); 
    # s
    # # differences will approximate deaths averted, but totals are at some arbitrary cutoff, plus there's the problem that ending time is later for baseline
    
    # tally up time alive: # need to cut all scenarios off at the same point, say 3 years:
    (c <- lapply(outcomes, function(y) apply(y, 4, function(x) time.in.state(states = c("undiagnosed", "diagnosed", "treating", "treating_adr", "failed", "pendingrelapse", "cured"), course = x, cutofftime = 3*12, carryforward = T)) ))
    unlist(lapply(c,function(x) mean(x["Mean",]))); unlist(lapply(c,function(x) sd(x["Mean",])))
    # incremental months of life per patient: 
    unlist(lapply(c,function(x) summary(x["Mean",]-mean(c$baseline["Mean",])))); unlist(lapply(c,function(x) sd(x["Mean",])))
    # plot years of life gained within first [3] years:
    par(mfrow=c(1,1), mar=c(4,4,3,1), oma=c(1,1,1,1))
    boxplot( lapply(c[2:6],function(x) x["Mean",] - mean(c$baseline["Mean",])), main="Incremental months of life gained per patient,\nover 3 years after TB onset", col=colors[2:6])
    
    
    # time on treatment (of some kind)
    (c <- lapply(outcomes, function(y) apply(y, 4, function(x) time.in.state(states = c("treating", "treating_adr"), course = x, cutofftime = 3*12, carryforward = T)) ))
    # median and IQR:
    rbind(unlist(lapply(c,function(x) mean(x["1st Qu.",]))), unlist(lapply(c,function(x) mean(x["Median",]))), unlist(lapply(c,function(x) mean(x["3rd Qu.",]))))
    unlist(lapply(c,function(x) mean(x["Mean",]))); unlist(lapply(c,function(x) sd(x["Mean",])))
    boxplot(lapply(c,function(x) x["Mean",]), main="Average time on treatment (months)", col=colors, ylim=c(0,6))
    
    # time deceased (compare to time alive, should see inverse - not exact but pretty close)
    (c <- lapply(outcomes, function(y) apply(y, 4, function(x) time.in.state(states = c("deceased"), course = x, cutofftime = 3*12, carryforward = T)) ))
    unlist(lapply(c,function(x) mean(x["Mean",]))); unlist(lapply(c,function(x) sd(x["Mean",])))
    unlist(lapply(c,function(x) mean(c$baseline["Mean",] - x["Mean",]))); unlist(lapply(c,function(x) sd(x["Mean",])))
    
    
    # cures, as opposed to death or still on treatment after 2 rounds:
    lapply(outcomes, function(x) summary(apply(x[,"TBstate",10,]==statetypes$cured, 2, mean)))
    boxplot(lapply(outcomes, function(x) apply(x[,"TBstate",10,]==statetypes$cured, 2, mean)), main="Proportion cured within two rounds of treatment", col=colors)
    
    # cures, as opposed to death or still on treatment at xx years:
    t <- 36; 
    lapply(lapply(outcomes, function(y) apply(y, 4, function(x) mean((x[, "TBstate",][cbind(1:nrow(cohort), apply(x[,"eventtime",], 1, function(z) ifelse(max(z)>t, which.max(z>t)-1, length(z)) ))])==statetypes$cured))), summary)
    boxplot(lapply(outcomes, function(y) apply(y, 4, function(x) mean((x[, "TBstate",][cbind(1:nrow(cohort), apply(x[,"eventtime",], 1, function(z) ifelse(max(z)>t, which.max(z>t)-1, length(z)) ))])==statetypes$cured))), 
            main=paste0("Proportion with cure (vs ongoing TB or death) ",t/12," years after TB onset"), col=colors)
    # something's wrong here, lower than counting cures at a relatively early time point (2 rounds of treatment)
    
    
    # deaths within 3 years of TB onset:
    boxplot(lapply(outcomes, function(y) apply(y, 4, function(x) sum(x[,"eventtype",]==eventtypes$death & x[,"eventtime",]<36)/(N))), main=c("Mortality within three years of TB onset, proportion of cohort"), col=colors)
    
    
    
    # time to cure, if cured (where cure is assumed to happen when treatment stops):
    par(mar=c(3,3,3,1))
    boxplot(lapply(outcomes, function(x) apply(x, 4, function(y) mean(
      y[, "eventtime",][cbind(1:nrow(cohort), apply(y[,"TBstate",], 1, function(x) ifelse(mean(x==statetypes$cured)>0, which.max(x==statetypes$cured), NA)))], na.rm=T)
    )), main="Average time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured")
    # bpamz4 and fullnovel are faster than the others
    
    
    # time in infectious states (excluding on-treatment)
    lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(c("undiagnosed"), y))) # not much difference, somewhat longer for baseline
    lapply(outcomes, function(y) time.in.state(c("diagnosed"), y)) #none with parameters set to zero
    lapply(outcomes, function(y) time.in.state(c("failed"), y)) # median and 3rd quartile zero here, but means much different
    
    # # cutting off tails to reduce effects of outliers: 
    # limits <- 1:10
    # lapply(outcomes, function(y) time.in.state(c("undiagnosed"), y[,,limits]))
    # lapply(outcomes, function(y) time.in.state(c("diagnosed"), y[,,limits]))
    # lapply(outcomes, function(y) time.in.state(c("failed"), y[,,limits]))
    
    # total infectious time:
    # (l <- lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(c("undiagnosed", "diagnosed", "failed"), y))))
    # boxplot(lapply(l, function(y) y["Mean",]), main="Average infectious time (in months)")
    # using a cutoff time to avoid bias:
    summary(outcomes$bpamz4[,"eventtime",dim(outcomes$bpamz4)[3],])
    (l <- lapply(outcomes, function(x) apply(x, 4, function(y) time.in.state(c("undiagnosed", "diagnosed", "failed"), y, cutofftime = 4*12))))
    boxplot(lapply(l, function(y) y["Mean",]), main="Average infectious time (in months)", col=colors)
    
    
    ## with a cohort size of 5000 mostly DS, the random variation overwhelms the regimen-related, with the exception that baseline is most and bpamz6 is least as expected.
    ## when there's a regimen delay, bpamz4 can perform better than fullnovel here
    
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
    