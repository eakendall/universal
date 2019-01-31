date <- "20190109"
setting <- "SEA"
source("bpamz_cohort.R")
source("bpamz_result_utils.R")

load(file = paste0("impact.",date,".Rdata"))
c <- data.frame(impact$baseline[,,1,1])
# dst <- readRDS(file = paste0("dst.",date,".RDS"))

# #  average states after 1 and 2 treatment attemptfor the "include"+reweighted subset cohort of interest (here, RIF-R):
# par(mfcol=c(1,2))
# l1 <- loutcomeboot(function(patiententry) {patiententry[c("TBstate", "eventtime"),4] }, impact, c, include=(c$RIF==1), copies=100, desiredsize = 1e5)
# l1b <- loutcomeboot(function(patiententry) {patiententry[c("TBstate", "eventtime"),4] }, impact, c, include=(c$RIF==1), copies=100, desiredsize = 1e5,oncluster = T)
# barplot(array(unlist(lapply(l1, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
#         main="Status after 1st round of\ntreatment, RR-TB cases", yaxt='n')
# l2 <- loutcomeboot(step8status, impact, c, include=(c$RIF==1), copies=100, desiredsize = 1e4)
# barplot(array(unlist(lapply(l2, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
#         main="Status after 2ndround of\ntreatment, RR-TB cases", yaxt='n')
# legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1, bg="white")
# benefit for novelrrx over novelrr is more than the increase in Xpert coverage, because the new pts whose xpert increases most had longest TBmortality exposure
# additional benefit of pantB here is those whose TB or RR is missed by Xpert? ~10% for TB and another 5% rif


# #  proportion cured after 1st and 2nd round:
# par(mfrow=c(1,2), oma=c(0,0,1,1), mar=c(3,3,3,1))
l4 <- loutcomeboot(step4812cure, impact, c, include=(c$RIF==1), copies=50, desiredsize = 1e5)
save(l4, file=paste0("l4.",date,".Rdata"))
lapply(l4, function(y) mean(apply(y[,1,],2, mean))); lapply(l4, function(y) sd(apply(y[,1,],2, mean)))
l4s <- loutcomeboot(step4812cure, impact, c, include=(c$RIF==0), copies=50, desiredsize = 1e5)
save(l4s, file=paste0("l4s.",date,".Rdata"))
lapply(l4s, function(y) mean(apply(y[,1,],2, mean))); lapply(l4s, function(y) sd(apply(y[,1,],2, mean)))
lapply(l4s, function(y) mean(apply(y[,2,],2, mean)))
lapply(l4s, function(y) mean((y[,2,])[y[,1,]==0])) # among those not cured first round,cured second?
lapply(l4s, function(y) mean((y[,3,])[y[,2,]==0])) # among those not cured second round,cured third? (many are dead)

# total infectious time:
li <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*10, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 50, desiredsize = 1e5)
save(li, file=paste0("li.",date,".Rdata"))
liall <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*10, carryforward=T, 
                   simoutput=impact, c=c, include=(c$RIF<2), copies = 50, desiredsize = 1e5)
save(liall, file=paste0("liall.",date,".Rdata"))

lis <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*10, carryforward=F, 
                   simoutput=impact, c=c, include=(c$RIF==0), copies = 50, desiredsize = 1e5)
lipans <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*3, carryforward=T, 
                    simoutput=impact, c=c, include=(c$RIF==0&c$PZA==0&c$MOXI==0), copies = 50, desiredsize = 1e5)
lihr <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*10, carryforward=F, 
                       simoutput=impact, c=c, include=(c$RIF==0&c$INH==1), copies = 50, desiredsize = 1e5)
lapply(lihr, function(x) mean(apply(x[,1,],2,mean)))
lapply(lihr, function(x) sd(apply(x[,1,],2,mean)))
# larger novelpantb benefit as it should be, but not sure why int time is increasing (worsening) with novelrr coverage, should do better when rr acquired in second round, andn should otherwise have no effect.

lapply(li, function(x) mean(apply(x[,1,],2,mean))); lapply(li, function(x) sd(apply(x[,1,],2,mean)))
lapply(liall, function(x) mean(apply(x[,1,],2,mean))); lapply(liall, function(x) sd(apply(x[,1,],2,mean)))
lapply(lis, function(x) mean(apply(x[,1,],2,mean)))
lapply(lipans, function(x) mean(apply(x[,1,],2,mean)))
# issue here, (slightly) higher cure rates each round with novelpantb among rs, but longer infectious time vs baseline or novelrr. 
# shorter treatment duration means those unsuccessfully treated relapse sooner, therefore contribute moe inf time before a given cutoff.
# got rid of carryforward (cut off patients who reach the end of data still infecitous), but should also have option to follow for a given number of treatment cycles rather than a time? So among the tail with multiple recurrences, those who recur sooner get followed for less total time. 
# and why does novelrr affect them at all?  the few who started RR but acquired rif resistant (0.5% of all, 10% of recurrence) will have their poor outcomes cut by >2x, 
# and if each recurrence results in ~6 infectious months, then reducing average rec time by 6*0.5%*1/2 ~.01 months, not 0.1
load(file="moreimpact.20190109.Rdata")
# getting rid of effects of 
allimpact <- c(impact, moreimpact); rm(moreimpact)
lipans <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,5,9), cutofftime=12*4, carryforward=T, 
                       simoutput=allimpact, c=c, include=(c$RIF==0&c$PZA==0&c$MOXI==0), copies = 50, desiredsize = 1e4)
lapply(lipans, function(x) mean(apply(x[,1,],2,mean)))

# a recurrence extends infectiousnes sby ~50%, so 0.5% reccurnce with 
# so overall, so overall should improve by 0.1%
# and effect more than tripled by novelrrx because it ~triples the number who get good treatment the first round, vs 2nd or 3rd or beyond
boxplot(lapply(l, function(x) x[,1,1]), main="Total infectious time, RR-TB cases")
## also consider by smear status here?? (but note that smear pos at dx doesn't mean always smear pos)

# time on treatment (of some kind)
followyears <- 20
ltall <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(statetypes$treating, statetypes$treating_adr), cutofftime=12*followyears, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF<2), copies = 50, desiredsize = 1e5)
save(ltall, file=paste0("ltall.",date,".Rdata"))
lapply(ltall, function(x) mean(apply(x[,1,],2,mean)))
lapply(ltall, function(x) sd(apply(x[,1,],2,mean)))

# boxplot(lapply(l, function(x) x[,1,1]), main="Total time on treatment per RIF-R patient", col=colors, ylim=c(0,20), ylab="months")
# boxplot(lapply(l, function(x) apply(x[,1,], 2, mean)), main="Average time on treatment for RIF-R patients", col=colors, ylim=c(0,20), ylab="months")


# # boxplot(lapply(l4, function(y) apply(y[,1,]- apply(l4$baseline[,1,],2, mean), 2, mean)),
# #         main="Incremental proportion of RIF-R cohort cured after first treatment,\ncompared to baseline average", col=colors[2:6])
# require(reshape2); require(ggplot2)
# l4.m <- lapply(l4, function(x) melt(x, ))
# boxplot(lapply(l4, function(y) apply(y[,1,], 2, mean)), ylim=c(0,1),
#         main="Proportion of RIF-R cohort cured after first treatment", col=colors[2:6])
# boxplot(lapply(l4, function(y) apply(y[,2,], 2, mean)), ylim=c(0,1),
#         main="Proportion of RIF-R cohort cured after second treatment", col=colors[2:6])

# # can look at step4status for a different subset:
# l3 <- loutcomeboot(step4status, impact, c, include=(c$MOXI==1), copies=5, desiredsize = 1e4)
# #and then tabulate across the cohort and plot as before:
# par(mfrow=c(1,1), oma=c(0,0,1,1))
# barplot(array(unlist(lapply(l3, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
#         main="Status after 1st treatment,\nFQ-r patients:\n", yaxt='n')
# legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1, bg='white')
# #for Moxi R, Rif S patients, novel regimen looks same or slightly worse
# l31 <- loutcomeboot(step4status, impact, c, include=(c$MOXI==1&c$RIF==0), copies=5, desiredsize = 1e4)
# #and then tabulate across the cohort and plot as before:
# par(mfrow=c(1,1), oma=c(0,0,1,1))
# barplot(array(unlist(lapply(l31, FUN = function(x) tabulate(x[,1,]))), dim=c(length(statetypes),length(impact))), names.arg = names(impact), col=colors,
#         main="Status after 1st treatment,\nFQ-r, RIF-s patients:\n", yaxt='n')
# legend("right", legend = rev(names(statetypes)), fill=rev(colors), cex=1, bg='white')


# # as a heat map for status over time (don't group by state)
# par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
# times <- 0:30
# whichpts<-c$RIF==1 #** can change here, should change title to match
# mapply(FUN = function(a, b) { y <- a[whichpts,,,1]; 
# image( 
#   do.call('rbind', lapply(times, function(t) y[, "TBstate",][cbind(1:dim(y)[[1]], apply(y[,"eventtime",], 1, function(x) 
#     ifelse(max(x)>t, which.max(x>t) -1, length(x)) ))])),
#   col=colors, xaxt='n', yaxt='n', main=paste0(length(whichpts)," patient types, ",b));
# axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
# },
# impact, names(impact))
# legend("bottomright", legend = rev(names(statetypes)), fill=rev(colors), cex=0.9)
# # or better, could make this 100 reps of a common or important patient type:
# w <- which.max(c$RIF==1 & c$MOXI==1)
# mapply(FUN = function(a, b) { y <- a[w,,,1:100]; # * set number to plot here
# # mapply(FUN = function(a, b) { y <- a[w,,,]; 
# image( 
#   do.call('rbind', lapply(times, function(t) y["TBstate",,][cbind(apply(y["eventtime",,], 2, function(x) 
#     ifelse(max(x)>t, which.max(x>t) -1, length(x)) ), 1:dim(y)[[3]])])),
#   col=colors, xaxt='n', yaxt='n', main=paste0(reps," duplicated patients, ",b));
# axis(side = 1, labels = times, at = seq(0,1,length=length(times)))
# },
# impact, names(impact))
# legend("bottomright", legend = (names(statetypes)), fill=(colors), cex=0.8)
# mtext(text = "months", side = 1, outer=T, line=-1)
# mtext(text = "Each horizontal line represents one realization of the patient's course under the specified scenario", side = 2, outer=T, line=-1)

# # for one run, sorted:
# par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
# times <- seq(0,30,by=1)
# y <- loutcomeboot(function(x) x[c("eventtime","TBstate"),], impact, c, include=(c$RIF==1), copies=1, desiredsize = 1e5)$baseline
# z <- abind(y[,seq(1,45,by=2),1], y[,seq(2,46,by=2),1], along = 3, new.names = list("patient"=1:1e5, "step"=1:23, "a"=c("time","state")))
# barplot(xlab = "months elapsed", ylab="patients in cohort", border = NA, 
#         t(array(do.call('rbind', lapply(times, function(t) 
#           tabulate(unlist( z[, ,"state"][
#             cbind(1:dim(z)[[2]], apply(z[,,"time"], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
#         )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes)))), 
#         col=colors, beside=F)
# legend("topright", legend = rev(names(statetypes)), fill=rev(colors), cex=0.8)
# mtext("Status of RIF-R patients over time, in baseline scenario (no novel regimen)", side=3)

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
yl <- loutcomeboot(function(x) x[c("eventtime","TBstate"),], impact, c, include=(c$RIF==1), copies=1, desiredsize = 1e5)
zl <- lapply(yl, function(y) abind(y[,seq(1,35,by=2),1], y[,seq(2,36,by=2),1], along = 3, new.names = list("patient"=1:1e5, "step"=1:18, "a"=c("time","state"))))
save(yl, zl, file=paste0("fig1RRstatestatus.",date,".Rdata"))
maxtime <- 30; times <- seq(0,maxtime,by=1)
unlist(statetypes)
fullstatenames <- c(" \nTB not yet treated\n ", "shoulodn't show up", "On treatment", " \nOn treatment &\nacquired resistance",
                    "Treatment failed", "Will relapse", "Cured", "Deceased", "Relapsed TB")
plotorder <- c(1,5,9,6,3,4,7,8)
pdf("Fig2.pdf", width=10, height=7)
layout(matrix(c(1,3,2,4,5,5), ncol=3), widths = c(4,4,2))
par(mar=c(3,3,1,2), oma=c(2,2,2,2))
mapply( function(z, title)
{ statestatus <- t(array(do.call('rbind', lapply(times, function(t) 
      tabulate(unlist( z[, ,"state"][
        cbind(1:dim(z)[[1]], apply(z[,,"time"], 1, function(x) ifelse(max(x)>t, which.max(x>t)-1, length(x)) ))]), nbins = length(statetypes))
            )), dim=c(length(times), length(statetypes)), dimnames=list("time"=times, "state"=names(statetypes))))
  
  b <- barplot(xlab = "", ylab="", yaxt='n', xaxt='n', border = NA, space=0,
          statestatus[plotorder,], 
          col=colors, beside=F)
  mtext(title, side=3, font=2,line=0.5, cex=1)
  axis(side = 2, at = seq(0,1e5,by=2e4),labels = seq(0,100,by=20), cex.axis=0.8, las=2)
  axis(side = 1, at = b[seq(1,length(b),by=3)],labels = seq(0,maxtime,by=3), cex.axis=0.8)
  label <- statestatus[plotorder,maxtime]/1e3 >= 0.5 
  text(x = (b[length(b)]+c(0,3,rep(0,sum(label)-2))), y = (cumsum(statestatus[plotorder,maxtime])- statestatus[plotorder,maxtime]/2)[label],pos=4, xpd=NA,
       labels = paste0(round(statestatus[plotorder,maxtime]/1e3,1)[label], "%"), cex=0.9)
  segments(b[length(b)],(cumsum(statestatus[plotorder,maxtime])[2]- statestatus[plotorder[2],maxtime]/2),
           b[length(b)]+3,(cumsum(statestatus[plotorder,maxtime])[2]- statestatus[plotorder[2],maxtime]/2), xpd=NA)
          }, 
          zl, c("Current practice", "Novel regimen for RIF-R", "Novel RIF-R regimen + expanded RIF DSTS", "Novel regimen for all"))
mtext("Months elapsed since TB onset", 1, outer=T, cex=1, adj = .35)
mtext("% of RIF-R TB cohort", 2, outer=T, cex=1, line=-0.5)
# mtext("Status of RIF-R cohort over time", side=3, font=2, outer=T, line=0.5)
plot.new()
legend("center", legend = rev(fullstatenames[plotorder]), fill=rev(colors[1:8]), cex=1, xpd=NA,title = expression(bold("Status")))
dev.off()

# can also look at above differences for other subsets e.g. +- moxi

# could change smear to 0 when treated (when on any treatment) if want it to reflect reduced transmission during that time. 
# As currently coded, need to combine TB status and (original, for the current or pendinng episode) smear status to judge current infectiousness. 


# tally up time alive: # need to cut all scenarios off at the same point, say 3 years - no, need a longer time window, see below:

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


# time deceased i.e. YLL (compare to time alive, should see inverse)
# but don't trust this because of tail of life expectancy
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=8, cutofftime=12*followyears, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,mean))
lapply(l, function(x) apply(x[,1,]- l$baseline[,1,],2,mean)/12)


# cures, as opposed to death or still on treatment at xx months:
t <- 48; 
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
## it's possible the problem was different maxtimes between the scenarios, since I was setting mortality limit for each internally based on the max time recorded, 
## meaning longer maxtime and more opportunity for natural mortality in the one with longer time steps -- but that should cause more mortality for baseline. 
# with larger samples it's clear the mortality reduction is small but ~2%. which would make sense of >50% of mortality is before 1st treatment, another 50%+ is after improper treatment, and after MDR treatment the new regimen reduces relapse/failure from 10% to 2% (10% reduction in recurrence prevents death in 5% of this 25%?)
# 5x greater impact for novelrr. will want to do sensitivity analysis around xpert coverage params, make sure this responds as expected. 

# time to cure, if cured (where cure is assumed to happen when treatment stops):
par(mar=c(3,3,3,1))
l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
                   simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply(x[,1,],2,summary))
boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRR-TB patients only")
l <-  loutcomeboot(individualoutcomefunction =time.and.courses.to.cure, 
                   simoutput=impact, c=c, include=(c$RIF==0), copies = 3, desiredsize = 1e5)
boxplot(lapply(l, function(x) x[,1,1]), main="Time (in months) from TB onset to\nsuccessful treatment completion, for those ultimately cured,\nRS-TB patients only")


# and for dsT:
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, 
                  simoutput=dst, c=c, include=(c$MOXI==1), copies = 3, desiredsize = 1e4)
lapply(l, function(x) apply(x[,1,],2,summary))
lapply(l, function(x) apply(x[,1,],2,mean)) 
boxplot(lapply(l, function(x) x[,1,1]), main="Total infectious time, RR-TB cases")


# estimated reduction in infectious time and force of infection (not accounting for infectiousness during treatment, changes in a person's infectiousness over time, etc): 
# RIF-R
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, simoutput=impact, c=c, include=(c$RIF==1), copies = 3, desiredsize = 1e5)
lapply(l, function(x) apply((l$baseline[,1,]-x[,1,]), 2, sum)/apply(l$baseline[,1,],2,sum))
#(going 8-38-50%)
# and consider RIF-S:
l <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, 
                  simoutput=impact, c=c, include=(c$RIF==0), copies = 3, desiredsize = 1e4)
lapply(l, function(x) apply((l$baseline[,1,]-x[,1,]), 2, sum)/apply(l$baseline[,1,],2,sum))
#(going 0-0-3%)
# and for all TB, going ~0 to 2o to 6%
## ** why is novelrr increasing total TB infectious time (increasing infectious time among RS's)? 
## let's rerun and check this too - done, and there's no diff (within chance) for novelrr and novelrrx (as should be) and a 2-5% reduction in infectious time for pannovel

reduction <- lapply(impact, function(x)
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
# updated version of below:
lm <- loutcomeboot(individualoutcomefunction = time.in.state, simoutput = dst, c = c, include = c$MOXI==1, copies = 5, desiredsize = 1e4, states=c(1,2,5), characteristics=c("RIF"), characteristicvals=c(1))
lr <- loutcomeboot(individualoutcomefunction = time.in.state, simoutput = dst, c = c, include = c$RIF==1, copies = 5, desiredsize = 1e4, states=c(1,2,5), characteristics=c("RIF"), characteristicvals=c(1))
lb <- loutcomeboot(individualoutcomefunction = time.in.state, simoutput = dst, c = c, include = c$MOXI==1, copies = 5, desiredsize = 1e4, states=c(1,2,5), characteristics=c("BDQ"), characteristicvals=c(1))
lp <- loutcomeboot(individualoutcomefunction = time.in.state, simoutput = dst, c = c, include = c$MOXI==1, copies = 5, desiredsize = 1e4, states=c(1,2,5), characteristics=c("PA"), characteristicvals=c(1))

par(mfrow=c(3,2), oma=c(1,2,3,1))
layout(array(c(1,1,3,2,2,4), dim=c(3,2)))
boxplot(lapply(lr, function(z) apply(z[,1,], 2, mean)), main="Rifampin", col=colors, ylim=c(0,1))
boxplot(lapply(lm, function(z) apply(z[,1,], 2, mean)), main="Moxifloxacin", col=colors, ylim=c(0,1))
boxplot(lapply(lb, function(z) apply(z[,1,], 2, mean)), main="Bedaquiline", col=colors, ylim=c(0,1))
boxplot(lapply(lp, function(z) apply(z[,1,], 2, mean)), main="Pretomanid", col=colors, ylim=c(0,1))
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



# showing with and without regimen delays:
# should look at the same population I'm considering in the main aalysis. Probably all patietns or RIF==1.
# I initially consider proportion cured in a given round, but total infectious time will be more important
delays <- readRDS(file = paste0("delays.",date,".RDS"))
# d <- lapply(loutcomeboot(step4812cure, delays, c, include=(rep(1, length=nrow(c))), copies=50, desiredsize = 1e5), function(y) apply(y[,1,], 2, mean))
# l <- lapply(loutcomeboot(step4812cure, impact, c, include=(rep(1, length=nrow(c))), copies=50, desiredsize = 1e5), function(y) apply(y[,1,], 2, mean))
dinf <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, 
                     simoutput=delays, c=c, include=(c$RIF==1), copies = 50, desiredsize = 1e5)
linf <- loutcomeboot(individualoutcomefunction =time.in.state, states=c(1,2,6), cutofftime=12*10, carryforward=T, 
                     simoutput=impact, c=c, include=(c$RIF==1), copies = 50, desiredsize = 1e5)
di <- lapply(dinf, function(x) apply(x[,1,], 2, mean))
li <- lapply(linf, function(x) apply(x[,1,], 2, mean))
str(di)
ddf <- data.frame(di); ldf <- data.frame(li)
ddf$params <- "Delay for alternate regimens"
ldf$params <- "No delays"
require(dplyr)
require(data.table)
m <- bind_rows(ldf, ddf)
mm <- melt(m, id.vars="params")
colnames(mm) <- c("Parameters", "Scenario", "Inftime")
mm$Parameters <- factor(mm$Parameters, levels=c("No delays", "Delay for alternate regimens"))
require(ggplot2)
ggplot(data = mm, aes(x=Scenario, y=Inftime)) +
  geom_boxplot(aes(fill=Parameters)) +
  ggtitle("Effect of regimen delays") + theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Average time with untreated active TB (months)") + ylim(0,30)
# rm(delays)
