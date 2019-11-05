# setwd("C:/Users/EAK/Google Drive/TPPs and population impacts/universal regimen/bpamz")
require(binom)
require(arm)


allnotifs <- read.csv("TB_notifications_2019-01-07.csv")

region <- "SEA"
notifs <- subset(allnotifs, g_whoregion==region)
notifs2017 <- subset(notifs, year==2017)
notifs2016 <- subset(notifs, year==2016)
sum(notifs2017[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2017[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)
# summary(rowSums(notifs2017[,c("new_ep","new_labconf","new_clindx")])/
#   rowSums(notifs2017[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")]))
sum(notifs2016[,c("new_ep","new_labconf","new_clindx")])/
  sum(notifs2016[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")])
sum(notifs2016[notifs2016$country=="India",c("new_ep","new_labconf","new_clindx")])/
  sum(notifs2016[notifs2016$country=="India",c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")])
sum(notifs2017[notifs2017$country=="India",c("new_ep","new_labconf","new_clindx")])/
  sum(notifs2017[notifs2017$country=="India",c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")])

notifs <- subset(allnotifs, country=="South Africa")
notifs2017 <- subset(notifs, year==2017)
notifs2016 <- subset(notifs, year==2016)
notifs2015 <- subset(notifs, year==2015)
notifs2014 <- subset(notifs, year==2014)
sum(notifs2017[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2017[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)
sum(notifs2016[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2016[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)
sum(notifs2015[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2015[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)
sum(notifs2014[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2014[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)

# SEA moxi resistance:

#RR
h1 <- 2.0*1.8
h2 <- 2.0+(0.5*(13.2-2.0))
i <- (sqrt(h1*h2))/13.2*16.3
(i*307+sqrt(h1*h2)*156)/(307+156)

#RS
h1 <- 0.5*1.8
h2 <- 0.5+(0.5*(6.2-0.5))
(sqrt(h1*h2))
i <- (sqrt(h1*h2))/6.2*2.3
(i*4651+sqrt(h1*h2)*2274)/(4651+2274)

# South Africa
h1 <- 2.1*1.8
h2 <- 2.1+(0.5*(10.1-2.1))
sqrt(h1*h2)

h1 <- 0.15*1.8
h2 <- 0.15+(0.5*(0.3-0.15))
sqrt(h1*h2)

# Moxi and PZA conf ints:

binom.confint(344,6701)
binom.confint(136,473)
binom.confint(247,6925)
binom.confint(71,469)

binom.confint(69,6921)
binom.confint(11,463)
binom.confint(252,2978)
binom.confint(255,1673)

# SAF PZA confints:
binom.confint(round(0.013*838+0.013*657), 838+657)
binom.confint(round(0.39*39+0.49*34), 39+34)
# SAf moxi contints
#high
binom.confint(round(0.003*908+0*620), 908+620)
binom.confint(round(0.038*40+0*31), 40+31)
#low
binom.confint(round(0.005*910+0.003*621), 910+621)
binom.confint(round(0.084*41+0.122*33), 41+33)

sqrt(5.6*6.2)
0.002*1.8
0.002+(0.0039-0.002)/2
sqrt(.0036*.00295)

# MOXI in SAf, editing original supplement written for SEA:
Prevalence of MOXI resistance in SAf is also based on pooling data from two sites (Zignol 2016).
Because MOXI critical concentrations and clinical breakpoints have been lowered since the MGIT cutoff were selected for those drug resistance surveys,
we used other data on the distribution of MOXI MICs among clinical MDR-TB isolates (Rigouts JAC 2016) 
to estimate what fraction of isolates in the South African drug resistance surveys would have been resistant at these lower breakpoints. 
To account for MIC differences between MICs in MGIT media (used in the drug resistance surveys) versus on LJ media (used by Rigouts et al), 
we translated all LJ MICs down one dilution (e.g. treating an MIC of 0.5 on LJ as equivalent to a 0.25 on MGIT).                                                                                                                                            

Of MDR-TB isolates with MIC of 4 ug/ml or above on LJ (2 or above on LJ), there were 1.8 times more with an MIC of 2 or above on LJ. 
Conversely, of the isolates with an MIC of 1 or 2 ug/ml on LJ (i.e., resistant at 0.5 but not at 2.0 ul/ml on MGIT), 
50% had an MIC of 2 rather than 1, and were taken to reflect MOXI resistance above the clinical breakpoint. 

Applying the former ratio to the pooled prevalence of 2.8% for MOXI-R at a MGIT MIC cutoff of 2.0 ug/ml among RIF-R TB in South Africa, 
we estimated a prevalence of 2.8%*1.8 = 5.6%. Conversely, applying the latter 50% proportion to the pooled prevalence of 9.5% 
  for MOXI-R among RIF-R at a cutoff of 0.5 ug/ml in South Africa, 
we estimated prevalence of 2.8%+1/2*(9.5%-2.8%) = 6.2% MOXI-R at a cutoff of 1.0 ug/ml on MGIT. 
We took a geometric mean of these for a final estimate of 5.9% prevalence of high-level MOXI-R among RIF-R. 
Because of limited data on the prevalence of MOXI-R at MIC cutoffs at of 0.25 versus 0.5, we used the MOXI-R prevalence at cutoff of 0.5 ug/ml,
9.5%, as a conservative estimate of the prevalence of total MOXI-R with MIC at or above the MGIT critical concentration of 0.25 ug/ml. 
This is likely an underestimate, but not by a large amount since the reported prevalence of levofloxacin and ofloxacin resistance were similar, respectively. 
The same process led to estimates for MOXI-R among RIF-S TB in South Africa of 0.3% high-level resistance (above the clinical breakpoint at which moxifloxacin could still add some anti-TB activity) and 0.4% total resistance. 
Ranges used in sensitivity analysis are the broadest of (a) a binomial confidence interval based on number of resistant isolates in the drug resistance surveys, (b) the range in point estimates over the sites surveyed, (c) the range between estimation methods. 


write.csv(paste0(allparams[,6],"-",allparams[,7]), file = "tables1 ranges.csv")


# Back-calculating prorporion 2 mo culture positivity:
#     All HRZE: invlogit((logit(.06278)-2.5289+2.5018*log(6))/.4399) = 15.4%
#     BPaMZ: invlogit((logit(.06278)-2.5289+2.5018*log(4))/.4399) = 1.8%
#     INH-S: invlogit((logit(.06278*.8896)-2.5289+2.5018*log(6))/.4399) = 12.1%
#   INH-R: invlogit((logit(3*.06278*.8896)-2.5289+2.5018*log(6))/.4399) = 69.8%?
# 
# Back-calculating failure
# approx(c(log(2), log(6)), c(0.154, 0.012), c(log(3), log(4), log(5)), method="linear")
# curve <- lm(t~f)
# curve
# curve(log(3))
