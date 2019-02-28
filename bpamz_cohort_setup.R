# setwd("C:/Users/EAK/Google Drive/TPPs and population impacts/universal regimen/universal")
require(binom)
require(arm)

region <- "SEA"

allnotifs <- read.csv("../TB_notifications_2019-01-07.csv")
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
sum(notifs2017[,c("new_ep","new_labconf","new_clindx")], na.rm=T)/
  sum(notifs2017[,c("new_ep","new_labconf","new_clindx","ret_rel_labconf","ret_rel_clindx","ret_rel_ep","ret_nrel")], na.rm=T)


# Moxi and PZA conf ints:
binom.confint(344,6701)
binom.confint(136,473)
binom.confint(247,6925)
binom.confint(71,469)
binom.confint(18,6921)
binom.confint(11,463)
binom.confint(252,2978)
binom.confint(255,1673)

Back-calculating prorporion 2 mo culture positivity:
    All HRZE: invlogit((logit(.06278)-2.5289+2.5018*log(6))/.4399) = 15.4%
    BPaMZ: invlogit((logit(.06278)-2.5289+2.5018*log(4))/.4399) = 1.8%
    INH-S: invlogit((logit(.06278*.8896)-2.5289+2.5018*log(6))/.4399) = 12.1%
  INH-R: invlogit((logit(3*.06278*.8896)-2.5289+2.5018*log(6))/.4399) = 69.8%?

Back-calculating failure
approx(c(log(2), log(6)), c(0.154, 0.012), c(log(3), log(4), log(5)), method="linear")
curve <- lm(t~f)
curve
curve(log(3))
