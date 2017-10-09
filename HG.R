setwd("D:/R Working Directory/Others statistics/Hakan G")

data<-read.csv("Hakan data.csv", sep=";", dec=",", header=T)
head(data, 3)
str(data)

data<- within(data, Condition<- relevel(Condition, ref = "Fresh"))

##Libraries
library(ggplot2)
library(nlme)
library(xlsx)
library(Hmisc)

t<-poly(unique(data$time), 3)
data[, paste("ot", seq_along(1:3), sep="")]<-t[data$time, seq_along(1:3)] ## did not work!!!!!


## Built-up of growth curves
##-----------------------------------For pmai----------------------------------------------------------
basic<-gls(data=data, pmai~1, method="ML", na.action=na.exclude) ##Built a basic model

randominter<-lme(data=data, pmai~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=100))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID, control=list(maxIter=1000, opt="optim"))
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition,control=list(msMaxIter=100, opt="optim"))
m7<-update(m6, .~poly(time, 3)*Condition)


##Check with ANOVA up to which polynomial there is a significant effect of time
a<-anova(basic, randominter,
      time1, time2, time3,
      randomslope1, randomslope2,randomslope3,
      m4, m5, m6, m7) 
a 
write.xlsx(a[ , 2:9], "pmai.xlsx") ## create and Excel output

##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
randomslope1<-update(time1, .~., random=~time|sampleID,
                     control=list(maxIter=1000, opt="optim", msMaxIter=100))
randomslope2<-update(time2, .~., random=~poly(time, 2)|sampleID,
                     control=list(maxIter=1000, opt="optim", msMaxIter=100))
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2,
      randomslope1, randomslope2,
      m4, m5, m6)

summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "pmai_t.xlsx")

plot(m6)
qqnorm(m6)


##Plot model
graph<-ggplot(data=data, aes(time, pmai, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Plasma membrane- and acrosome-intact sperm (%)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For as----------------------------------------------------------
basic<-gls(data=data, as~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, as~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID, control=list(maxIter=1000, opt="optim"))
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
      time1, time2, time3,
      randomslope1, randomslope2,randomslope3,
      m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "as.xlsx") ## create and Excel output

##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
randomslope2<-update(time2, .~., random=~poly(time, 2)|sampleID,
                     control=list(maxIter=1000, opt="optim", msMaxIter=100))
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2,
      randomslope1, randomslope2,
      m4, m5, m6)

summary(m6)
intervals(m6)

t<-summary(m6)$tTable
t
write.xlsx(t, "as_t.xlsx")

plot(m6)
qqnorm(m6)

##Plot model
graph<-ggplot(data=data, aes(time, as, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Acrosome-defected sperm (%)")+
  xlab("Time post thaw (hours)")

##-----------------------------------For hmmp----------------------------------------------------------
basic<-gls(data=data, hmmp~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, hmmp~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=100, opt="optim", msMaxIter=100))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))

randomslope1<-update(time1, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
##When I tried to add a random effect on time^3, an error occured.
##So this model goes up to polynomial order 2.


a<-anova(basic, randominter,
         time1, time2,
         randomslope1, randomslope2,
         m4, m5, m6)
a
write.xlsx(a[ , 2:9], "hmmp.xlsx")

summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "hmmp_t.xlsx")

plot(m6)
qqnorm(m6)

##Plot model
graph<-ggplot(data=data, aes(time, hmmp, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Sperm with high mitochondrial membrane potential (%)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For at----------------------------------------------------------
basic<-gls(data=data, at~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, at~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=100, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2,randomslope3,
         m4, m5, m6, m7)

a
write.xlsx(a[ , 2:9], "at.xlsx")

summary(m7)
intervals(m7)

t<-summary(m7)$tTable
t
write.xlsx(t, "at_t.xlsx")

plot(m7)
qqnorm(m7)

##Plot model
graph<-ggplot(data=data, aes(time, at, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m7), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Mean DFI")+
  xlab("Time post thaw (hours)")

##-----------------------------------For sd----------------------------------------------------------
basic<-gls(data=data, sd~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, sd~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2,randomslope3,
         m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "sd.xlsx")

##Remove the excessive random effect of bull on the non-significant polynomial (of order 2 and 3),
##so as to have valid intervals calculation
time1<-update(randominter, .~.+time, correlation=corCAR1())
randomslope1<-update(time1, .~., random=~time|sampleID)
m4<-update(randomslope1, .~. + Condition)
m5<-update(m4, .~time*Condition)

a<-anova(basic, randominter,
      time1, randomslope1,m4, m5)
a

summary(m5)
intervals(m5)

t<-summary(m5)$tTable
t
write.xlsx(t, "sd_t.xlsx")

plot(m5)
qqnorm(m5)

##Plot model
graph<-ggplot(data=data, aes(time, sd, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m5), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("SD of DFI")+
  xlab("Time post thaw (hours)")

##-----------------------------------For dfi----------------------------------------------------------
basic<-gls(data=data, dfi~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, dfi~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2, randomslope3,
         m4, m5, m6, m7)
a ## m6 fitted best
write.xlsx(a[ , 2:9], "dfi.xlsx")


##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
randomslope1<-update(time1, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2, 
      randomslope1, randomslope2,
      m4, m5, m6)

summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "dfi_t.xlsx")

plot(m6)
qqnorm(m6)

##Plot model
graph<-ggplot(data=data, aes(time, dfi, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("DNA fragmentation index (%)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For no----------------------------------------------------------
basic<-gls(data=data, no~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, no~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=100, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2, randomslope3,
         m4, m5, m6, m7)
a

write.xlsx(a[ , 2:9], "no.xlsx")

summary(m7)
intervals(m7) ##not working!!!

##COefficient of the 2nd and 3rd polynomial time term in m7 is not significant,
##so I try to reduce it down to a 1st order polynomial
randomslope1<-update(time1, .~., random=~time|sampleID)
m4<-update(randomslope1, .~. + Condition)
m5<-update(m4, .~time*Condition)
anova(basic, randominter,
      time1, 
      randomslope1,
      m4, m5)
summary(m5)
intervals(m5)

t<-round(summary(m5)$tTable, 2)
t
write.xlsx(t, "no_t.xlsx")

plot(m5)
qqnorm(m5)


##Plot model
graph<-ggplot(data=data, aes(time, no, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m5), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Nitride Oxide (channels)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For dcfh----------------------------------------------------------
basic<-gls(data=data, dcfh~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, dcfh~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=100, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2, randomslope3,
         m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "dcfh.xlsx")


summary(m7)
intervals(m7) ##not working!!!

##COefficient of the 3rd polynomial time term in m7 is not significant,
##so I try to reduce it down to a 2nd order polynomial
randomslope1<-update(time1, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2, 
      randomslope1, randomslope2,
      m4, m5, m6)
summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "dcfh_t.xlsx")

plot(m6)
qqnorm(m6)

##Plot model
graph<-ggplot(data=data, aes(time, dcfh, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("DCFH (channels)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For dhr----------------------------------------------------------
basic<-gls(data=data, dhr~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, dhr~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope1, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2,randomslope3,
         m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "dhr.xlsx")

##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
randomslope1<-update(time2, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2, 
      randomslope1, randomslope2,
      m4, m5, m6)

##Try removing the 2nd order polynomial too
randomslope1<-update(time1, .~., random=~time|sampleID)
m4<-update(randomslope1, .~. + Condition)
m5<-update(m4, .~time*Condition)
anova(basic, randominter,
      time1, 
      randomslope1,
      m4, m5)
summary(m5)

summary(m5)
intervals(m5)

t<-round(summary(m5)$tTable, 2)
t
write.xlsx(t, "dhr_t.xlsx")

plot(m5)
qqnorm(m5)

##Plot model
graph<-ggplot(data=data, aes(time, dhr, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m5), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("DHR (channels)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For mitosox ----------------------------------------------------------
basic<-gls(data=data, mitosox~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, mitosox~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2, randomslope3,
         m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "mitosox.xlsx")

##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
randomslope1<-update(time1, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
anova(basic, randominter,
      time1, time2, 
      randomslope1, randomslope2,
      m4, m5, m6)


summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "mitosox_t.xlsx")

plot(m6)
qqnorm(m6)

##Plot model
graph<-ggplot(data=data, aes(time, mitosox, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("MitoSOX (channels)")+
  xlab("Time of incubation (hours)")

##-----------------------------------For progr----------------------------------------------------------
basic<-gls(data=data, progr~1, method="ML", na.action=na.exclude)

randominter<-lme(data=data, progr~1, random=~1|sampleID, method="ML",
                 na.action=na.exclude, control=list(maxIter=1000, msMaxIter=100, opt="optim"))

time1<-update(randominter, .~.+time, correlation=corCAR1())
time2<-update(time1, .~ poly(time, 2))
time3<-update(time2, .~ poly(time, 3))
randomslope1<-update(time3, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
randomslope3<-update(randomslope2, .~., random=~poly(time, 3)|sampleID)
m4<-update(randomslope3, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)
m7<-update(m6, .~poly(time, 3)*Condition)


a<-anova(basic, randominter,
         time1, time2, time3,
         randomslope1, randomslope2, randomslope3,
         m4, m5, m6, m7)
a
write.xlsx(a[ , 2:9], "progr.xlsx")


##Remove the excessive random effect of bull on the non-significant polynomial (order 3),
##so as to have valid intervals calculation
time2<-update(time1, .~ poly(time, 2))
randomslope1<-update(time1, .~., random=~time|sampleID)
randomslope2<-update(randomslope1, .~., random=~poly(time, 2)|sampleID)
m4<-update(randomslope2, .~. + Condition)
m5<-update(m4, .~time*Condition)
m6<-update(m5, .~poly(time, 2)*Condition)

anova(basic, randominter,
      time1, time2,
      randomslope1, randomslope2,
      m4, m5, m6)

summary(m6)
intervals(m6)

t<-round(summary(m6)$tTable, 2)
t
write.xlsx(t, "progr_t.xlsx")

graph<-ggplot(data=data, aes(time, progr, shape=Condition))
graph+stat_summary(fun.y=mean, geom="point", size=3)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2)+
  stat_summary(aes(y=fitted(m6), linetype=Condition), fun.y=mean, geom="line")+
  scale_x_continuous(breaks=c(0, 3, 6, 12, 24))+
  ylab("Progressively motile sperm (%)")+
  xlab("Time of incubation (hours)")



##____________________________________________________________________________________##

##______________________ Correlations _______________________________________________

dataFresh<-data[data$Condition=="Fresh", 6:16]
dataCryo<-data[data$Condition=="Cryopreserved", 6:16]


dataFresh<-as.matrix(dataFresh)
dataCryo<-as.matrix(dataCryo)

CorFresh<-rcorr(dataFresh, type="spearman")
CorCryo<-rcorr(dataCryo, type="spearman")

write.xlsx(CorFresh$r, "Corr Coef for Fresh.xlsx")
write.xlsx(CorCryo$r, "Corr Coef for Cryopreserved.xlsx")

write.xlsx(CorFresh$P, "Corr P for Fresh.xlsx")
write.xlsx(CorCryo$P, "Corr P for Cryopreserved.xlsx")



##Scatterplots

graph<-ggplot(data=data, aes(no, at, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("NO")+ylab("Mean DFI")

graph<-ggplot(data=data, aes(dcfh, at, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DCFH")+ylab("Mean DFI")

graph<-ggplot(data=data, aes(dhr, at, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DHR")+ylab("Mean DFI")

graph<-ggplot(data=data, aes(mitosox, at, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("MitoSOX")+ylab("Mean DFI")

graph<-ggplot(data=data, aes(no, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("NO")+ylab("DFI%")

graph<-ggplot(data=data, aes(dcfh, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DCFH")+ylab("DFI% (%)")

graph<-ggplot(data=data, aes(dhr, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DHR")+ylab("DFI% (%)")

graph<-ggplot(data=data, aes(mitosox, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("MitoSOX")+ylab("DFI% (%)")

graph<-ggplot(data=data, aes(no, sd, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("NO")+ylab("SD of DFI")

graph<-ggplot(data=data, aes(dcfh, sd, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DCFH")+ylab("SD of DFI")

graph<-ggplot(data=data, aes(dhr, sd, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("DHR")+ylab("SD of DFI")

graph<-ggplot(data=data, aes(mitosox, sd, colour=Condition))
graph+geom_point()+
  geom_smooth(method="lm", se=F, size=1)+
  xlab("MitoSOX")+ylab("SD of DFI")


##___________________________________________________________________________________________


## Correlations between SCSA and ROS parameters
graph<-ggplot(data=data, aes(dcfh, at, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("Mean DFI")+xlab("DCFH positive sperm(%)")

graph<-ggplot(data=data, aes(dhr, at, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("Mean DFI")+xlab("DHR positive sperm(%)")

graph<-ggplot(data=data, aes(no, at, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("Mean DFI")+xlab("NO")

graph<-ggplot(data=data, aes(mitosox, at, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("Mean DFI")+xlab("Mito-SOX levels")

graph<-ggplot(data=data, aes(dcfh, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("DNA fragmentation index (%)")+xlab("DCFH positive sperm(%)")

graph<-ggplot(data=data, aes(dhr, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("DNA fragmentation index (%)")+xlab("DHR")

graph<-ggplot(data=data, aes(no, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("DNA fragmentation index (%)")+xlab("NO")

graph<-ggplot(data=data, aes(mitosox, dfi, colour=Condition))
graph+geom_point()+
  geom_smooth(se=F, method="lm", colour="black")+
  facet_grid(time~Condition, scales="free")+
  ylab("DNA fragmentation index (%)")+xlab("Mito-SOX")
