rm(list = ls(all = TRUE))


                ####   INDEX TABLE
### load required packages 
### 1. Transport of root aphid eggs and ant larvae by Lasius flavus and other ants ###
### 2. The nature of the cue triggering the transport of root aphid eggs by Lasius flavus ###
### 3. Protection of root aphids eggs by Lasius flavus against pathogens ###
### 4. Protection of root aphids eggs by Lasius flavus against predators ###
### 5. Hatching of root aphid eggs ###



### load required packages
library(coxme)
library(survival)
library(ggplot2)
library(dplyr)
library(emmeans)
library(lme4)
library(RVAideMemoire)
library(effects)
library(car)
library(tidyr)
library(survminer)
library(multcomp)
library(riskRegression)
library(adjustedCurves)
library(pammtools)



######################################################################################
### 1. Transport of root aphid eggs and ant larvae by Lasius flavus and other ants ###
######################################################################################

datasetpref<-read.table("retrieval.txt",header=TRUE)  
table(datasetpref$colony,datasetpref$trophobiont)

# 1.a. Transport colony LF_B
####################################
datasetpref_LF_B<-datasetpref[datasetpref$colony=="LF_B",]
fit_LF_B<-coxme(Surv(time, status) ~ trophobiont+(1|set),datasetpref_LF_B)
anova(fit_LF_B)
ph.test <- cox.zph(fit_LF_B) #proportional hazards assumption ok

#plot, accounting for random variable set
cox_mod_B <- coxph(Surv(time, status) ~ trophobiont +frailty(set),
                   data=datasetpref_LF_B, x=TRUE)
predict_fun <- function(...) {
  predictRisk(...)
}
datasetpref_LF_B$trophobiont<-as.factor(datasetpref_LF_B$trophobiont)
adjsurv_B <- adjustedsurv(data=datasetpref_LF_B,
                          variable="trophobiont",
                          ev_time="time",
                          event="status",
                          method="direct",
                          bootstrap=TRUE,
                          n_boot=500,
                          outcome_model=cox_mod_B,
                          predict_fun=predict_fun)
plotLF_B<-plot(adjsurv_B, conf_int=T, use_boot=TRUE,custom_colors=c("black","gray80"))+theme(legend.position="top")+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))

#survival plot, not accounting for random variable set
fit_LF_B_s <- survfit( Surv(time, status) ~ trophobiont, data = datasetpref_LF_B)
ggsurvplot_facet(fit_LF_B_s, datasetpref_LF_B, facet.by = "colony",fun="event", palette = c("black","gray80","blue","red"),pval = F,conf.int=T)+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))

# 1.b. Transport colony LF_C
###################################
datasetpref<-read.table("retrieval.txt",header=TRUE)  
datasetpref_LF_C<-datasetpref[datasetpref$colony=="LF_C",]

fit_LF_C<-coxme(Surv(time, status) ~ trophobiont+ (1|set), datasetpref_LF_C)
anova(fit_LF_C)
ph.test <- cox.zph(fit_LF_C) #proportional hazards assumption ok 

#plot, accounting for random variable set
cox_mod_C <- coxph(Surv(time, status) ~ trophobiont +frailty(set),
                   data=datasetpref_LF_C, x=TRUE)
predict_fun <- function(...) {
  predictRisk(...)
}
datasetpref_LF_C$trophobiont<-as.factor(datasetpref_LF_C$trophobiont)
adjsurv_C <- adjustedsurv(data=datasetpref_LF_C,
                          variable="trophobiont",
                          ev_time="time",
                          event="status",
                          method="direct",
                          bootstrap=TRUE,
                          n_boot=500,
                          outcome_model=cox_mod_C,
                          predict_fun=predict_fun)
plotLF_C<-plot(adjsurv_C, conf_int=T, use_boot=TRUE,custom_colors = c("black","gray80"))+theme(legend.position="top")+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))

#survival plot, not accounting for random variable set
fit_LF_C_s <- survfit( Surv(time, status) ~ trophobiont, data = datasetpref_LF_C )
ggsurvplot_facet(fit_LF_C_s, datasetpref_LF_C, facet.by = "colony",fun="event", palette = c("black","gray80","blue","red"),pval = F,conf.int=T)+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))


# 1.c. Transport colony LF_A
#################################
datasetpref<-read.table("retrieval.txt",header=TRUE) 
datasetpref_LF_A<-datasetpref[datasetpref$colony=="LF_A",] #select only colony LF_A

#compare aphid_egg vs ant_larvae
datasetpref_eggsvslarvae <- datasetpref_LF_A[datasetpref_LF_A$trophobiont %in% c("aphid_egg", "ant_larva"), ]
fit_LF_A<-coxme(Surv(time, status) ~ trophobiont+ (1|set),datasetpref_eggsvslarvae)
anova(fit_LF_A)
summary(fit_LF_A)
ph.test <- cox.zph(fit_LF_A) #proportional hazards assumption violated

datasetpref_eggsvslarvae$status[datasetpref_eggsvslarvae$time>15] <- 0 #to prepare for binomial analysis: if time >15 min than status 0
datasetpref_eggsvslarvae$time[datasetpref_eggsvslarvae$time>15]<-15 #to prepare for binomial analysis: if time >15 min set time 15
fitbinom_eggsvslarvae<-glmer((status)~trophobiont+(1|set),family="binomial",datasetpref_eggsvslarvae) #constant hazard violated --> binom model after particular time
overdisp.glmer(fitbinom_eggsvslarvae)
plot(allEffects(fitbinom_eggsvslarvae),type="response")
Anova(fitbinom_eggsvslarvae)


#############################################################################################
### 2. The nature of the cue triggering the transport of root aphid eggs by Lasius flavus ###
#############################################################################################

datasetpref<-read.table("retrieval.txt",header=TRUE) 
datasetpref_LF_A<-datasetpref[datasetpref$colony=="LF_A",] #select only colony LF_A

#compare aphid_egg vs aphid_egg_hexane_treated 
datasetpref_eggsvseggshx <- datasetpref_LF_A[datasetpref_LF_A$trophobiont %in% c("aphid_egg", "egg_hexane_treated"), ]
fit_eggsvseggshx<-coxme(Surv(time, status) ~ trophobiont+ (1|set),datasetpref_eggsvseggshx)
anova(fit_eggsvseggshx)
summary(fit_eggsvseggshx)
ph.test <- cox.zph(fit_eggsvseggshx) #proportional hazards assumption violated

datasetpref_eggsvseggshx$status[datasetpref_eggsvseggshx$time>15] <- 0 #to prepare for binomial analysis: if time >15 min than status 0
datasetpref_eggsvseggshx$time[datasetpref_eggsvseggshx$time>15]<-15 #to prepare for binomial analysis: if time >15 min set time 15
fitbinom_eggsvseggshx<-glmer((status)~trophobiont+(1|set),family="binomial",datasetpref_eggsvseggshx) #constant hazard violated --> binom model after particular time
overdisp.glmer(fitbinom_eggsvseggshx)
plot(allEffects(fitbinom_eggsvseggshx),type="response")
Anova(fitbinom_eggsvseggshx)

#survival plot, accounting for random variable set
datasetpref<-read.table("retrieval.txt",header=TRUE) 
datasetpref_LF_A<-datasetpref[datasetpref$colony=="LF_A",] #select only colony LF_A
datasetpref_LF_A<-datasetpref_LF_A[datasetpref_LF_A$trophobiont!="egg_hexane_methanol_treated",]
cox_mod_A <- coxph(Surv(time, status) ~ trophobiont +frailty(set),
                   data=datasetpref_LF_A, x=TRUE)
predict_fun <- function(...) {
  predictRisk(...)
}
datasetpref_LF_A$trophobiont<-as.factor(datasetpref_LF_A$trophobiont)

adjsurv_A <- adjustedsurv(data=datasetpref_LF_A,
                          variable="trophobiont",
                          ev_time="time",
                          event="status",
                          method="direct",
                          bootstrap=TRUE,
                          n_boot=500,
                          outcome_model=cox_mod_A,
                          predict_fun=predict_fun)
plotLF_A<-plot(adjsurv_A, conf_int=T, use_boot=TRUE,custom_colors=c("black","gray80","blue"))+theme(legend.position="top")+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))

# survivalplot without considering random term
fit_LF_A_s <- survfit(Surv(time, status) ~ trophobiont, data = datasetpref_LF_A)
ggsurvplot_facet(fit_LF_A_s, datasetpref_LF_A, facet.by = "colony",fun="event", palette = c("black","gray80","blue","red"),pval = F,conf.int=T)+ scale_x_continuous(breaks=seq(0,120,30))+scale_y_continuous(expand = c(0,0),limits = c(0, 1))+scale_x_continuous(expand = c(0, 0))



############################################################################
### 3. Protection of root aphids eggs by Lasius flavus against pathogens ###
############################################################################

datasetfungus<-read.table("fungus.txt",header=TRUE)  

    #parametric test: assumptions violated 
    #fit<-glm(cbind(end,5-end)~treatment+colony,binomial,datasetfungus)
    # Anova(fit)
    # library(DHARMa)
    # simulationOutput <- simulateResiduals(fittedModel = fit, plot = T)

testfungusA<-wilcox.test(end~treatment,datasetfungus[datasetfungus$colony=="A",]) #Wilcoxon rank sum test for colony A
testfungusB<-wilcox.test(end~treatment,datasetfungus[datasetfungus$colony=="B",]) #Wilcoxon rank sum test for colony B



############################################################################
### 4. Protection of root aphids eggs by Lasius flavus against predators ###
############################################################################

datasetpred<-read.table("predation.txt",header=TRUE)  
testpred<-wilcox.test(5-end~treatment,datasetpred) #Wilcoxon rank sum test 

datasetpred %>% group_by(treatment) %>%
  summarize( mean = mean(5-end),sd= sd(5-end)) #mean eggs damaged min and plus beetle treatment



######################################
### 5. Hatching of root aphid eggs ###
######################################

datasethatch<-read.table("hatching.txt",header=TRUE)  
datasethatch$replicate<-as.factor(datasethatch$replicate)
datasethatch$treatment<-as.factor(datasethatch$treatment)

fit1<-coxme(Surv(time_after_start, hatched) ~ treatment+ (1|replicate), datasethatch)
anova(fit1)
ph.test <- cox.zph(fit1)  #ok


# 5.a. survivalplot adjusted for random term
############################################

  #https://stackoverflow.com/questions/77654561/plotting-adjusted-survival-curve-for-a-mixed-effects-cox-regression-and-or-time
cox_mod <- coxph(Surv(time_after_start, hatched) ~ treatment +frailty(replicate),
                 data=datasethatch, x=TRUE)
anova(cox_mod)
predict_fun <- function(...) {
  predictRisk(...)
}
datasethatch$treatment<-as.factor(datasethatch$treatment)
adjsurv <- adjustedsurv(data=datasethatch,
                        variable="treatment",
                        ev_time="time_after_start",
                        event="hatched",
                        method="direct",
                        bootstrap=TRUE,
                        n_boot=100,
                        outcome_model=cox_mod,
                        predict_fun=predict_fun)
plot(adjsurv, conf_int=T, use_boot=TRUE,xlim=c(0, 30))


# 5.b. survivalplot without considering random term
##################################################
datasethatch<-read.table("hatching.txt",header=TRUE)  
fit <- survfit( Surv(time_after_start, hatched) ~ treatment, data = datasethatch )
ggsurvplot(fit, datasethatch,fun="event", palette = c("black","gray80","red","blue"),pval = F,conf.int=F)
