library(survival)
library(survminer)
library(tableone)
library(psych)
library(dplyr)

WHS <- read.csv(file="whas500.csv")

WHS$GENDER<-as.factor(WHS$GENDER)
WHS$CVD <- as.factor(WHS$CVD)
WHS$AFB<-as.factor(WHS$AFB)
WHS$SHO<-as.factor(WHS$SHO)
WHS$CHF<-as.factor(WHS$CHF)
WHS$AV3<-as.factor(WHS$AV3)
WHS$MIORD<-as.factor(WHS$MIORD)
WHS$MITYPE<-as.factor(WHS$MITYPE)
WHS$YEAR <- as.character(WHS$YEAR)

head(WHS)
View(WHS)
# Fitting Kaplan Meier Estimate for Males 
attach(WHS)
survfit.male<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS, subset = {GENDER==0})
summary(survfit.male)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.male, main = "Kaplan-Meier Survival Curve for Males", xlab = "Follow-up-Time", ylab = "Survival Probability", col = "red")
quantile(survfit.male)

# Fitting Kaplan Meier Estimate for Males on a log transformation
survfit.male.log<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS, conf.type = "log", subset = {GENDER==0})
summary(survfit.male.log)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.male.log, main = "Kaplan-Meier Survival Curve for Males\n on a log transformation", xlab = "Follow-up-Time", ylab = "Survival Probability", col = "red")
quantile(survfit.male.log)


# Fitting Kaplan Meier Estimate for Males on a Complementary log-log transformation
survfit.male.loglog<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS, conf.type = "log-log", subset = {GENDER==0})
summary(survfit.male.loglog)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.male.loglog, main = "Kaplan-Meier Survival Curve for Males\n on a log log transformation", xlab = "Follow-up-Time", ylab = "Survival Probability", col = "red")
quantile(survfit.male.loglog)

# Fitting Kaplan Meier Estimate for Females 
survfit.female<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS, subset = {GENDER==1})
summary(survfit.female)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.female, main = "Kaplan-Meier Survival Curve for Females", xlab = "Follow-up-Time", ylab = "Survival Probability",
     col = "blue")
quantile(survfit.female)

# Fitting Kaplan Meier Estimate for Females on a  log transformation  
survfit.female.log<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS,conf.type = "log", subset = {GENDER==1})
summary(survfit.female.log)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.female.log, main = "Kaplan-Meier Survival Curve for Females\n on a log scale", xlab = "Follow-up-Time",
     ylab = "Survival Probability", col = "blue")
quantile(survfit.female.log)

# Fitting Kaplan Meier Estimate for Females on a Complementary log log  transformation  
survfit.female.loglog<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS,conf.type = "log-log", subset = {GENDER==1})
summary(survfit.female.loglog)
#Plotting a Kaplan- Meier Survival Curve for Male 
plot(survfit.female.loglog, main = "Kaplan-Meier Survival Curve for Females\n on a log log scale",
     xlab = "Follow-up-Time", ylab = "Survival Probability", col = "blue")
quantile(survfit.female.loglog)

## From the Kaplan Meier Estimates, it appears that one of the confidence intervals is exactly one. The log transformation 
## could not take take of that. However one of the advantages of using the log-log transformation is that
### they provide upper confidence interval not exactly 1 but within a convincing limit. 


# Plotting a smooth survival curve for the data 
library(muhaz)
result.smooth<- muhaz(LENFOL,FSTAT)

#Plotting a hazard plot
plot(result.smooth$haz.est~result.smooth$est.grid, type = "l", xlab = "Time", ylab = "Hazard", main= "Hazard Plot")

result.km<- survfit(Surv(LENFOL,FSTAT)~1, data = WHS)
plot(result.km, main = "Smooth Survival Curve", xlab = "Time", ylab = "Survival Probability")
haz <- result.smooth$haz.est[-length(result.smooth$haz.est)]  
times <- result.smooth$est.grid
surv <- exp(-cumsum(haz*diff(times)))
lines(surv ~ times[-length(times)], lwd=2, col="blue")

# Assuming that the survival times follow an exponential distribution, find the maximum likelihood estimator lambda hat = d/V 
#Then plot the exponential survival curve.  

surv.km <- Surv(LENFOL,FSTAT)
number_of_death <- sum(FSTAT)
person_weeks <- sum(LENFOL)
maximum_likelihood_estimate <- number_of_death/person_weeks
exp_surv <- exp(-(maximum_likelihood_estimate*LENFOL))
plot(exp_surv, ylim = c(0,1), xlim = c(0,13),main = "Exponential Survival Curve",
     xlab = "Time (in weeks)", ylab = "Survival Probability", type = "l", col= "red")

## We can plot hazard plot for male and female as in Homework 2 question 5


# Nelson-Aalen Harrington-Fleming Survival Curve 

Fleming.Harrington<- survfit(Surv(LENFOL,FSTAT)~1, stype = 2, data = WHS)
summary(Fleming.Harrington)
plot(Fleming.Harrington, main = "Nelson Aaalen Harrington Fleming Survival Curve", xlab = "Time", 
     ylab = "Survival Probability")

#Fitting a Weibul Survival Distribution of the Fleming Harrington Estimates 

yy<- log(-log(Fleming.Harrington$surv))
xx<- log(Fleming.Harrington$time)
plot(yy~xx, pch= 20, xlab= "Log time", ylab = "Survival", main ="Weibull Distribution")
result.lm<- lm(yy~xx)
result.lm  
abline(result.lm, col = "red")  

#Estimating the parameters of a Weibull distribution that fits the Whas500 data
result.lm
#Call:
lm(formula = yy ~ xx)

#Coefficients:
#  (Intercept)           xx  
# -3.6584       0.4325  

# y = -3.6584 + 0.4325 
# slope = 0.4325    Intercept = ln(lambda) = -3.6584 Therefore lambda = 0.0257



# The weighted log rank test 
survdiff(Surv(LENFOL,FSTAT)~GENDER, rho = 0, data = WHS)

# Prentice modified Gehan test
survdiff(Surv(LENFOL,FSTAT)~GENDER, rho = 1, data = WHS)

#Log-rank test gives equal weights to every time points (i.e., wi = 1 for all i’s). The Gehan test assigns more
#weight on the early events.


wilcox.test(LENFOL~FSTAT, data = WHS)

#Wilcoxon test is more similar to the Prentice Modification of the Gehan test. They are both weighted tests,
#plus, Wilcoxon test is a special case of the Prentice Modification of the Gehan test.


# Fitting a cox proportional hazard model with History of Cardiovascular Disease (CVD) as a predictor
library(splines)
result.coxph.model<- coxph(Surv(LENFOL, FSTAT)~CVD, data = WHS)
summary(result.coxph.model)


# Constructing a 95% confidence interval for the hazard ration  using the parameter estimate and 
# its standard error. 

beta = result.coxph.model$coefficients
standard.error = 0.1680   #from the output above 
confidence.interval<- c(beta-1.96*standard.error, beta+1.96*standard.error)
confidence.interval




library(survival)
library(survminer)
library(tableone)
library(psych)
library(dplyr)


WHS <- read.csv(file="whas500.csv")

WHS$GENDER<-as.factor(WHS$GENDER)
WHS$CVD <- as.factor(WHS$CVD)
WHS$AFB<-as.factor(WHS$AFB)
WHS$SHO<-as.factor(WHS$SHO)
WHS$CHF<-as.factor(WHS$CHF)
WHS$AV3<-as.factor(WHS$AV3)
WHS$MIORD<-as.factor(WHS$MIORD)
WHS$MITYPE<-as.factor(WHS$MITYPE)
WHS$YEAR <- as.character(WHS$YEAR)

head(WHS)

smoothSEcurve <- function(yy, xx) {
  # use after a call to "plot"
  # fit a lowess curve and 95% confidence interval curve
  xx.list <- min(xx) + ((0:100)/100)*(max(xx) - min(xx)) # make list of x values
  
  # Then fit loess function through the points (xx, yy) at the listed values
  yy.xx <- predict(loess(yy ~ xx), se=T, newdata=data.frame(xx=xx.list))
  
  lines(yy.xx$fit ~ xx.list, lwd=2)
  lines(yy.xx$fit - qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
  lines(yy.xx$fit + qt(0.975, yy.xx$df)*yy.xx$se.fit ~ xx.list, lty=2)
}


null_model<-coxph(Surv(LENFOL, FSTAT)~1, data = WHS)
summary(null_model)

par(mfrow=c(2,2))


WHS$cVD <- as.factor(WHS$CVD)
rr<-residuals(null_model, type = "martingale")
attach(WHS)
plot(rr~CVD) ##CVD is constant against residuals. 
lines(lowess(rr~WHS$CVD), col="maroon")
title("Martingale residuals vs CVD")



plot(rr~AGE) ##Age is fine. Likely doesn't need to be transformed. 
smoothSEcurve(rr,WHS$AGE)
title ("Martingale residuals vs Age")

logage <- log(AGE)
plot(rr~logage)
smoothSEcurve(rr,logage)
title("Martingale residuals vs logage")

cox.cVD<- coxph(Surv(LENFOL, FSTAT)~CVD, data=WHS)
summary(cox.cVD)

survreg.CVD<-survreg(Surv(LENFOL, FSTAT)~CVD, data=WHS, dist = "weibull")
summary(survreg.CVD)

0.6159/2.14

##Deciding on forms of variables to be included in model. plot variable vs martingale residuals.
plot(rr~WHS$GENDER)
smoothSEcurve(rr, WHS$GENDER)
lines(lowess(rr~WHS$GENDER), col="maroon") ##linear relationship 

##Heart Rate 
plot(rr~WHS$HR)
smoothSEcurve(rr, WHS$HR)

logHR<-log(WHS$HR)
plot(rr~logHR)
smoothSEcurve(rr, logHR) ##log transformation is better

attach(WHS)
plot(SYSBP, DIASBP)
lines(lowess(SYSBP, DIASBP), col="red") ##SYSBP, DIASBP correlated. 

plot(rr~AFB)  ##AFB
lines(lowess(rr~AFB), col="maroon") ##linear r/nship

plot(rr~SHO)
lines(lowess(rr~SHO), col="maroon") ##linear

plot(rr~CHF)
lines(lowess(rr~CHF), col="maroon") ##linear

plot(rr~AV3)
lines(lowess(rr~AV3), col="maroon") ##linear

plot(rr~MIORD)
lines(lowess(rr~MIORD), col="maroon") #linear

plot(rr~MITYPE)
lines(lowess(rr~MITYPE), col="maroon") #linear

plot(rr~YEAR)
lines(lowess(rr~YEAR), col="maroon") ##linear or constatnt

summary(SYSBP)
summary(DIASBP)

##1 = normal/elevated, 2=stage 1 hypertension, 3=stage 2 hypertension

WHS2 <- WHS %>%
  mutate(hyp = case_when(SYSBP >= 140 | DIASBP >=90 ~ '3',
                         SYSBP >=130 | DIASBP >=80 ~ '2',
                         SYSBP <= 129 & DIASBP <80 ~ '1'))

null_model2<-coxph(Surv(LENFOL, FSTAT)~1, data = WHS2)
rr2<-residuals(null_model2, type="martingale")

attach(WHS2)
WHS2$hyp<-as.character(WHS2$hyp)
plot(rr2~WHS2$hyp)
lines(lowess(rr~WHS2$hyp), col="maroon") ##linear 

plot(rr~BMI)
lines(lowess(rr~BMI), col="maroon")

logBMI <- log(BMI)
plot(rr~logBMI)
lines(lowess(rr~logBMI), col="maroon") ##Better after transformation. 



##None need to be transformed. May need to use new hypertension variable instead of systolic and diastolic bp. 
##Univariate models maybe?
age_model<- coxph(Surv(LENFOL, FSTAT)~AGE, data=WHS2)
summary(age_model) ##significant

gend_model<-coxph(Surv(LENFOL, FSTAT)~GENDER, data=WHS2)
summary(gend_model) ##significant

HR_model<-coxph(Surv(LENFOL, FSTAT)~logHR, data=WHS2)
summary(HR_model) ##significant

hyp_model <- coxph(Surv(LENFOL, FSTAT)~hyp, data=WHS2)
summary(hyp_model) ##hyp3 significant 

BM_model <-coxph(Surv(LENFOL, FSTAT)~BMI, data=WHS2)
summary(BM_model) ##significant

CVD_model <-coxph(Surv(LENFOL, FSTAT)~CVD, data=WHS2)
summary(CVD_model) ##p=0.0922

AFB_model <-coxph(Surv(LENFOL, FSTAT)~AFB, data=WHS2)
summary(AFB_model) ##significant 

SHO_model<-coxph(Surv(LENFOL, FSTAT)~SHO, data=WHS2)
summary(SHO_model) ##significant 

CHF_model<-coxph(Surv(LENFOL, FSTAT)~CHF, data=WHS2)
summary(CHF_model) ##significant

##**AV3_model<- coxph(Surv(LENFOL, FSTAT)~AV3, data=WHS2)
summary(AV3_model) ##p=0.219 ***********

ORD_model<-coxph(Surv(LENFOL, FSTAT)~MIORD, data=WHS2)
summary(ORD_model) ##significant 

type_model <- coxph(Surv(LENFOL, FSTAT)~MITYPE, data=WHS2)
summary(type_model)  ##significant 

yr_model<- coxph(Surv(LENFOL, FSTAT)~YEAR, data=WHS2)
summary(yr_model)   ##year3 significant 


WHS2$YEAR <- as.character(WHS2$YEAR)
Sat_model <- coxph(Surv(LENFOL, FSTAT)~ AGE + GENDER + logHR + hyp + logBMI + CVD + AFB + SHO + CHF + MIORD + MITYPE + YEAR,  data=WHS2)
summary(Sat_model)

result.step <- step(Sat_model, scope = list(upper=~ AGE + 
                                              GENDER + logHR + hyp + logBMI + 
                                              CVD + AFB + SHO + CHF + MIORD + MITYPE +
                                              YEAR, lower = ~AGE) )




prelim.maineffects_model<- coxph(Surv(LENFOL, FSTAT)~ AGE + GENDER + logHR + logBMI + SHO + CHF + MITYPE + YEAR , data=WHS2)
summary(prelim.maineffects_model)

##Check for interactions. Perhaps age and logHR have an interaction. Initial AMI type could affect complications while hospitalized?
em_model<- coxph(Surv(LENFOL, FSTAT)~ AGE + GENDER + logHR + logBMI + SHO + CHF + MITYPE + YEAR + AGE*logHR + MITYPE*SHO + MITYPE*CHF, data=WHS2)
summary(em_model)

sat.em.model<- coxph(Surv(LENFOL, FSTAT)~ AGE + GENDER + logHR + logBMI + SHO + CHF + MITYPE + YEAR + AGE*logHR + MITYPE*SHO + MITYPE*CHF + AGE*GENDER + 
                       AGE*logBMI + AGE*SHO + AGE*CHF + AGE*MITYPE + AGE*YEAR + GENDER*logHR + GENDER*logBMI + GENDER*SHO + GENDER*CHF +
                       GENDER*MITYPE + GENDER*YEAR + logHR*logBMI + logHR*SHO + logHR*CHF + logHR*MITYPE + logHR*YEAR + logBMI*SHO +
                       logBMI*CHF + logBMI*MITYPE + logBMI*YEAR + SHO*CHF + SHO*YEAR + CHF*YEAR + MITYPE*YEAR, data=WHS2)

summary(sat.em.model)

result.step.em<-step(sat.em.model, scope = list(upper=~ AGE + GENDER + logHR + logBMI + SHO + CHF + MITYPE + YEAR +
                                                  AGE*logHR + MITYPE*SHO + MITYPE*CHF + AGE*GENDER + AGE*logBMI + AGE*SHO + AGE*CHF + AGE*MITYPE + AGE*YEAR + GENDER*logHR + GENDER*logBMI + 
                                                  GENDER*SHO + GENDER*CHF + GENDER*MITYPE + GENDER*YEAR + logHR*logBMI + logHR*SHO + logHR*CHF + logHR*MITYPE + logHR*YEAR + logBMI*SHO +
                                                  logBMI*CHF + logBMI*MITYPE + logBMI*YEAR + SHO*CHF + SHO*YEAR + CHF*YEAR + MITYPE*YEAR, lower = ~AGE + GENDER + YEAR))

result.em.model<-coxph(Surv(LENFOL, FSTAT) ~ AGE + GENDER + logHR + logBMI + SHO + CHF + MITYPE + YEAR + SHO:MITYPE + AGE:GENDER + AGE:CHF + GENDER:SHO + 
                         logHR:MITYPE + logBMI:SHO + logHR:CHF, data=WHS2)

summary(result.em.model)



















