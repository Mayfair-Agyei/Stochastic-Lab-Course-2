install.packages("survival")
install.packages("raster")
install.packages("survminer")
install.packages("ggfortify")
library(survival)
library(raster)
library(survminer)
library(tidyverse)
library(ggfortify)

# 8a

Thoracic <- read.delim("R/Stochastic Lab Course/Thoracic.txt", sep = " ")

colnames(Thoracic) = c("DGN", "PRE4", "PRE5", "PRE6","PRE7","PRE8","PRE9","PRE10",
              "PRE11","PRE14","PRE17","PRE19","PRE25","PRE30","PRE32","AGE","Risk1Y")


Survival.Thoracic <- Thoracic %>% select(PRE30, AGE, Risk1Y)
Survival.Thoracic$PRE30 <- as.integer(as.logical(Survival.Thoracic$PRE30))
Survival.Thoracic$Risk1Y <- as.integer(as.logical(Survival.Thoracic$Risk1Y ))

#Defining time,event and groups
Time <- Survival.Thoracic$AGE 
Event <- Survival.Thoracic$Risk1Y 
Groups <- Survival.Thoracic$PRE30  

#Computing nonparametric estimators of the survivor function: Kaplan-Meier
Kaplan.Meier <- survfit(Surv(Time, Event) ~ 1, data = Thoracic, type = "kaplan-meier")

#Computing nonparametric estimators of the survivor function: Flemming-Harrington
Fleming.Harrington  <- survfit(Surv(Time, Event) ~ 1, data = Thoracic, type = "fleming-harrington")

Fit <- list(Kaplan.Meier = Kaplan.Meier, Fleming.Harrington = Fleming.Harrington)


plot(Kaplan.Meier, col = "red", ylab ="Survival Probabilty", xlab="Time")
lines(Fleming.Harrington, col="yellow")
legend("bottomleft", inset=.04,
     legend = c("Kaplan-Meier", "Fleming-Harrington"),
        lwd=c(1,1), lty=c(1,1),
        col=c("red", "yellow"))


#Fitting the exponential model to the data
ExponentialModel.Fit <- survreg(Surv(Time, Event) ~ 1, data = Thoracic, dist="exponential")
Lam.Exponential <- exp(-ExponentialModel.Fit$coefficients)
Exponential.Fit <- exp(-Lam.Exponential*c(0:90))

#Fitting the Weibull model to the data
WeibullModel.Fit <- survreg(Surv(AGE, Risk1Y) ~ 1, data = Thoracic, dist="weibull")
Lam.Weibull <- exp(-WeibullModel.Fit$coefficients)
Alpha <- 1/WeibullModel.Fit$scale
Weibull.Fit <- exp(-(Lam.Weibull * c(0:90))^Alpha)
par(mfrow=c(1,1))

# plot
plot(Exponential.Fit, col = "red") 
lines(Weibull.Fit, col = "yellow") 
lines(Kaplan.Meier, col = "green") 
legend("bottomleft", inset=.04,
       legend =  c("Exponential", "Weibull", "Kaplan-Meier"),
       col =  c("red", "yellow", "green"),
       lwd=c(2,2), lty=c(2,2),
       box.lty=0) +
title(main="Parametric estimators of the survivor function with Kaplan-Meier")

#Checking if weibull is adequate for the data
plot(log(Kaplan.Meier$time) ,log(-log(Kaplan.Meier$surv)), col = "red")
abline(a = Alpha*log(Lam.Weibull), b = Alpha, col = "green" )

# 8b
# COmputing the Kaplan-Meier estimators for each group
fit <- survfit(Surv(Time, Event) ~ PRE30, data = Thoracic)
ggsurvplot(fit, data = Thoracic)

#Fitting the Weibull model to smokers' group
Smokers <- Thoracic %>% filter(PRE30 == "TRUE")
WeibullSmokers.Fit <- survreg(Surv(Time, Event) ~ 1, data = Smokers, dist="weibull")
Lam.WeibullSmokers <- exp(-WeibullSmokers.Fit$coefficients)
Alpha.Smokers <- 1/WeibullSmokers.Fit$scale
WeibullSmokers.FitData <- exp(-(Lam.WeibullSmokers*c(0:90))^Alpha.Smokers)

#Fitting the Weibull model to nonsmokers' group
Nonsmokers <- Thoracic %>% filter(PRE30 == "FALSE")
WeibullNonsmokers.Fit <- survreg(Surv(Time, Event) ~ 1, data = Nonsmokers, dist="weibull")
Lam.WeibullNonsmokers<- exp(-WeibullNonsmokers.Fit$coefficients)
Alpha.Nonsmokers <- 1/WeibullNonsmokers.Fit$scale
WeibullNonsmokers.FitData <- exp(-(Lam.WeibullNonsmokers*c(0:90))^Alpha.Nonsmokers)
sum(Thoracic$PRE30)/nrow(Thoracic) * 100

#Computing nonparametric estimators of the survivor function: Kaplan-Meier
Kaplan.MeierSmokers <- survfit(Surv(Time, Event) ~ 1, data = Smokers, type = "kaplan-meier")
Kaplan.MeierNonsmokers <- survfit(Surv(Time, Event) ~ 1, data = Nonsmokers, type = "kaplan-meier")

#plot
plot(WeibullSmokers.FitData, col = "yellow", type = "l")
lines(WeibullNonsmokers.FitData, col='red')
lines(Kaplan.MeierSmokers, col = "green")
lines(Kaplan.MeierNonsmokers, col = "blue")
legend("bottomleft", inset=.0001,
       legend =  c("Weibull smokers", "Weibull nonsmokers", "Kaplan-Meier smokers", "Kaplan-Meier nonsmokers"),
       col =  c("yellow", "red", "green", "blue"),
       lwd=c(1,1), lty=c(1,1),
       box.lty=0)
title(main="Nonparametric & Parametric estimators \n of the survivor function by group")

#Checking if the Weibull model an appropriate assumption in both groups
plot(log(Kaplan.MeierSmokers$time) ,log(-log(Kaplan.MeierSmokers$surv)), col = "red")
abline(a = Alpha.Smokers*log(Lam.WeibullSmokers), b = Alpha.Smokers, col = "blue" )
title(main="Checking if weibull is appropriate \n for smokers' group")
plot(log(Kaplan.MeierNonsmokers$time) ,log(-log(Kaplan.MeierNonsmokers$surv)), col = "red")
abline(a = Alpha.Nonsmokers*log(Lam.WeibullNonsmokers), b = Alpha.Nonsmokers, col = "blue" )
title(main="Checking if weibull is appropriate \n for both group")

