#load the necessary libraries
library(survival)
library(SurvRegCensCov)
library(Hmisc)
library(mice)
library(rms)
library(survminer)
library(ggplot2)
library(ggfortify)

###################################################### Initial data preparation ######################################################################

#read data into object
mgus<-read.csv(file.choose(), header=T)

#confirm object is dataframe
is.data.frame(mgus)

#drop columns that will not be used
vars <- names(mgus) %in% c("id", "X", "death", "futime")
mgus <- mgus[!vars]

#get summary statistics
str(mgus)
summary(mgus)

# Impute missing 
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(mgus,2,pMiss)

mgusImputed <- mice(mgus,m=5,maxit=50,meth='pmm')
summary(mgusImputed)

# Verify imputed data is valid based on distribution of original data points
xyplot(mgusImputed, creat ~ sex+age+hgb+mspike+ptime+pstat,pch=18,cex=1)
xyplot(mgusImputed, hgb ~ sex+age+creat+mspike+ptime+pstat,pch=18,cex=1)
xyplot(mgusImputed, mspike ~ sex+age+hgb+creat+ptime+pstat,pch=18,cex=1)

densityplot(mgusImputed)

stripplot(mgusImputed, pch = 20, cex = 1.2)

############################################### Build initial KM models ####################################################################

#fit survival curve without explanatory variables
mgus.fit<-with(mgusImputed, survfit(Surv(ptime, pstat)~1, error="greenwood", conf.type="log-log", conf.int=0.95))
summary(mgus.fit)
summary(mgus.fit,times=c(81,240,373))

#fit survival curve stratifying by gender
mgus.sex.fit<-with(mgusImputed, survfit(Surv(ptime, pstat)~sex, error="greenwood", conf.type="log-log", conf.int=0.95))
summary(mgus.sex.fit)

#Plot the curves
autoplot(mgus.fit$analyses[1], xlab="Time until progression to a PCM (months)", ylab="Survival Probability")
autoplot(mgus.sex.fit$analyses[1], xlab="Time until progression to a PCM (months)", ylab="Survival Probability")

#Perform the log-rank test with equal weights 
mgus.sex.test0<-with(mgusImputed, survdiff(Surv(ptime, pstat)~sex, rho=0))
mgus.sex.test0

##################################### Visually inspect PH assumptions for categorical variables ###################################################

#Plot a km model with imputed data for gender
plot(mgus.sex.fit$analyses[[1]]$time, log(-log(mgus.sex.fit$analyses[[1]]$surv)),
xlab="Time (months)", ylab="Log-Cumulative Hazard Function for Sex",type="l",lty=1:2) 
legend(350, -5.5, legend=levels(mgus$sex), lty=1:2)


##################################### Build initial model ###########################################################################

#Fit coxph model with imputed data sets and pool the results
mgus.ph.initial <- with(mgusImputed, coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike))
mgus.ph.initial.pooled <- pool(mgus.ph.initial)
summary(mgus.ph.initial.pooled)

###################################### Check for zero correlation with time ##########################################################

#Test each covariate for zero correlation between the time points and the associated sequence of estimates for the regression coefficient
tmp <- with(mgusImputed, cox.zph(coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike)))
tmp

# Plot Schoenfeld test for each covariate
ggcoxzph(tmp$analyses[[1]])

#Check martingale residuals for each covariate
m.resid<-with(mgusImputed, resid(coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike), "mart"))
par(mfrow=c(1,1))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$sex, xlab="Sex", ylab="Martingale residuals")
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$age, xlab="Age", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$age, m.resid$analyses[[1]]))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$creat, xlab="Creatinine Levels", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$creat, m.resid$analyses[[1]]))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$hgb, xlab="Hgb Levels", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$hgb, m.resid$analyses[[1]]))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$mspike, xlab="Monoclonal Serum Levels", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$mspike, m.resid$analyses[[1]]))

# Check Schoenfeld residuals for each covariate
s.resid<-with(mgusImputed, resid(coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike), "scho"))
par(mfrow=c(1,1))
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,1], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for sex")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,2], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for age")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,3], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for hgb")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,4], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for creat")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,5], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for mspike")

# Check for influential outliers
ggcoxdiagnostics(mgus.ph.initial$analyses[[1]], type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(mgus.ph.initial$analyses[[1]], type = "deviance", linear.predictions = FALSE, ggtheme = theme_bw())

############################################### Preform variable selection ########################################################################
# Backwards stepwise variable selection
mgus.ph.var <- with(mgusImputed, step(coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike)))
mgus.ph.var.pooled <- pool(mgus.ph.var)
summary(mgus.ph.var.pooled)

############################################### Check Assumptions again for selected variables ####################################
#Test each covariate for zero correlation between the time points and the associated sequence of estimates for the regression coefficient
tmp <- with(mgusImputed, cox.zph(coxph(Surv(ptime, pstat)~age+hgb+mspike)))
tmp

# Plot Schoenfeld test for each covariate
ggcoxzph(tmp$analyses[[1]])

#Check martingale residuals for each covariate
m.resid<-with(mgusImputed, resid(coxph(Surv(ptime, pstat)~age+hgb+mspike), "mart"))
par(mfrow=c(1,1))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$age, xlab="Age", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$age, m.resid$analyses[[1]]))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$hgb, xlab="Hgb Levels", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$hgb, m.resid$analyses[[1]]))
plot(m.resid$analyses[[1]]~complete(mgusImputed,1)$mspike, xlab="Monoclonal Serum Levels", ylab="Martingale residuals")
lines(lowess(complete(mgusImputed,1)$mspike, m.resid$analyses[[1]]))

# Check Schoenfeld residuals for each covariate
s.resid<-with(mgusImputed, resid(coxph(Surv(ptime, pstat)~sex+age+hgb+creat+mspike), "scho"))
par(mfrow=c(1,1))
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,1], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for age")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,2], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for hgb")
plot(as.numeric(dimnames(s.resid$analyses[[1]])[[1]]), s.resid$analyses[[1]][,3], xlab="Observed Failure Time", ylab="Schoenfeld Residuals for mspike")


############################################### Build final model #################################################################################
#Build final model with imputations
mgus.ph.final <- with(mgusImputed, coxph(Surv(ptime, pstat)~age+hgb+mspike))
mgus.ph.final.pooled <- pool(mgus.ph.final)
summary(mgus.ph.final.pooled)
summary(survfit(mgus.ph.final$analyses[[1]]),times=c(81,240,373))
autoplot(survfit(mgus.ph.final$analyses[[1]]), xlab="Time until progression to a PCM (months)", ylab="Survival Probability")

#Build model without imputations
no_imputations <- coxph(Surv(ptime, pstat)~age+hgb+mspike, data=mgus)
no_imputations.fit <- survfit(no_imputations, error="greenwood", conf.type="log-log", conf.int=0.95)
summary(no_imputations.fit)
summary(no_imputations.fit, times=c(81,240,373))
summary(no_imputations)

autoplot(no_imputations.fit, xlab="Time until progression to a PCM (months)", ylab="Survival Probability")
