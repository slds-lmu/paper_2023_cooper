### penalizedCSC-predictions.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 11 2022 (18:39)
## Version:
## Last-Updated: May 11 2022 (18:48)
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
##
### Commentary:
## Here we hack our own coefficients into the Cox regression models
## of a CSC object for prediction of the risk of event 1 in the presence
## of a competing risk
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

library(riskRegression)
library(survival)
set.seed(8)
dlarge <- sampleData(1000,outcome="competing.risk")
dsmall <- sampleData(100,outcome="competing.risk")
dtest <- sampleData(400,outcome="competing.risk")

F0 = CSC(Hist(time, event) ~ X1+X2+X7+X9, data = dlarge, singular.ok = TRUE)
F1 = CSC(Hist(time, event) ~ X1+X2+X7+X9, data = dsmall, singular.ok = TRUE)
F2 = F1

# coefficients from large model
F1$models[[1]]$coefficients <- F0$models[[1]]$coefficients
F1$models[[2]]$coefficients <- F0$models[[2]]$coefficients
# funny coefficients
F2$models[[1]]$coefficients <- c(-.3,-2,0,-1)
F2$models[[2]]$coefficients <- c(.1,.1,.1,.1)

# F0: has the coefficients from the large dataset and
#    the baseline hazard from the large dataset
# F1: has the coefficients from the large dataset and
#    the baseline hazard from the small dataset
# F2: has the funny coefficients and
#    the baseline hazard from the small dataset
predict(F0,times = 5,newdata = dtest[1:5,])
predict(F1,times = 5,newdata = dtest[1:5,])
predict(F2,times = 5,newdata = dtest[1:5,])

predictRisk(F1,times = 5,newdata = dtest[1:5,])
predictRisk(F2,times = 5,newdata = dtest[1:5,])

# evaluation of predicted risks in test data
x <- Score(list(F0 = F0,F1 = F1,F2 = F2),formula = Hist(time,event)~1,data = dtest)
# we see that the AUC is independent of the baseline hazard and hence
# the same for F0 and F1, but the Brier score is smaller for F1
summary(x)

x$Brier$score
x$Brier$contrasts

######################################################################
### penalizedCSC-predictions.R ends here
