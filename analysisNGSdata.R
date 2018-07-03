##############################################
## ---------------------------------------- ##
## ---- Analysis of NGS data Nivolumab ---- ##
## ---------------------------------------- ##
##############################################


## Set path to working directory

path.wd <- "my.working.directory"
setwd(path.wd)

## Load required packages

library(xlsx)
library(survminer)
library(survivalROC)
library(glmnet)
library(survival)
library(compareC)

## Import data

myData <- read.xlsx(file = "dataNGS.xlsx", row.names = T,
                    sheetIndex = 1, stringsAsFactors = F)
attach(myData)


## ------------------------------------------
## Analysis of clinical data ----------------
## ------------------------------------------


## Sex

table(Sex)
prop.table(table(Sex))

table(Sex, bestRECIST)
prop.table(table(Sex, bestRECIST),2)

chisq.test(table(Sex, bestRECIST))

## Age

summary(AgeNivo)
sd(AgeNivo)

summary(AgeNivo[bestRECIST=="PD"])
sd(AgeNivo[bestRECIST=="PD"])

summary(AgeNivo[bestRECIST!="PD"])
sd(AgeNivo[bestRECIST!="PD"])

wilcox.test(AgeNivo~bestRECIST)

## Smoker

table(Smoker)
prop.table(table(Smoker))

table(Smoker, bestRECIST)
prop.table(table(Smoker, bestRECIST),2)

fisher.test(table(Smoker, bestRECIST))

## WHO performance status

table(WHO.performance.status)
prop.table(table(WHO.performance.status))

table(WHO.performance.status, bestRECIST)
prop.table(table(WHO.performance.status, bestRECIST),2)

chisq.test(table(WHO.performance.status, bestRECIST))

## Histology

table(Histology)
prop.table(table(Histology))

table(Histology, bestRECIST)
prop.table(table(Histology, bestRECIST),2)

chisq.test(table(Histology[Histology%in%c("non-squamous cell", "squamous cell")], 
                 bestRECIST[Histology%in%c("non-squamous cell", "squamous cell")]))

## Tumor stage

table(Tumor.stage)
prop.table(table(Tumor.stage))

table(Tumor.stage, bestRECIST)
prop.table(table(Tumor.stage, bestRECIST),2)

fisher.test(table(Tumor.stage, bestRECIST))

## EGFR mutation

table(EGFR.mutation)
prop.table(table(EGFR.mutation))

table(EGFR.mutation, bestRECIST)
prop.table(table(EGFR.mutation, bestRECIST),2)

fisher.test(table(EGFR.mutation, bestRECIST))

## KRAS mutation

table(KRAS.mutation)
prop.table(table(KRAS.mutation))

table(KRAS.mutation, bestRECIST)
prop.table(table(KRAS.mutation, bestRECIST),2)

fisher.test(table(KRAS.mutation, bestRECIST))

## ALK rearrangement

table(ALK.rearrangement)
prop.table(table(ALK.rearrangement))

table(ALK.rearrangement, bestRECIST)
prop.table(table(ALK.rearrangement, bestRECIST),2)

## PD-L1 expression (Sp142 antibody)

table(PDL1_Sp142)
prop.table(table(PDL1_Sp142))

table(PDL1_Sp142, bestRECIST)
prop.table(table(PDL1_Sp142, bestRECIST),2)

chisq.test(table(PDL1_Sp142, bestRECIST))


## ------------------------------------------
## Analysis of NGS data ---------------------
## ------------------------------------------


## Univariate analyses ----------------------
## ------------------------------------------

## TMB per Mb

fit <- survfit(Surv(PFS, evtPFS)~TMB.per.Mb, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~TMB.per.Mb, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 1A

fit <- survfit(Surv(PFS, evtPFS)~Signature.1A, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.1A, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 1B

fit <- survfit(Surv(PFS, evtPFS)~Signature.1B, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.1B, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 3

fit <- survfit(Surv(PFS, evtPFS)~Signature.3, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.3, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 4

fit <- survfit(Surv(PFS, evtPFS)~Signature.4, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.4, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 6

fit <- survfit(Surv(PFS, evtPFS)~Signature.6, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.6, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 16

fit <- survfit(Surv(PFS, evtPFS)~Signature.16, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.16, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Signature 21

fit <- survfit(Surv(PFS, evtPFS)~Signature.21, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~Signature.21, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Number of LOH > 15mB

fit <- survfit(Surv(PFS, evtPFS)~nbLOH15, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~nbLOH15, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Number of deletions > 15mB

fit <- survfit(Surv(PFS, evtPFS)~nbDel15, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~nbDel15, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Ploidy

fit <- survfit(Surv(PFS, evtPFS)~ploidy, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~ploidy, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Number of CNA clusters

fit <- survfit(Surv(PFS, evtPFS)~nbCNAclusters, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~nbCNAclusters, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## MAPK-PI3K-AKT pathway alterations

fit <- survfit(Surv(PFS, evtPFS)~MAPK.PI3K.AKT.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~MAPK.PI3K.AKT.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## WNT pathway alterations

fit <- survfit(Surv(PFS, evtPFS)~WNT.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~WNT.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## DNA repair pathway alterations

fit <- survfit(Surv(PFS, evtPFS)~DNA.repair.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~DNA.repair.pathway, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Number of neopeptides

fit <- survfit(Surv(PFS, evtPFS)~nbNeopeptides, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~nbNeopeptides, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Number of TCR clusters

fit <- survfit(Surv(PFS, evtPFS)~as.factor(nbTCRclusters=="high"), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(PFS, evtPFS)~as.factor(nbTCRclusters!="low"), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fit <- survfit(Surv(OS, evtOS)~as.factor(nbTCRclusters!="low"), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)


## Multivariate analyses --------------------
## ------------------------------------------

detach(myData)
myData <- na.omit(myData[,-30])
attach(myData)

## Clinical model

fit.clinical <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+WHO.performance.status+Histology+Tumor.stage+Smoker+PDL1_Sp142)
val_pred <- fit.clinical$linear.predictors
time_dep.clinical <- survivalROC(PFS, evtPFS, marker=val_pred,
                         predict.time=6, method="KM")
time_dep.clinical$AUC

## TMB model

fit.tmb <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+WHO.performance.status+Histology+Tumor.stage+Smoker+PDL1_Sp142+TMB.per.Mb)
val_pred <- fit.tmb$linear.predictors
time_dep.tmb <- survivalROC(PFS, evtPFS, marker=val_pred,
                        predict.time=6, method="KM") 
time_dep.tmb$AUC

## Exome derived model

myData[myData=="low"] <- 0
myData[myData=="diploid"] <- 0
myData[myData=="monoclonal"] <- 0
myData[myData=="high"] <- 1
myData[myData=="aneuploid"] <- 1
myData[myData=="polyclonal"] <- 1
myData[myData=="none"] <- 1

set.seed(1)
fit.glmnet <- cv.glmnet(x = data.matrix(myData[,1:21]), 
                        y = Surv(myData$PFS, myData$evtPFS),
                        family = "cox")
plot(fit.glmnet)
coef(fit.glmnet, s = "lambda.1se")

val_pred <- predict(fit.glmnet, newx = data.matrix(myData[,1:21]), s = "lambda.1se")
time_dep.glmnet <- survivalROC(PFS, evtPFS, marker=val_pred,
                               predict.time=6, method="KM") 
time_dep.glmnet$AUC


## Compare C indices ------------------------
## ------------------------------------------


compareC(PFS, evtPFS,
         fit.clinical$linear.predictors,
         fit.tmb$linear.predictors)
compareC(PFS, evtPFS,
         fit.clinical$linear.predictors,
         val_pred)
compareC(PFS, evtPFS,
         fit.tmb$linear.predictors,
         val_pred)

## Cutoff -----------------------------------
## ------------------------------------------

fit <- survfit(Surv(PFS, evtPFS)~as.factor(val_pred>-0.3383), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)
