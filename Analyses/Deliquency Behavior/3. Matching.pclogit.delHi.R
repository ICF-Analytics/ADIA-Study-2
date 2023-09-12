options(digits = 3)
options(scipen = 1000000)
set.seed(0203)

library(tidyverse)
library(survival)


setwd("C:/Users/21983/OneDrive - ICF/ADIA/study 2")
outv <- "RB_deliq_hi"
#' RB_deliq and RB_deliq_Mod do not enough outcome for match control

#===============================================================================
###II. Exploring possible protective factors
#===============================================================================
sdat <- readRDS(file = paste0("data/", outv, ".m.Rds"))

#1. define a training sample
length(unique(sdat$id))
set.seed(0203)
ts <- sample(unique(sdat$id), length(unique(sdat$id)) / 2)
sdat$ts <- as.numeric(sdat$id %in% ts)
table(sdat$ts, useNA = "ifany")

#2.- Protective factors
dim(sdat)
summary(sdat)
#ALL PCE indicators
zn <- c("anysupadu_1",  	"anysupparent_1",  	"anysuprelative_1"
  , "anysupnonfam_1",  	"fam_sat_1",  	"home_safety_1"
  , "prrelation_1",  	"neighborhood_exp_1",  	"school_safety_t_1"
  , "srvc_use_1", "childcare_1", "neighborhood_safety_2"
  , "neighborhood_exp_2", "school_safety_y_2", "school_safety_t_2"
  , "srvc_use_2", "anysupadu_3", "anysupparent_3"
  , "anysuprelative_3", "anysupnonfam_3", "prrelation_3"
  , "bestfriend_3", "socialpart_3", "parent_involv_3"
  , "resid_stab_3", "neighborhood_safety_3", "neighborhood_exp_3"
  , "srvc_use_3")

summary(sdat[, zn])


###2.2- centering

#select varaibles and create dummies
sdat$X <-
  model.matrix(~ .
  , model.frame(~ ., sdat[, zn], na.action = na.pass))[, -1] %>%
  data.frame()
head(sdat$X)


library(missRanger)
###single imputation
sdat$iX <- missRanger(sdat$X, pmm.k = 3)
lapply(sdat$iX, table)


###multiple imputation 
#5  for quick fits 20 for final run
filled <- replicate(5, {
    ix <- missRanger(data.frame(sdat$X)
    , num.trees = 50, pmm.k = 5)
    }, simplify = FALSE)

summary(filled)



#===============================================================================
###3. - Run tree algorithm on training data
#===============================================================================
df0 <- #sdat
  sdat %>%
  filter(ts == 1)
dim(df0)
dim(sdat)
set.seed(0203)


#3.3. - penalized conditional logistic
library(clogitL1)
cll1 <- with(df0, clogitL1(iX, y, id, numLambda = 100
  , minLambdaRatio = 0.000001, switch = 0, alpha = 1))

par(mfrow = c(1, 1))
plot(cll1, logX = TRUE)

clcv <- cv.clogitL1(cll1)
plot(clcv)
summary(clcv)

colnames(sdat$X[summary(clcv)$beta_minCV != 0])
colnames(sdat$X[summary(clcv)$beta_minCV1se != 0])

#===============================================================================
###.-inference (95% CI; p-values)
#===============================================================================
df1 <- #sdat
  sdat %>%
  filter(ts == 0)
dim(df1)
dim(sdat)


#condtional logitic for penalized cl
cfit2 <- clogit(y ~
 + school_safety_t_1 + neighborhood_exp_2
 + strata(id)
 , data = sdat)

summary(cfit2)




#using multiple imputed covariates
mean(is.na(df1$school_safety_t_1))
mean(is.na(df1$neighborhood_exp_2))

models <- lapply(filled, \(x)
    clogit(y ~
    + school_safety_t_1 + neighborhood_exp_2
    + strata(id)
    , data = data.frame(y = df1$y, id = df1$id, x[sdat$ts == 0, ])
    ))

library(mice)
summary(pooled_fit <- pool(models))