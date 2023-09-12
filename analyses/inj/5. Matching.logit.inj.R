options(digits = 3)
options(scipen = 1000000)
set.seed(0203)

library(tidyverse)
library(survival)


setwd("C:/Users/21983/OneDrive - ICF/ADIA/study 2")
outv <- "inj_victim"
#' RB_deliq and RB_deliq_Mod do not enough outcome for match control

#===============================================================================
###II. Exploring possible protective factors
#===============================================================================
sdat <- readRDS(file = paste0("data/", outv, ".2.m.Rds"))

#1. define a training sample
length(unique(sdat$id))
set.seed(0203)
ts <- sample(unique(sdat$id), length(unique(sdat$id)) / 2)
sdat$ts <- as.numeric(sdat$id %in% ts)
table(sdat$ts, useNA = "ifany")

#2.- Protective factors
dim(sdat)
summary(sdat)
#PCE across time indicators
zn <- c('anysupadu', 'anysupparent', 'anysuprelative', 'anysupnonfam', 'fam_sat', 'home_safety', 'prrelation',
                  'bestfriend', 'socialpart', 'parent_involv', 'resid_stab', 'neighborhood_safety', 'neighborhood_exp',
                  'school_safety_y', 'school_safety_t','srvc_use', 'childcare')

#zn <- c('anysupadu_s', 'anysupparent_s', 'anysuprelative_s', 'anysupnonfam_s', 'fam_sat_s', 'home_safety_s', 'prrelation_s',
#        'bestfriend', 'socialpart_s', 'parent_involv_s', 'resid_stab', 'neighborhood_safety_s', 'neighborhood_exp_s',
#        'school_safety_y_s', 'school_safety_t_s', 'srvc_use_s', 'childcare')
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
###.-inference (95% CI; p-values)
#===============================================================================
#df1 <- #sdat
#  sdat %>%
#  filter(ts == 0)
#dim(df1)
dim(sdat)


#condtional logitic for penalized cl

vars  <- paste(paste(zn,collapse='+',sep='+'), 'strata(id)', sep='+'); vars
fm  <- formula(paste('y ~', vars))

#fit <- glm(fm,family = quasibinomial, data=sdat,)
#summary(fit)

cfit <- clogit(fm, data = sdat)

summary(cfit)
round(summary(cfit)$coef, 4)

#take out srvc_use and 3 variables that have large se
cfit1 <- clogit(y ~anysupparent+anysupnonfam+anysuprelative+anysupnonfam                 
                +home_safety+bestfriend+socialpart+parent_involv
                +resid_stab+school_safety_y+childcare
                +strata(id)
                , data = sdat)
summary(cfit1)

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