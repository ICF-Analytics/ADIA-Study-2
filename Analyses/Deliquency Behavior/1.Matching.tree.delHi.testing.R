options(digits = 3)
options(scipen = 1000000)
set.seed(0203)

library(tidyverse)
library(dplyr)
library(designmatch)
# Zubizarreta et al. (2011) https://doi.org/10.1198/tas.2011.11072
# Kelz et al. (2013) https://doi.org/10.1097/sla.0b013e31829654f3
library(cobalt) # to asesss/plot balance
library(evtree)
library(rpart)
library(partykit)
library(rpart.plot)
library(survey)


setwd("C:/Users/21983/OneDrive - ICF/ADIA/study 2")
outv <- "RB_deliq_hi"
#' RB_deliq and RB_deliq_Mod do not enough outcome for match control

#===============================================================================

dat <- read.csv("data/LS_all.csv")
head(dat)
colnames(dat)
summary(dat)

#===============================================================================
### I.- Developing sample
#===============================================================================
### 1.- We focus on population expericing ACE during childhood
# (and before the otucome was asssess)
dat$anyACE <- as.numeric(apply(
  dat[, c("M", "MP", "MS", "MN", "MF", "ML", "ME", "MM", "MD")], 1, sum) >= 1)
table(dat$anyACE) #0=383 1=971
dat <- dat[dat$anyACE == 1, ] #n=971

###' 2.- In that population, we want to campare
#' those who develop a certain outcome to those who do not
#' "case-control" design)
#' table(dat$RB_deliq, useNA="ifany") #0=317 1=318 na=336
#' table(dat$RB_deliq_Mod, useNA="ifany") #0=413 1=222 na=336
table(dat$RB_deliq_hi, useNA = "ifany") #0=513 1=121 na=337
dat$y  <- dat$RB_deliq_hi
table(dat$y, useNA = "ifany")
dat <- dat[!is.na(dat$y), ]


### 3.- "holding constant" some characteristics
dat$cohort <- as.numeric(substr(dat$ChildDOB, 6, 9))
summary(dat$cohort)

xn <- c("center", "childrace_bsl",
  "cohort", "childsex_bsl", "caregiver_married16", "hh_income16")
#3.0 treat missing as a separate category
#(so we do not lose cases becosue of missing covariates here)
dat[, xn] <- lapply(dat[, xn], addNA, ifany = TRUE) #
summary(dat[, xn])
apply(is.na(dat[, xn]), 2, table)




#-------------------------------------------------------------------------------
### 4.- Finding mathes

#i.- Group indicator
# in general "treatment" indicator
# in "case-control" studies, outcome indicator
# note that the data needs to be sorted in decreasing order
dat <- dat[order(dat$y, decreasing = TRUE),]
t_ind = dat$y

#ii.- Distance matrix-----------------------------------------------------------
# A measure of similarity/disimilarity in the covariates
# to shose individual matches

#x <- c("center", "childrace_bsl", "cohort",
# "childsex_bsl", "caregiver_married16", "hh_income16")
x <- c("center", "childrace_bsl", "childsex_bsl",
  "caregiver_married16", "hh_income16")
#' take out cohort from pair-matching, still in mom bal
X         <- model.matrix(~ ., data = dat[, x])[, -1]
dist_mat = distmat(t_ind, X)

#iii.  Moment balance-----------------------------------------------------------
# constrain differences in means to be at most .05 standard deviations apart
# this influence the marginal distribution not the individual matches

z <- c("center", "childrace_bsl", "cohort",
  "childsex_bsl", "caregiver_married16", "hh_income16")
Z         <- model.matrix(~ ., data = dat[, z])[, -1]
mom_tols = round(absstddif(Z, t_ind, .05), 2)
#' if there is a rare factors (e.g., NA),
#' either take it out from X, or
#' change tol .05 to larger number .1, .15

mom = list(covs = Z, tols = mom_tols)


#iv.- Subset matching weight----------------------------------------------------
# this is a tunning parameter that regulate tredeoff between
# number of matches and quality of the matches
subset_weight =  median(dist_mat)
#' this parameter regulate matches,
#' if qulaity getting bad or drop a lot of cases, 
#' we can change

#! Ye: there are some other parameters that you can adjust
#! see help(bmatch) and/or reference
#! e.g., if we go for diffrent waves, we may want to match exactly on
#! hisotrical values of the outcome (e.g., y12 if we are looking at y14)

# Fine balance, same as moment blance, but for factors--------------------------

#' fine_covs =  model.matrix(~ ., data=dat[,"center"])
#' fine = list(covs = fine_covs)

# Exact matching
#' exact_covs = cbind(black)
#' exact = list(covs = exact_covs)

#v.- Solver options-------------------------------------------------------------
t_max = 60 * 5
solver = "glpk"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace = 0)

#vi.-  Match
out <- bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,
              mom = mom,  solver = solver)


### 5.-  assess balance---------------------------------------------------------
baltab <- bal.tab(out, treat = t_ind, covs = X
  , s.d.denom = "treated", un = TRUE, disp = c("means", "sds"))
print(baltab, un = FALSE, disp.v.threshold = FALSE)

lplot <- love.plot(out, treat = t_ind, covs = X
  , thresholds = c(m = .1), binary = "std", s.d.denom = "treated")
lplot

png(paste0("output/", outv, ".love.plot.png")
  , width = 480 * 4, heigh = 480 * 4, res = 300)
lplot
dev.off()


#proportion of variance explained by matching
wv <- \(x, id) {
  x <- as.numeric(x)
  dx <- x - ave(x, id, FUN = mean)
  ssb <- sum((ave(x, id, FUN = mean) - mean(x))^2)
  ssw <- sum(dx^2)
  dfb <- length(unique(id)) - 1
  dfw <- length(unique(id))
  msb <- ssb / dfb
  msw <- ssb / dfw
  ssb / (ssb + ssw)
  }

with(sdat, wv(childsex_bsl, id))
with(sdat, wv(cohort, id))
with(sdat, wv(hh_income16, id))

table(out$group_id, useNA = "ifany")
sdat <- dat[c(out$t_id, out$c_id), ]
dim(sdat)
sdat$id <- out$group_id
saveRDS(sdat, file = paste0("data/", outv, ".m.Rds"))