options(digits = 3)
options(scipen = 1000000)

library(dplyr)
library(designmatch)
# Zubizarreta et al. (2011) https://doi.org/10.1198/tas.2011.11072
# Kelz et al. (2013) https://doi.org/10.1097/sla.0b013e31829654f3
library(cobalt) # to asesss/plot balance
library(tidyverse)
library(rpart)
library(partykit)
library(rpart.plot)
library(evtree)
library(survey)
library(gt)


# Ye: un-comment this next line to change path to your own
setwd('C:/Users/21983/OneDrive - ICF/ADIA/study 2') 

outv <- "RB_su16"

#===============================================================================

dat <- read.csv('data/LS_all_c1.csv')
head(dat)
colnames(dat)
summary(dat)
#===============================================================================
### I.- Developing sample
#===============================================================================
### 1.- compute variables
dat$cohort <- as.numeric(substr(dat$ChildDOB,6,9))
summary(dat$cohort)

### 2.- We focus on population expericing ACE during childhood
# (and before the otucome was asssess)
dat$anyACE <- as.numeric(apply(
  dat[, c('M','MP','MS','MN','MF','ML','ME','MM','MD')], 1, sum) >= 1)
table(dat$anyACE) #0=383 1=971
dat <- dat[dat$anyACE==1,] #n=971

### 3.- In that population, we want to campare those who develop a certain outcome to those who do not
# ()'case-control' design)

dat$y <- dat$RB_su #exposed and no missing on outcome 
table(dat$y,useNA='ifany') #0=256 1=481 na=234
dat <- dat[!is.na(dat$y), ] #n=737

### 4.- 'holding constant' some characteristics

x.n <- c('center', 'childrace_bsl','cohort','childsex_bsl','caregiver_married16','hh_income16')
#4.0 treat missing as a separate category
#(so we do not lose cases becosue of missing covariates here)
dat[,x.n] <- lapply(dat[,x.n] ,addNA,ifany =T) #
summary(dat[,x.n])
apply(is.na(dat[,x.n]),2,table)

#-------------------------------------------------------------------------------
### 5.- Finding matches

#i.- Group indicator
# in general 'treatment' indicator
# in 'case-control' studies, outcome indicator
# note that the data needs to be sorted in decreasing order
dat <- dat[order(dat$y,decreasing =T),]
t_ind = dat$y

#ii.- Distance matrix-----------------------------------------------------------
# A measure of similarity/disimilarity in the covariates
# to shose individual matches

#x <- c('center', 'childrace_bsl','cohort','childsex_bsl','caregiver_married16','hh_income16')
x <- c('center', 'childrace_bsl','childsex_bsl','caregiver_married16','hh_income16')#take out cohort from pair-matching, still in mom bal
X         <- model.matrix(~ .,data=dat[,x])[,-1]
dist_mat = distmat(t_ind, X)

#iii.  Moment balance-----------------------------------------------------------
# constrain differences in means to be at most .05 standard deviations apart
# this influence the marginal distribution not the individual matches

z <- c('center', 'childrace_bsl','cohort','childsex_bsl','caregiver_married16','hh_income16')
Z         <- model.matrix(~ .,data=dat[,z])[,-1]
mom_tols = round(absstddif(Z, t_ind, .05), 2)#if there is a rare factors (e.g., NA), either take it out from X, or change .05 to larger number .1, .15
mom = list(covs = Z, tols = mom_tols)


#iv.- Subset matching weight----------------------------------------------------
# this is a tunning parameter that regulate tredeoff between
# number of matches and quality of the matches
subset_weight =  median(dist_mat) #this parameter regulate matches, if qulaity getting bad or drop a lot of cases, can change

#! Ye: there are some other parameters that you can adjust
#! see help(bmatch) and/or reference
#! e.g., if we go for diffrent waves, we may want to match exactly on
#! hisotrical values of the outcome (e.g., y12 if we are looking at y14)

# Fine balance, same as moment blance, but for factors--------------------------

#fine_covs =  model.matrix(~ ., data=dat[,'center'])
#fine = list(covs = fine_covs)

# Exact matching
#exact_covs = cbind(black)
#exact = list(covs = exact_covs)

#v.- Solver options-------------------------------------------------------------
t_max = 60*5
solver = "glpk"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace = 0)

#vi.-  Match
out <- bmatch(t_ind = t_ind, dist_mat = dist_mat, subset_weight = subset_weight,
              mom = mom,  solver = solver)


### 5.-  assess balance---------------------------------------------------------
baltab<-bal.tab(out,treat=t_ind,covs=X,s.d.denom='treated',un = TRUE,disp = c("means", "sds"))
print(baltab, un = FALSE, disp.v.threshold = FALSE)


lplot <- love.plot(out,treat=t_ind,covs=X,thresholds = c(m = .1), binary = "std",s.d.denom='treated')
lplot
png(paste0("output/",outv,".love.plot.png"), width = 480*4,heigh=480*4,res=300)
lplot
dev.off()


table(out$group_id,useNA='ifany')
sdat <- dat[c(out$t_id,out$c_id	),]
dim(sdat)
sdat$id <- out$group_id
saveRDS(sdat, file=paste0('data/', outv, '.m.Rds'))
#C:/Users/21983/OneDrive - ICF/ADIA/study 2
#write.csv(sdat,file=file.path(pathi,'anySUD_matched.csv'),na='')
#Ye: you may watn to save files and plots

