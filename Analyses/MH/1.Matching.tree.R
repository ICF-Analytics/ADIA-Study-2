options(digits = 3)
options(scipen = 1000000)
set.seed(0203)

library(dplyr)
library(designmatch)
# Zubizarreta et al. (2011) https://doi.org/10.1198/tas.2011.11072
# Kelz et al. (2013) https://doi.org/10.1097/sla.0b013e31829654f3
library(cobalt) # to asesss/plot balance
library(rpart)
library(partykit)
library(rpart.plot)
library(survey)

getwd()
# Ye: un-comment this next line to change path to your own
setwd('C:/Users/21983/OneDrive - ICF/ADIA/study 2') 

dat <- read.csv('data/LS_MH_PCE_anyt.csv')
head(dat)
colnames(dat)
summary(dat)
#===============================================================================
### I.- Developing sample
#===============================================================================
### 1.- We focus on population expericing ACE during childhood
# (and before the otucome was asssess)
dat$anyACE <- apply(
  dat[, c('M','MP','MS','MN','MF','ML','ME','MM','MD')], 1, sum) >= 1
table(dat$anyACE)
dat <- dat[dat$anyACE,]

### 2.- In that population, we want to campare those who develop a certain outcome to those who do not
# ()'case-control' design)
#dat$y <- dat$anyMH
dat$y <- dat$anysud
table(dat$y,useNA='ifany')
dat <- dat[!is.na(dat$y),]


### 3.- 'holding constant' some characteristics
dat$cohort <- as.numeric(substr(dat$ChildDOB,6,9))
summary(dat$cohort)

x.n <- c('center', 'childrace_bsl','cohort','childsex_bsl','caregiver_married16','hh_income16')
#3.0 treat missing as a separate category
#(so we do not lose cases becosue of missing covariates here)
dat[,x.n] <- lapply(dat[,x.n] ,addNA,ifany =T) #
summary(dat[,x.n])
apply(is.na(dat[,x.n]),2,table)




#-------------------------------------------------------------------------------
### 4.- Finding mathes

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
bal.tab(out,treat=t_ind,covs=X,s.d.denom='treated',un = TRUE,disp = c("means", "sds"))

lplot <- love.plot(out,treat=t_ind,covs=X,thresholds = c(m = .1), binary = "std",s.d.denom='treated')
lplot
png('output/anySUD.love.plot.png', width = 480*4,heigh=480*4,res=300)
lplot
dev.off()


table(out$group_id,useNA='ifany')
sdat <- dat[c(out$t_id,out$c_id	),]
dim(sdat)
sdat$id <- out$group_id
saveRDS(sdat, file='data/anySUD.m.Rds')
#C:/Users/21983/OneDrive - ICF/ADIA/study 2
#write.csv(sdat,file=file.path(pathi,'anySUD_matched.csv'),na='')
#Ye: you may watn to save files and plots

#===============================================================================
###II. Exploring possible protective factors
#===============================================================================
sdat <- readRDS(sdat, file='data/anySUD.m.Rds')

#1. define a training sample
length(unique(sdat$id))
set.seed(0203)
ts <- sample(unique(sdat$id), length(unique(sdat$id)) / 2)
sdat$ts <- sdat$id %in% ts


#2.- Protective factors
dim(sdat)
summary(sdat)
#PCE across time indicators
z.n  <- c('anysupadu', 'anysupparent', 'anysuprelative', 'anysupnonfam', 'fam_sat', 'home_safety', 'prrelation',
          'bestfriend', 'socialpart', 'parent_involv', 'resid_stab', 'neighborhood_safety', 'neighborhood_exp',
          'school_safety_y', 'srvc_use', 'childcare')

summary(sdat[, z.n])


###2.2- centering
# we need to center covariates at the pair-mean
# Stanfill et al. (2019) https://doi.org/10.1177/1179597219858954

#'old school' way
#sdat   <-  sdat[order(sdat$id),]
#sdat$Z <-  model.matrix(~ ., model.frame(~ ., sdat[,z.n], na.action=na.pass))[,-1]
#sdat$Z <-  do.call(rbind,lapply(split(data.frame(Z),sdat$id), function(d) apply(d,2,scale,scale=F) ))

#tider
library(tidyverse)
sdat$Z <-
  sdat[, z.n] %>%
  with(., model.matrix(~ .,
                       model.frame(~ ., data=., na.action = na.pass))[,-1]) %>%
  data.frame() %>%
  group_by(sdat$id) %>%
  mutate(across(everything(), scale, scale = FALSE)) %>%
  ungroup %>%
  .[, -dim(.)[2]]

#-------------------------------------------------------------------------------
###Ye: this is a small simulation test
# It generates a predictor X
# and check that the procedure captures the predictor
#sdat$x <- ifelse(runif(length(sdat$y ))>0.25,!sdat$y,sdat$y)
#table(sdat$x,sdat$y)

#sdat$Z <-
#  sdat[,c(z.n,'x')] %>%
#  with(., model.matrix(~ ., model.frame(~ ., data=.,na.action=na.pass))[,-1]) %>%
#  data.frame ()%>%
#  group_by(sdat$id) %>%
#  mutate(across(everything(),scale,scale=F)) %>%
#  ungroup %>%
#  .[,-dim(.)[2]]
#-------------------------------------------------------------------------------

###3. - Run tree algorithm on training data

#3.1. - evolutionary tree (just for comparisson)
library(evtree)
efit <- with(sdat[sdat$id %in% ts, ], evtree(factor(y) ~.,
                                             data = na.omit(data.frame(y, Z))))
plot(efit)


###3.2.- calssical tress (does not requiere listwise deletion)
#a.- grow large tree
#library(rpart)
set.seed(0203)
tree <- with(sdat[sdat$ts, ], rpart(factor(y) ~ . ,
                                    data = data.frame(y, Z), cp = -1, xval = 10,
                                    method = "class"))

plotcp(tree, col = "red")
#help(rpart)
#A good choice of cp for pruning is often 
#the leftmost value for which the mean lies below the horizontal line.

#b.- prune based on complexity
opcp <- tree$cptable[, "CP"][which.min(tree$cptable[, "xerror"])]
ptree <- prune(tree, cp = opcp)


#library(rpart.plot)
rpart.plot(ptree, roundint = FALSE)


png("output/anySUD.classical.tree.plot.png", width = 480 * 4, heigh = 480 * 4, res = 300)
rpart.plot(ptree, roundint = FALSE)
dev.off()

#===============================================================================
###.-inference (95% CI; p-values)
#===============================================================================

sdat$node.rtree <- factor(predict(as.party(ptree),
                                  type = "node", newdata = sdat$Z))
table(sdat$node.rtree)


#library(survival)
cfite <- clogit(y ~ node.rtree  + strata(id), data = sdat[!sdat$ts, ])
anova(cfite, update(cfite, . ~ 1 +  strata(id)))
summary(cfite)


#mean by group in training and validation data
t(with(sdat[sdat$ts,], tapply(y, node.rtree, mean)))
t(with(sdat[!sdat$ts,], tapply(y, node.rtree, mean)))
table(sdat$node.rtree)