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


setwd('C:/Users/21983/OneDrive - ICF/ADIA/study 2') 
outv <- 'RB_deliq_hi' # RB_deliq and RB_deliq_Mod do not enough outcome for match control

#===============================================================================

dat <- read.csv('data/LS_all.csv')
head(dat)
colnames(dat)
summary(dat)
#===============================================================================
### I.- Developing sample
#===============================================================================
### 1.- We focus on population expericing ACE during childhood
# (and before the otucome was asssess)
dat$anyACE <- as.numeric(apply(
  dat[, c('M','MP','MS','MN','MF','ML','ME','MM','MD')], 1, sum) >= 1)
table(dat$anyACE) #0=383 1=971
dat <- dat[dat$anyACE==1,] #n=971

### 2.- In that population, we want to campare those who develop a certain outcome to those who do not
# ()'case-control' design)
#table(dat$RB_deliq, useNA='ifany') #0=317 1=318 na=336
#table(dat$RB_deliq_Mod, useNA='ifany') #0=413 1=222 na=336
table(dat$RB_deliq_hi, useNA='ifany') #0=513 1=121 na=337
dat$y  <- dat$RB_deliq_hi  
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

#===============================================================================
###II. Exploring possible protective factors
#===============================================================================
sdat <- readRDS(sdat, file=paste0('data/', outv, '.m.Rds'))

#1. define a training sample
length(unique(sdat$id))
set.seed(0203)
ts <- sample(unique(sdat$id), length(unique(sdat$id)) / 2)
sdat$ts <- as.numeric(sdat$id %in% ts)
table(sdat$ts, useNA='ifany')

#2.- Protective factors
dim(sdat)
summary(sdat)
#ALL PCE indicators
z.n  <- c('anysupadu_1',  	'anysupparent_1',  	'anysuprelative_1',  	'anysupnonfam_1',  	'fam_sat_1',  	'home_safety_1',  	
  'prrelation_1',  	'neighborhood_exp_1',  	'school_safety_t_1',  	'srvc_use_1',  	'childcare_1',  	
  'neighborhood_safety_2',  	'neighborhood_exp_2',  	'school_safety_y_2',  	'school_safety_t_2',  	
  'srvc_use_2',  	'anysupadu_3',  	'anysupparent_3',  	'anysuprelative_3',  	'anysupnonfam_3',  	'prrelation_3',  	
  'bestfriend_3',  	'socialpart_3',  	'parent_involv_3',  	'resid_stab_3',  	'neighborhood_safety_3',  	'neighborhood_exp_3',  	
  'srvc_use_3')  

#PCE across time indicators
#z.n  <- c('anysupadu', 'anysupparent', 'anysuprelative', 'anysupnonfam', 'fam_sat', 'home_safety', 'prrelation',
#          'bestfriend', 'socialpart', 'parent_involv', 'resid_stab', 'neighborhood_safety', 'neighborhood_exp',
#          'school_safety_y', 'srvc_use', 'childcare')

summary(sdat[, z.n])


###2.2- centering
# we need to center covariates at the pair-mean
# Stanfill et al. (2019) https://doi.org/10.1177/1179597219858954

#'old school' way
#sdat   <-  sdat[order(sdat$id),]
#sdat$Z <-  model.matrix(~ ., model.frame(~ ., sdat[,z.n], na.action=na.pass))[,-1]
#sdat$Z <-  do.call(rbind,lapply(split(data.frame(Z),sdat$id), function(d) apply(d,2,scale,scale=F) ))

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
df0 <-
  sdat %>%
  filter(ts == 1)
dim(df0)
dim(sdat)
set.seed(0203)
#3.1. - evolutionary tree (just for comparisson)

efit <- with(df0, evtree(factor(y) ~.,
                                             data = na.omit(data.frame(y, Z))))
plot(efit)


###3.2.- calssical tress (does not requiere listwise deletion)
#a.- grow large tree
#library(rpart)

tree <- with(df0, rpart(factor(y) ~ . ,
                                    data = data.frame(y, Z), cp = -1, xval = 10,
                                    method = "class"
                                    ,control=list(minbucket = 15)
                        )) #adjust to avoid produce too many groups

plotcp(tree, col = "red")
tree$cptable
#> tree$cptable
#CP nsplit rel error xerror   xstd
#1  0.1333      0     1.000  1.300 0.0871
#2  0.0500      1     0.867  1.100 0.0908
#3  0.0333      3     0.767  0.983 0.0913
#4 -1.0000      4     0.733  0.983 0.0913

#A good choice of cp for pruning is often 
#the leftmost value for which the mean lies below the horizontal line.

#b.- prune based on complexity

opcp <- tree$cptable[, "CP"][which.min(tree$cptable[, "xerror"])]

#alternative pruneing
#opcp <- 0.983 + 1*0.0913 # 1 SE rule
ptree <- prune(tree, cp = opcp)

rpart.plot(ptree, roundint = FALSE)
plot(as.party(ptree))
prp(ptree, type = 4, # left and right split labels (see Figure 2)
    clip.right.labs = FALSE, # full right split labels
    extra = 101, # show nbr of obs and percentages (see Figure 3)
    under = TRUE, # position extra info _under_ the boxes
    under.cex = 1, # size of text under the boxes (default is .8)
    fallen.leaves = TRUE, # put leaves at the bottom of plot
    box.palette = "GnYlRd", # color of the boxes
    branch = .3, # branch lines with narrow shoulders and down slopes
    round = 0, # no rounding of node corners i.e. use rectangles
    leaf.round = 9, # round leaf nodes (for leaves, this supersedes the round arg)
    #prefix = "Depression Score\n", # prepend this string to the node labels
    main = "High deliquent behavior classical tree", # main title
    cex.main = 1.0, # use big text for main title
    branch.col = "gray", # color of branch lines
    branch.lwd = 2, # line width of branch lines
    roundint=FALSE)


png(paste0("output/", outv, ".tree.plot.classical.png")
    , width = 480 * 4, heigh = 480 * 4, res = 300)
rpart.plot(ptree, roundint = FALSE)
plot(as.party(ptree))
prp(ptree, type = 4, # left and right split labels (see Figure 2)
    clip.right.labs = FALSE, # full right split labels
    extra = 101, # show nbr of obs and percentages (see Figure 3)
    under = TRUE, # position extra info _under_ the boxes
    under.cex = 1, # size of text under the boxes (default is .8)
    fallen.leaves = TRUE, # put leaves at the bottom of plot
    box.palette = "GnYlRd", # color of the boxes
    branch = .3, # branch lines with narrow shoulders and down slopes
    round = 0, # no rounding of node corners i.e. use rectangles
    leaf.round = 9, # round leaf nodes (for leaves, this supersedes the round arg)
    #prefix = "Depression Score\n", # prepend this string to the node labels
    main = "High deliquent behavior classical tree", # main title
    cex.main = 1.0, # use big text for main title
    branch.col = "gray", # color of branch lines
    branch.lwd = 2, # line width of branch lines
    roundint=FALSE)
dev.off()


rules <- partykit:::.list.rules.party(as.party(ptree))
sdat$node.rtree <- factor(predict(as.party(ptree),
                                  type = "node", newdata = sdat$Z), labels=rules)
table(sdat$node.rtree)

#===============================================================================
###.-inference (95% CI; p-values)
#===============================================================================
df1 <-
  sdat %>%
  filter(ts == 0)
dim(df1)
dim(sdat)


#library(survival)
cfite <- clogit(y ~ node.rtree  + strata(id), data = df1)
anova(cfite, update(cfite, . ~ 1 +  strata(id)))
summary(cfite)


#mean by group in training and validation data
t(with(sdat[sdat$ts==0,], tapply(y, node.rtree, mean)))
t(with(sdat[sdat$ts==1,], tapply(y, node.rtree, mean)))
table(sdat$node.rtree)
table(df1$node.rtree)