options(digits = 3)
options(scipen = 1000000)
options(warn=0)

library(tidyverse)
library(gt)
library(cobalt)
library(causalTree)
library(survey)

setwd("C:/Users/21983/OneDrive - ICF/ADIA/study 2/")

source("R/find.weights.r")

out <- "inj_victim"

#===============================================================================

dat <- read.csv("data/LS_all.csv")
dim(dat)
head(dat)
colnames(dat)

#' YE: As you said, in Study 2 we are interested in effect modification.
#' This suggest an alternative to the "case-control" design
#' (i.e., compare expousure to PC among kids exposed to ACEs
#' with [case] and without [control] an outcome)
#' We can instead use a "cohort" desing
#' (i.e., compare an outcome between those exposed and not to ACE
#' and check if the estimated effect vary by PC expousure

#===============================================================================
### . Outcome varaible
#===============================================================================


dat$y <- dat$inj_victim 
summary(dat$y)
dat <- dat[!is.na(dat$y), ] #n=864
#===============================================================================
### . Compute some varaibles
#===============================================================================
dat$cohort <- as.numeric(substr(dat$ChildDOB, 6, 9))
dat$cohort <- factor(dat$cohort)
summary(dat$cohort)

dat$anyACE <- as.numeric(apply(
  dat[, c("M", "MP", "MS", "MN", "MF", "ML", "ME", "MM", "MD")], 1, sum) >= 1)
table(dat$anyACE, useNA = "ifany")

#for the causal tree, one solution is to collapse certainly unexposed with exposure unknown
dat$anyACE[is.na(dat$anyACE)] <- 0

#===============================================================================
### . Find balancing weights among those with and without any ACEs
#===============================================================================
#Approximate Balancing Weights (these are an alternative to entropy balancing)
#Reference:  https://doi.org/10.48550/arXiv.2007.09056

# In orther to balance the missing pattern we need to
# for categorical variables, create an NA category (addNA)
# this should inlcude binary variables (even if not declared as factors)
# for continuous, add indicator is.na and impute mean

x <- c("center", "childrace_bsl", "cohort"
       , "childsex_bsl", "caregiver_married16", "hh_income16")

dat$X <-
  dat[, x] %>%
  mutate(across(where(is.factor), addNA, ifany = TRUE))  %>%
  mutate(across(
    where(~ is.numeric(.x) && any(is.na(.x)) && n_distinct(.x)==2)
    , ~addNA(factor(.x)))) %>%
  mutate(across(
    where(~ is.numeric(.x) && any(is.na(.x)))
    , is.na, .names = "NA_{.col}")) %>%
  mutate(across(where(is.numeric)
                , ~ replace(., is.na(.), mean(., na.rm = TRUE)))) %>%
  model.matrix(~., .) %>%
  .[, -1]


#weights add up to 1 in each group
tgt  <- colMeans(dat$X)
sd2  <- apply(dat$X, 2, \(v) if(mean(v) < .05) sqrt(.05 * .95) else sd(v))
phiX <- scale(dat$X, center = tgt, scale = sd2)
dat$w[dat$anyACE == 1] <- find.weights(phiX[dat$anyACE == 1, ], lambda = 0.01)
dat$w[dat$anyACE == 0] <- find.weights(phiX[dat$anyACE == 0, ], lambda = 0.01)
#' Ye: lambda=0 give maximumn similarity
#' at the cost of more varaibility of the weights
#' (i.e.,less effective sample size)
#' you coudl increase lambda a bit, e.g. 0.05,
#' to reduce weight  variability,
#' if you still get good balance

tapply(dat$y, dat$anyACE, mean)
tapply(dat$y * dat$w, dat$anyACE, sum)


with(dat, data.frame(X, anyACE, w)) %>%
  group_by(anyACE) %>%
  summarize(across(everything() & !w, list(
    mean,
    ~ weighted.mean(., w))
  )) %>%
  t

# Ye - how does this work with a continuous outcome?
bal.tab(dat$X, treat = dat$y, weights = dat$w,
        s.d.denom = "treated", un = TRUE, disp = c("means", "sds"))
love.plot(dat$X, treat = dat$anyACE, weights = dat$w,
          thresholds = c(m = .1), binary = "std", s.d.denom = "treated")

#===============================================================================
### . define a training sample
#===============================================================================
#stratfied random sample
dat <-
  dat %>%
  mutate(strata = paste0(x, anyACE)) %>%
  group_by(strata) %>%
  mutate(training.sample = if_else(id %in%
                                     sample(unique(id), round(n_distinct(id) / 2)), 1, 0)) %>%
  ungroup
table(dat$training.sample)
table(dat$strata, dat$training.sample)

file_name <- paste0("data/", out, ".LS.w.Rds")
saveRDS(dat, file = file_name)
#Ye: ALL the code above needs to be run once (it appplies to all outcomes).
# the code below, in contrast, vary by outome

#===============================================================================
### . find effect modifiers
#===============================================================================
dat <- readRDS(paste0("data/", out, ".LS.w.Rds"))

#ALL 28 PCE indicators
pf <- c('anysupadu_1',  	'anysupparent_1',  	'anysuprelative_1',  	'anysupnonfam_1',  	'fam_sat_1',  	'home_safety_1',  	
'prrelation_1',  	'neighborhood_exp_1',  	'school_safety_t_1',  	'srvc_use_1',  	'childcare_1',  	
'neighborhood_safety_2',  	'neighborhood_exp_2',  	'school_safety_y_2',  	'school_safety_t_2',  	
'srvc_use_2',  	'anysupadu_3',  	'anysupparent_3',  	'anysuprelative_3',  	'anysupnonfam_3',  	'prrelation_3',  	
'bestfriend_3',  	'socialpart_3',  	'parent_involv_3',  	'resid_stab_3',  	'neighborhood_safety_3',  	'neighborhood_exp_3',  	
'srvc_use_3')  


#PCE across time indicators
#pf <-  c("anysupadu", "anysupparent", "anysuprelative"
#         , "anysupnonfam", "fam_sat", "home_safety", "prrelation"
#         , "bestfriend", "socialpart", "parent_involv", "resid_stab"
#         , "neighborhood_safety", "neighborhood_exp"
#         , "school_safety_y", "srvc_use", "childcare")

summary(dat[, pf])
dat$PF <- dat[, pf]
# sdat$PF <- ifelse(runif(length(sdat$y ))>0.25,!sdat$y,sdat$y)
# table(sdat$PF, sdat$y)


###case weights
summary(dat$w)
tapply(dat$w, dat$anyACE, sum)
dat$cw <- NA
dat$cw[dat$anyACE == 1] <- dat$w[dat$anyACE == 1] * sum(dat$anyACE)
dat$cw[dat$anyACE == 0] <- dat$w[dat$anyACE == 0] * sum(!dat$anyACE)
tapply(dat$cw, dat$anyACE, sum)


df0 <-
  dat %>%
  filter(training.sample == 1)
dim(df0)
dim(dat)

# causal tree
# Reference https://doi.org/10.48550/arXiv.1504.01132

set.seed(0203)
tree <- with(df0, causalTree(factor(y) ~ .
                             , data = cbind(y, PF)
                             , treatment = anyACE
                             , weights = cw
                             , split.Rule =  "TOT" #' Try "CT", "fit" or "TOT"
                             , cv.option  =  "TOT" #' Try "CT", "fit" or "TOT"
                             , split.Honest = FALSE
                             , cv.Honest = FALSE
                             , minsize = 20 #manually fix to 20
)
)
plotcp(tree, col = "red")
tree$cptable

# tree$cptable
#> tree$cptable
#CP nsplit rel error xerror   xstd
#1 0.04824      0     1.000  1.000 0.0465
#2 0.01334      1     0.952  0.996 0.0463
#3 0.00484      2     0.938  0.996 0.0463
#4 0.00347      3     0.934  0.997 0.0464
#5 0.00000      4     0.930  0.998 0.0465


opcp <- tree$cptable[, 1][which.min(tree$cptable[, 4])]
ptree_causal <- prune(tree, cp = opcp)

#plot
rpart.plot(ptree_causal, roundint = FALSE)
plot(as.party(ptree_causal))
prp(ptree_causal, type = 4, # left and right split labels (see Figure 2)
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
    main = "Injury causal tree", # main title
    cex.main = 1.0, # use big text for main title
    branch.col = "gray", # color of branch lines
    branch.lwd = 2, # line width of branch lines
    roundint=FALSE) 

#save plot
png(paste0("output/", out, ".tree.plot.causal.png")
    , width = 480 * 4, heigh = 480 * 4, res = 300)
plot(as.party(ptree_causal))
rpart.plot(ptree_causal, roundint = FALSE)
prp(ptree_causal, type = 4, # left and right split labels (see Figure 2)
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
    main = "Injury causal tree", # main title
    cex.main = 1.0, # use big text for main title
    branch.col = "gray", # color of branch lines
    branch.lwd = 2, # line width of branch lines
    roundint=FALSE)
dev.off()

rules <- partykit:::.list.rules.party(as.party(ptree_causal))
dat$node.cau  <- factor(
  predict(as.party(ptree_causal),
          type = "node", newdata = dat, labels=rules))
table(dat$node.cau)

#===============================================================================
### . Inferece
#===============================================================================

###1.- select confirmation sample
df1 <- dat %>% filter(training.sample == 0)
dim(df1)
dim(dat)

###2.- define subgroups and check balance
table(df1$anyACE)
df1$ACE <- factor(df1$anyACE
                  , levels = c(0, 1),  labels = c("NoACE", "ACE"))

df1$ace_pf <- paste(df1$ACE, df1$node.cau, sep = "_")

df1$ace_pf <- factor(df1$ace_pf)

table(df1$ace_pf, useNA='ifany')
sapply(split(df1, df1$ace_pf), \(D) with(D, colMeans(D$X)))


###3.- Find weights that balance across all groups
tgt  <- colMeans(df1$X)
sd2  <- apply(df1$X, 2, \(v) if(mean(v) < .05) sqrt(.05 * .95) else sd(v))
df1$phiX <- scale(df1$X, center = tgt, scale = sd2)
grp <- levels(df1$ace_pf)
weights <-
  lapply(grp, \(g) {
    with(df1[!is.na(df1$ace_pf) & df1$ace_pf == g, ],
         data.frame(id, find.weights(phiX, lambda = 0.01))
    )}) %>%
  bind_rows

colnames(weights)[2] <- "new_w"
df1  <- left_join(df1, weights, by = "id")
df1$new_w[is.na(df1$new_w)] <- 0

###4.- Check balance before and after
tb_r <-
  with(df1, data.frame(X, ace_pf)) %>%
  subset(!is.na(ace_pf)) %>%
  group_by(grp = ace_pf) %>%
  summarize(across(all_of(colnames(df1$X))
                   , mean), n = n()) %>%
  pivot_longer(!grp, names_to = "var", values_to = "freq") %>%
  pivot_wider(names_from = "grp", values_from = "freq")

tb_w <-
  with(df1, data.frame(X, new_w, ace_pf)) %>%
  subset(!is.na(ace_pf)) %>%
  group_by(grp = ace_pf) %>%
  summarize(across(all_of(colnames(df1$X))
                   , ~ weighted.mean(., new_w)), n = n()) %>%
  pivot_longer(!grp, names_to = "var", values_to = "freq") %>%
  pivot_wider(names_from = "grp", values_from = "freq")

left_join(tb_r, tb_w, by = "var") %>% data.frame()

#saving the output in nice format
left_join(tb_r, tb_w, by = "var", suffix = c(".raw", ".weighted")) %>%
  gt %>%
  tab_header(title = "Before and after weighting") %>%
  gtsave(paste0("output/", out, "_raw_weighted.html"))

###5.- Test interaction
sdw    <- svydesign(id = ~ id, weights = ~ new_w, data = df1)


fit0   <- svyglm(as.factor(y) ~ anyACE
                 , design = sdw, family = binomial(link="cloglog"))
cbind(coef(fit0), confint(fit0))[2,]

fit1 <- svyglm(as.factor(y) ~ anyACE:node.cau + node.cau  - 1,
               design = sdw, family = binomial(link="cloglog"))
tb   <- cbind(coef(fit1), confint(fit1))
tb[grep("anyACETRUE:", rownames(tb)), ]

anova(fit0, fit1)
cbind(unique(dat$node.cau),
             coef(fit1)[grep("anyACE:", names(coef(fit1)))],
             confint(fit1)[grep("anyACE:", names(coef(fit1))), ]
       )[order(unique(dat$node.cau)), ]

#adding code for node labels
tb <- data.frame(unique(dat$node.cau),
                 coef(fit1)[grep("anyACE:", names(coef(fit1)))],
                 confint(fit1)[grep("anyACE:", names(coef(fit1))), ]
)[order(unique(dat$node.cau)), ]
colnames(tb) <- c("Node", "Est.", "95%CI_LL", "95%CI_UL")

#saving the output in nice format
tb %>%
  gt %>%
  tab_header(title = "Predicted values") %>%
  gtsave(paste0("output/predicted_cau_", out, ".html"))

#to find sample size
table(df1$node.cau[complete.cases(df1[,all.vars(formula(fit1))])])

###6.- Doubly robust test
#' BUT listwise deletion !!
fit0dr   <- svyglm(as.factor(y) ~ anyACE
                   + center + childrace_bsl + cohort
                   + childsex_bsl + caregiver_married16 + hh_income16
                   , design = sdw, family = binomial(link="cloglog"))
cbind(coef(fit0dr), confint(fit0dr))[2,]

fit1dr <- svyglm(as.factor(y) ~ anyACE:node.cau + node.cau  - 1
                 + center + childrace_bsl + cohort
                 + childsex_bsl + caregiver_married16 + hh_income16
                 , design = sdw, family = binomial(link="cloglog"))
tb   <- cbind(coef(fit1dr), confint(fit1dr))
tb[grep("anyACETRUE:", rownames(tb)), ]


anova(fit0dr, fit1dr)
