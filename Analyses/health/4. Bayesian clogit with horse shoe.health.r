options(digits = 3)
options(scipen = 1000000)
library(tidyverse)
library(missRanger)

set.seed(0203)



#setwd("C:/Users/21983/OneDrive - ICF/ADIA/study 2")
outv <- "PH_Health"

#===============================================================================
###I. Imputing  protective factors
#===============================================================================
sdat <- readRDS(file = paste0("data/", outv, ".m.Rds"))
sdat <- sdat [order(sdat$id), ]

#.- Protective factors

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
cbind(
  apply(sdat[, zn], 2, min, na.rm = TRUE),
  apply(sdat[, zn], 2, max, na.rm = TRUE)
  )

##Ye: residential stability range form 0 to 3, it is the only one

#select varaibles
sdat$X <- sdat[, c(zn, "y", "id")]

#YE: this is only necesary if there are  factors
# in order to create dummies
#sdat$X <-
#  model.matrix(~ .
#  , model.frame(~ ., sdat[, zn], na.action = na.pass))[, -1] %>%
#  data.frame()
# head(sdat$X)



###multiple imputation
#5  for quick fits 100 for final run

filled <- replicate(20, {
    ix <- missRanger(data.frame(sdat$X)
    , num.trees = 50, pmm.k = 5)
    }, simplify = FALSE)



#===============================================================================
###II. Fitting Bayesian clogit with horsehose prior
#===============================================================================
library(posterior)
library(rstanarm)
options(mc.cores = parallel::detectCores())

post <-
  lapply(filled, \(x) {
    stan_clogit(y ~ . - id
     , strata = id
     , data = x
     , QR = TRUE
     , chains = 1, iter = 500 # for speed only
     , prior = hs()
     )
  })

post_all <-
  lapply(post, \(fit) {
    as.array(fit$stanfit) %>%
    as_draws_array
  }) %>%
  bind_draws(along = "chain")

tb <-
  summarize_draws(post_all
  , mean, ~quantile(.x, probs = c(0.05, 0.95))
  , default_convergence_measures()
  , .num_args = list(digits = 3, notation = "dec")) %>%
  print(n = 29)
tb

library(gt)
tb  %>%
  mutate_if(is.numeric, ~ as.numeric(.)) %>%
  gt %>%
  tab_header(title = "Difference in outcome by group") %>%
  fmt_number(2:4, 1:29, 3)  %>%
  fmt_number(5, 1:29, 1)  %>%
  fmt_number(6:7, 1:29, 0)  %>%
  gtsave(paste0("output/BVS_", outv, ".html"))
