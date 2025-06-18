# Forecast from INLA model by fixing hyperpars (legacy split) vs rolling training period (modern split)
# 2024-08-30

rm(list = ls())
library(ggplot2)
library(dplyr)
library(INLA)
library(purrr)
library(patchwork)
library(scoringRules)
library(readr)

theme_set(theme_classic())
source("scripts/inla-forecast-functions.R")
dat2 <- read_csv("data/2025-01-24-ATRC-data-for-review.csv")

species <- names(dat2)[c(4:7)]

# datcovar <- read.csv("Outputs/2024-01-06_temperature-timeseries-with-lags.csv")
datcovar <- read.csv("data/2025-01-10_enviro_PC-timeseries-with-lags.csv")
datcovar2 <- datcovar %>%
  select(year, PC1, PC2, PC1_1yrlag, PC1_2yrlag, PC2_1yrlag, PC2_2yrlag)

save_model_fits <- TRUE
prec.prior <- c(0.01, 0.01)

spp_names <- c("Trachinops caudimaculatus", "Tosia australis",
                "Haliotis rubra", "Jasus edwardsii")

model_types <- c("ar1", "rw1", "rw2")
years <- min(dat2$year):max(dat2$year)

#Forecast parameters
t1_train <- 1992 #first year of training
tn_train <- 2006 #last year of training
t0_test_vals <- 2006:(max(years)-1) #forecast origin
# t0_test_vals <- c(2006, 2015)
#stat to use for calibrating forecast model
central_stat_to_use <- "0\\.5quant" #mean, mode or "0\\.5quant"
#recomend 0.5 quant because of skew

# ------------------------------------
# ******************************
# Data prep
# ******************************
# ------------------------------------

n <- length(years)
regions <- unique(dat2$region)
nregions <- length(unique(dat2$region))
nsites <- length(unique(dat2$isite))

# ------------------------------------
# ******************************
# Loop over models, species and time periods to fit, calibrate and forecast
# ******************************
# ------------------------------------

#
#Fixed training period 
#
xall <- NULL

for (ispp in spp_names){
  for (imodel in model_types){
    print(ispp)
    print(imodel)
    xout <- fit_models(ispp, imodel, dat2, t1_train, tn_train, t0_test_vals, 
                       save_model_fits = TRUE, 
                       split_name = "fixed", save_plots = TRUE, 
                       central_stat_to_use = central_stat_to_use,
                       savedate = "2025-01-10")   
    xall <- bind_rows(xall, xout)
  }
}
xall <- xall %>%
  select(-contains(spp_names))

#
#Rolling training period
#

xall_rolling <- NULL

for (ispp in spp_names){
  for (imodel in model_types){
    # ispp <- spp_names[1]
    # imodel <- model_types[1]
    print(ispp)
    print(imodel)
    for (iorigin in t0_test_vals){
      t1_train_rolling <- t1_train +iorigin - min(t0_test_vals)
      print(t1_train_rolling:iorigin)
      xout <- fit_models(ispp, imodel, dat2, t1_train_rolling, iorigin, iorigin, 
                         save_model_fits = TRUE,
                         split_name = "rolling", save_plots = TRUE, 
                         central_stat_to_use = central_stat_to_use,
                         savedate = "2025-01-10")
      xall_rolling <- bind_rows(xall_rolling, xout)
    }
  }
}

#find column names in xall_rolling that match spp_names and remove those columns
xall_rolling <- xall_rolling %>%
  select(-contains(spp_names))

#save model fits
save(xall, xall_rolling, file = "Outputs/2025-01-10_inla-model-forecasts.rda")
