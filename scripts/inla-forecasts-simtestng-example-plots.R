# Forecast from INLA model by fixing hyperpars
# and making priors spikes. 
# Simulation testing
# 2024-11-30

rm(list = ls())
library(ggplot2)
library(dplyr)
library(INLA)
library(purrr)
library(patchwork)
library(scoringRules)
library(parallel)

theme_set(theme_classic())
source("scripts/score-functions.R")
source("scripts/inla-forecast-functions.R")


save_model_fits <- FALSE
prec.prior <- c(0.01, 0.01)
small_number <- 1e-4 #for adding to zero values of MASE
nsims <- 1 #number of simulation replicates for each scneario

# model_types <- c("rw1", "rw2")
model_types <- c("ar1")

#Forecast parameters
drop_year <- 32 #year of productivity change
years <- 1:55
n_years <- length(years)
t1_train <- 1 #first year of training
tn_train <- 25 #last year of training
t0_test_vals <- seq(tn_train, max(years)-1, by =5)#tn_train:(max(years)-1) #forecast origin
# t0_test_vals <- c(2006, 2015)
#stat to use for calibrating forecast model
central_stat_to_use <- "0\\.5quant" #mean, mode or "0\\.5quant"
#recomend 0.5 quant because of skew

#
# Gompertz process parameters
#
z_initial <- log(20)
r <- 0.25  # growth rate
k_initial <- log(20)   # carrying capacity
sigma <- 0.05  # process noise

k_drop <- c(0.25, 0.5, 0.75, 1)
spp_names <- paste0("Drop_", k_drop*100)

# ------------------------------------
# ******************************
# Data prep
# ******************************
# ------------------------------------


n <- length(years)


# ------------------------------------
# ******************************
# Loop over models, species and time periods to fit, calibrate and forecast
# ******************************
# ------------------------------------

xall_fixed <- NULL



#function to run one simulation for both split types
do_sim <- function(isim, ispp, imodel){
  #ispp <- spp_names[1]
  #imodel <- model_types[1]
  #isim <- 1
  # print(ispp)
  # print(imodel)
  cat("Rep proportion complete",isim/nsims)
  #
  # Simulate data for testing
  # # Simulate Gompertz process with Poisson observation errors
  
  # Initialize the process
  set.seed(isim)
  z <- numeric(n_years)
  z[1] <- z_initial  # initial population size in log scale
  k <- rep(k_initial, n_years)  # carrying capacity
  k[drop_year:n_years] <- k_initial*k_drop[idrop]  # carrying capacity
  theta <- (k-r)/k
  
  # Simulate the process
  for (t in 2:n_years) {
    z[t] <- r + theta[t]*z[t-1] + rnorm(1, 0, sigma)
  }
  # Convert to actual population size and add Poisson observation errors
  population_size <- exp(z)
  observed_data <- rpois(n_years, population_size)
  
  #plot(years, observed_data, type = "l", xlab = "Year", ylab = "Observed population size")
  #lines(years, population_size, col = "red")
  
  # Create a data frame for the simulated data
  dat2 <- data.frame(year = years)
  dat2[,ispp] <- observed_data
  dat2[, paste0(ispp, "_mean")] <- population_size

  
  #
  # Fit model - fixed training period 
  #
  
  xall_fixed <- fit_models_1site(ispp, imodel, dat2, t1_train, tn_train, t0_test_vals, 
                                 save_model_fits = FALSE, 
                                 split_name = "fixed", save_plots = TRUE, 
                                 central_stat_to_use = central_stat_to_use)   
  
  #
  # Fit model - rolling training period
  #
  xall_rolling <- NULL
  for (iorigin in t0_test_vals){
    t1_train_rolling <- t1_train +iorigin - min(t0_test_vals)
    print(t1_train_rolling:iorigin)
    xout <- fit_models_1site(ispp, imodel, dat2, t1_train_rolling, iorigin, iorigin, 
                             save_model_fits = FALSE,
                             split_name = "rolling", save_plots = TRUE, 
                             central_stat_to_use = central_stat_to_use)
    xall_rolling <- bind_rows(xall_rolling, xout)
  }
  
  xall2 <- 
    bind_rows(xall_fixed, xall_rolling) %>%
    filter(year > t0_test) %>%
    mutate(mase = abs(y2- med) / (ASE_scale + small_number),
           # calculate MASE
           horizon = year - t0_test, 
           mase2 = ifelse(mase > quantile(mase, 0.9, na.rm = TRUE), 
                          quantile(mase, 0.9, na.rm = TRUE), mase),
           log_mase2 = log10(mase2)) %>%
    select(-y)
  xall2[,ispp] <- NULL
  
  xall2$isim <- isim
  return(xall2)
}

#
# Loop over models and drops
#


#approx size for pre-allocation

datout <- data.frame(matrix(NA, nrow = 245*length(spp_names)*length(model_types)*nsims,
                            ncol = 15))
irow <- 1

for (ispp in spp_names[2]){
  idrop <- which(spp_names == ispp)
  for (imodel in model_types){
    print(idrop)
    print(imodel)

      xall2 <- mclapply(1:nsims, do_sim, ispp, imodel, mc.cores = 4)
    # system.time(xall2 <- lapply(1:nsims, do_sim, ispp, imodel))
    
    # for (isim in 1:nsims){
    
    xall2 <- do.call("rbind", xall2)
    #
    #save results
    #
    
    if(irow ==1) names(datout) <- names(xall2)
    datout[irow:(nrow(xall2)+irow-1),] <- xall2
    irow <- nrow(xall2)+irow
    
    rm(xall2)
    gc()
  }
}



