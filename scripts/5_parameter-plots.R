#Parameter plots and example plots 
# 2024-09-29

rm(list = ls())
library(ggplot2)
library(dplyr)
library(INLA)
library(purrr)
library(patchwork)
library(scoringRules)

theme_set(theme_classic())
source("scripts/inla-forecast-functions.R")
dat2 <- read_csv("data/2025-01-24-ATRC-data-for-review.csv")

species <- names(dat_merged)[c(5:49)]

datcovar <- read.csv("data/2025-01-10_enviro_PC-timeseries-with-lags.csv")

spp_names <- c("Trachinops caudimaculatus", "Tosia australis",
               "Jasus edwardsii", "Centrostephanus rodgersii")

# # # ------------------------------------
# ******************************
# Data prep
# ******************************
# ------------------------------------

years <- min(dat_merged$year):max(dat_merged$year)

n <- length(years)
regions <- unique(dat2$region)
nregions <- length(unique(regions))
nsites <- length(unique(dat2$isite))

#
# investiage one model
#

split_names <- c("rolling", "fixed")
imodel <- "ar1"
# t0_test_vals <- 2006:(max(years)-1) #forecast origin
t0_test_vals <- c(2006, 2010, 2015)
params_out <- NULL
params_out_fixed <- NULL
dat_plot <- NULL
for (ispp in spp_names){
  for (isplit in split_names){
    for (t0_test in t0_test_vals){
      # isplit <- "rolling"
      # ispp <- "Trachinops caudimaculatus"
      # t0_test <- 2006
      
      mtrain <- readRDS(
        paste0("Outputs/inla-forecasts/2025-01-10_inla-model-train_",
               isplit, "_", ispp, "_", imodel, "_", t0_test,".rds"))
      m_forecast <- readRDS(
        paste0("Outputs/inla-forecasts/2025-01-10_inla-model-forecast_",
               isplit, "_", ispp, "_", imodel, "_", t0_test, ".rds"))
      
      dattemp <- dat2 %>%
        mutate(y = get(ispp))
      years_test <- (t0_test+1):max(years) #years of testing
      
      dat_forecast <- dattemp %>%
        #set years in forecast period to NA
        mutate(y = ifelse(year %in% years_test, NA, y))
      
      # Extract hyperparameter estimates
      params_temp <- mtrain$summary.hyperpar
      params_temp$Split <- 
        #use title case
        paste0(substr(toupper(isplit), 1, 1), tolower(substr(isplit, 2, nchar(isplit))))
      params_temp$spp <- ispp
      params_temp$model <- imodel
      params_temp$t0_test <- t0_test
      params_temp$hyperpar <- rownames(params_temp)
      params_out <- rbind(params_out, params_temp)
      
      #extract fixed effect param estimates
      params_temp_fixed <- mtrain$summary.fixed
      params_temp_fixed$Split <- 
        #use title case
        paste0(substr(toupper(isplit), 1, 1), tolower(substr(isplit, 2, nchar(isplit))))
      params_temp_fixed$spp <- ispp
      params_temp_fixed$model <- imodel
      params_temp_fixed$t0_test <- t0_test
      params_temp_fixed$param <- rownames(params_temp_fixed)
      params_out_fixed <- rbind(params_out_fixed, params_temp_fixed)
      
      #
      # Compare forecasts 
      #
      
      #extract forecasts
      year_min <- min(mtrain$summary.random$year$ID)
      dattemp <- dat2 %>%
        filter(year >= year_min) %>%
        mutate(y = get(ispp)) %>%
        select(year, site_code, survey_id, isite, iregion, region, y)
      dattemp$median <- m_forecast$summary.fitted$`0.5quant`
      dattemp$lower <- m_forecast$summary.fitted$`0.025quant`
      dattemp$upper <- m_forecast$summary.fitted$`0.975quant`
      dattemp$Split <- paste0(substr(toupper(isplit), 1, 1), tolower(substr(isplit, 2, nchar(isplit))))
      dattemp$spp <- ispp
      dattemp$model <- imodel
      dattemp$t0_test <- t0_test
      
      dat_plot <- rbind(dat_plot, dattemp)
    }
  }
}

dp <- position_dodge(width = 1)

#plot rho
g1_list <- list()
g2_list <- list()
g_forecast <- list()
lw <- 1.2
ps <- 3

for (ispp in spp_names) {
    
  #
  #parameter plot
  #
  g1 <- params_out %>%
    filter(hyperpar == "Rho for year", spp == ispp) %>%
    mutate(Split = if_else(Split == "Fixed", "Legacy", "Modern")) %>%
    ggplot() +
    geom_point(aes(x = t0_test, y = `0.5quant`, color = Split), 
               position = dp, size = ps) +
    geom_linerange(
      aes(x = t0_test, ymin = `0.025quant`, ymax = `0.975quant`, color = Split),
      position = dp, linewidth = lw) +
    scale_x_continuous(breaks = t0_test_vals) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "", y = "Rho")
  
  g1_list[[ispp]] <- g1
  
  
  #plot precision
  g2 <- params_out %>%
    filter(hyperpar == c("Precision for year"), spp == ispp) %>% 
    mutate(Split = if_else(Split == "Fixed", "Legacy", "Modern")) %>%
    ggplot() +
    geom_point(aes(x = t0_test, y = `0.5quant`, color = Split), 
               position = dp, size = ps) +
    geom_linerange(
      aes(x = t0_test, ymin = `0.025quant`, ymax = `0.975quant`, color = Split),
      position = dp, linewidth = lw) +
    scale_x_continuous(breaks = t0_test_vals) +
    
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Forecast origin", y = "Precision")
  g2_list[[ispp]] <- g2
  
}

ghyper <- wrap_plots(g1_list) / wrap_plots(g2_list)+
  #gather legends
  plot_layout(guides = "collect", axis_titles = "collect_y") +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")",
                  tag_prefix = "(")
ghyper
# ggsave("Outputs/2025-01-10_hyperpar-plot.png", ghyper, 
      #  width = 7, height = 10)


#plot fixed effects
g3 <- params_out_fixed %>%
  filter(param == "(Intercept)") %>%
  ggplot() +
  geom_point(aes(x = t0_test, y = `0.5quant`, color = split), 
             position = dp) +
  geom_linerange(
    aes(x = t0_test, ymin = `0.025quant`, ymax = `0.975quant`, color = split),
    position = dp) +
  facet_wrap( ~ spp, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g4 <- params_out_fixed %>%
  filter(param == "PC1") %>%
  ggplot() +
  geom_point(aes(x = t0_test, y = `0.5quant`, color = Split), 
             position = dp) +
  geom_linerange(
    aes(x = t0_test, ymin = `0.025quant`, ymax = `0.975quant`, color = Split),
    position = dp) +
  facet_wrap( ~ spp, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g4


