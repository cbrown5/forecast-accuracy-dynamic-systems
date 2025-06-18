#Process forecast results and create plots
# 2024-09-02

rm(list = ls())
library(tidyverse)
library(mgcv)
library(patchwork)

load("Outputs/2025-01-10_inla-model-forecasts.rda")

dat2 <- read_csv("data/2025-01-24-ATRC-data-for-review.csv")
source("scripts/functions-sim-predictions.R")

small_number <- 1e-4 #for adding to zero values of MASE

theme_set(theme_classic())

spp_use2 <- c("Tosia australis", "Trachinops caudimaculatus", 
             "Jasus edwardsii",
             "Centrostephanus rodgersii")

#Calculate MASE

xall2 <- 
  bind_rows(xall, xall_rolling) %>%
  filter(ispp %in% spp_use) %>%
  filter(year > t0_test) %>%
  mutate(mase = abs(y2- med) / (ASE_scale + small_number),
         # calculate MASE
         horizon = year - t0_test, 
         mase2 = ifelse(mase > quantile(mase, 0.9, na.rm = TRUE), 
                        quantile(mase, 0.9, na.rm = TRUE), mase),
         log_mase2 = log10(mase2))


# xsum <- xall2 %>%
#   group_by(ispp, imodel, split_name, horizon, t0_test) %>%
#   summarize(mase = mean(mase, na.rm = TRUE),
#             n = n())



#select species
species <- unique(xall2$ispp)
g_horizon <- g_years <- g1sum <- g2sum <- NULL

for (ispp in species){
  print(ispp)
  # ispp <- species[2]
  spp_select <- ispp
  
  dat1 <- filter(xall2, ispp == spp_select) %>%
    mutate(model = factor(imodel),
           split_type = factor(split_name))
  
  dat1_sum <- dat1 %>%
    group_by(model, year, horizon, split_type) %>%
    summarize(mase = median(log(mase2)))
  
  g1 <- dat1 %>%
    filter(year %in% c(2012, 2015, 2020)) %>%
    ggplot() + 
    aes(x = horizon, y = mase, color = split_type) + 
    geom_point(alpha = 0.1) +
    stat_smooth() + 
    facet_grid(model~year, scales = "free_y") +
    labs(x = "Horizon (yrs)", y = "Median absolute squared error", 
         title = ispp) + 
    scale_y_log10()
  g1sum <- c(g1sum, list(g1))
  
  g2 <- dat1 %>%
    filter(horizon %in% c(1, 5, 8)) %>%
    ggplot() + 
    aes(x = year, y = mase, color = split_type) + 
    geom_point(alpha = 0.1) +
    stat_smooth() + 
    facet_grid(model~horizon, scales = "free_y") +
    labs(x = "Year", y = "Median absolute squared error", 
         title = ispp) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
    scale_y_log10()
  g2sum <- c(g2sum, list(g2))
  
}



wrap_plots(g1sum)

wrap_plots(g2sum)

#
# Summary plots 
#

xall3 <- xall2 %>%
  mutate(model = factor(imodel),
         split_type = factor(split_name),
         Period = case_when(
           year < 2012 ~ "Near-term",
           year >= 2012 & year < 2018 ~ "Mid-term",
           year >= 2015 ~ "Long-term"
         ),
         # rename species by their labels:
         spp_type = case_when(
           ispp == "Centrostephanus rodgersii" ~ "Increasing",
           ispp == "Tosia australis" ~ "Collapsing",
           ispp == "Trachinops caudimaculatus" ~ "Fluctuating",
           ispp == "Jasus edwardsii" ~ "Stable"
         )
  ) #%>%
# group_by(ispp, model, split_type, Period, horizon) %>%
# summarize(mase = mean(mase2, na.rm = TRUE))


pal <- rev(c("black", "darkblue", "hotpink"))

g1 <- xall3 %>%
  filter(model == "ar1") %>%
  filter(ispp %in% spp_use2) %>%
  mutate(split_period = paste(split_type, Period),
         split_type = factor(split_type, labels = c("Legacy", "Modern"))) %>%
  ggplot() + 
  aes(x = horizon, y = mase, group = split_period, color = Period,
      fill = Period,
      linetype = split_type) + 
  # geom_point(alpha = 0.1) +
  scale_linetype_manual("Splitting rule", values = c("solid", "dotted")) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = TRUE) + 
  facet_wrap(~spp_type, scales = "free_y", axes ="all_x") +
  scale_color_manual("Period", values = pal) +
  scale_fill_manual("Period", values = pal) +
  scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  labs(x = "Horizon (yrs)", y = "Median absolute scaled error") +
  # ylim(0, NA) + 
  coord_cartesian(ylim=c(0,3.1)) + 
  xlim(0, 16) 
g1

ggsave("figures/2025-10-01_fig4.png", g1, 
       width = 8, height = 6)

g2 <- xall3 %>%
  #set values > 3 to 3
  mutate(mase = ifelse(mase > 1000, 1000, mase)) %>%
  filter(ispp %in% spp_use2) %>%
  mutate(split_period = paste(split_type, Period),
         split_type = factor(split_type, labels = c("Legacy", "Modern"))) %>%
  ggplot() + 
  aes(x = horizon, y = mase, group = split_period, color = Period,
      fill = Period,
      linetype = split_type) + 
  # geom_point(alpha = 0.1) +
  scale_linetype_manual("Splitting rule", values = c("solid", "dotted")) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = TRUE) + 
  facet_grid(model~spp_type, scales = "free_y", axes ="all_x") +
  scale_color_manual("Period", values = pal) +
  scale_fill_manual("Period", values = pal) +
  scale_x_continuous(breaks = seq(0, 16, by = 2)) +
  labs(x = "Horizon (yrs)", y = "Median absolute scaled error") +
  ylim(0, NA) + 
  xlim(0, 16) 
g2
ggsave("figures/2025-10-01_all-spp-all-models-forecast-horizons.png", g2, 
       width = 10, height = 8)

