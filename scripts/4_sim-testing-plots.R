#Plots for simulation testing results

library(dplyr)
library(ggplot2)
library(mgcv)

# Load the simulation results
load("Outputs/AR1-sim-results-2025-01-010.rda")
datout1 <- datout2
load("Outputs/RW1-RW2-sim-results-2025-01-010.rda")

# Bind the rows of the two datasets
datout_combined <- bind_rows(datout1, datout2) %>%
    #reorder ispp factor levels
        filter(t0_test %in% c(25, 35, 45)) %>%
    mutate(Scnr = 
    factor(ispp, 
        levels = c("Drop_25", "Drop_50", "Drop_75", "Drop_100"),
        labels = c("25% K", "50% K", "75% K", "100% K")),
        #relabel t0_test
        Epoch = factor(t0_test,labels = c("Before", "+3 yrs", "+12 yrs")),
        Split_Name = factor(split_name, labels = c("Legacy", "Modern"))
        )


# Create summary plots
datout_combined %>%
    mutate(grps = paste(split_name, imodel)) %>% 
    ggplot() + 
    #plot with horizon on x, mase2 on y, color by t0_test, fill by t0_test, group by split_period, 
    #linetype by split_name, facet grid by Scnr and imodel
    aes(x = horizon, y = mase2, color = imodel, 
        fill = imodel, group = grps,
        linetype = Split_Name) +
    #add a smooth line using a GAM model with 4 degrees of freedom
    stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = TRUE) +
    #add a manual scale for linetype
    scale_linetype_manual("Split method", values = c("solid", "dotted")) +
    #facet the plot by Scnr and imodel
    facet_grid(Scnr~Epoch) +
    #set fill and color scales
    scale_fill_brewer("Model", palette = "Set2") +
    scale_color_brewer("Model", palette = "Set2") +
    #set plot labels
    labs(x = "Horizon (yrs)", y = "Median absolute squared error") + 
    theme_bw()

#save teh plot
ggsave("figures/combined-sim-testing-plots.png", width = 6, height = 4)




# plot simpilified

datout_combined %>%
    mutate(grps = paste(split_name, imodel)) %>% 
    filter(imodel == "ar1" & Epoch == "+3 yrs" & Scnr == "50% K") %>%
    ggplot() + 
    #plot with horizon on x, mase2 on y, color by t0_test, fill by t0_test, group by split_period, 
    #linetype by split_name, facet grid by Scnr and imodel
    aes(x = horizon, y = mase2, 
         group = grps,
        linetype = Split_Name, color = Split_Name) +
    #add a smooth line using a GAM model with 4 degrees of freedom
    stat_smooth(method = "gam", formula = y ~ s(x, k = 4), se = TRUE) +
    #add a manual scale for linetype
    #facet the plot by Scnr and imodel
    #set fill and color scales
    #set plot labels
    labs(x = "Years from origin", y = "Error") + 
    theme_bw() + 
    scale_color_manual(values = c("red", "blue")) + 
    #remove axes labels, increase axes name size
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          #remove legend
          legend.position = "none")
ggsave("figures/forecast-horizon-plot.png", 
width = 4, height = 2)
