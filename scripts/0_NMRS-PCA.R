#PCA of environmental time-series from NMRS data 
# 2025-01-10


library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(mgcv)

datvars <- read_csv("data/2025-01-24-NMRS-data-for-review.csv")
datsst <- read_csv("data/2025-01-24-SST-data-for-review.csv")

#make dataframe to keep predictions for each year in january
datpred <- data.frame(year = 1986:2024, month = 1)



#loop over each variable
for (var in c("Nitrate_umolL", "Silicate_umolL", "Salinity")){
  #filter the data
  dattemp <- datvars %>%
    select(date, year, days_since_start, month,SampleDepth_m, !!sym(var))
  dattemp$y <- pull(dattemp[,var])
  
  #fit the model
  fit <- gam(y ~ s(year, k = 30) + s(month, k = 4, bs = "cs"), 
    data = dattemp)
  #predict the model
  datpred[,var] <- predict(fit, newdata = datpred)
  
}

#
# Plot covariates
#

#join temperatue data
datc2 <- left_join(datpred, datsst, by = "year") %>%
  filter(year > 1988 & year < 2024)
#
# Make a PCA of the environmental covaraites
#

#pick the variables
datpca <- datc2 %>%
  select(Nitrate_umolL, Silicate_umolL, Salinity, mean_temp, max_temp, min_temp)
#convert to matrix
datpca <- as.matrix(datpca)
#scale the data
datpca <- scale(datpca)

#Plot the covarites (scaled)
datplot <- datpca %>%
  data.frame() 
datplot$Year <- 1989:2023
g1 <- datplot %>%
  tidyr::pivot_longer(-Year, names_to = "Variable", values_to = "Value") %>%
  ggplot() + 
  aes(x = Year, y = Value, color = Variable) + 
  geom_hline(yintercept = 0) +
  labs(y = "Value (scaled by S.D.)") +
  geom_line()
  g1
# ggsave(g1, file = "outputs/figures-methods-ms/physical-variables.png")

#perform PCA
pca <- prcomp(datpca, center = TRUE, scale = TRUE)
#plot the results
summary(pca)
#plot proportion of variance as barchart
vardat <- data.frame(var = pca$sdev^2/sum(pca$sdev^2), 
                     PC = 1:6)

#plot the components as different coloured lines in ggplot
datc2$PC1 <- pca$x[,1]
datc2$PC2 <- pca$x[,2]
datc2$PC3 <- pca$x[,3]

g2 <- ggplot(datc2) + 
  geom_line(aes(x = year, y = PC1, color = "PC1")) +
  geom_line(aes(x = year, y = PC2, color = "PC2")) +
  geom_line(aes(x = year, y = PC3, color = "PC3")) +
  theme_bw() +
  labs("Value")

g3 <- ggplot(vardat) + 
  aes(x = PC, y = var) + 
  geom_col() +
  labs(y = "Proportion of variance")

library(patchwork)
gall <- g1 + g2 + g3 + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(",
                  tag_suffix = ")", 
                  ) &
  theme(plot.tag = element_text(face = 'bold'))
gall
  
ggsave(gall, file = "figures/PCA-results.png", 
       width = 10, height = 3)

#make a biplot
biplot(pca)

#make a covariates file with the first two PCA components
# for each year and for 1 and 2 year lags
datcovar_1yrlag <- datc2 %>%
  select(year, PC1, PC2, PC3) %>%
  mutate(year = year + 1) %>%
  rename(PC1_1yrlag = PC1, PC2_1yrlag = PC2, PC3_1yrlag = PC3)
datcovar_2yrlag <- datc2 %>%
  select(year, PC1, PC2, PC3) %>%
  mutate(year = year + 2) %>%
  rename(PC1_2yrlag = PC1, PC2_2yrlag = PC2, PC3_2yrlag = PC3)

datcovar <- left_join(datc2, datcovar_1yrlag, by = "year") %>%
  left_join(datcovar_2yrlag, by = "year")

#save the file
write.csv(datcovar, "data/2025-01-10_enviro_PC-timeseries-with-lags.csv", row.names = FALSE)
