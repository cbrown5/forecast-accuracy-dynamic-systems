# Categories of species dynamics
#
# CJ Brown 
# 2025-03-22

#
rm(list = ls())
library(ggplot2)
library(dplyr)
library(INLA)
library(patchwork)
library(mgcv)
library(tidyr)

load(file = "Outputs/covariates/2024-08-31_dat2.rda")
spp_names <- c("Trachinops caudimaculatus", 
    "Tosia australis",
    "Jasus edwardsii",
    "Centrostephanus rodgersii")

    #Didn't end up using GAM as gets too complicated documenting decisions about
    # how to split long-term trend from
    # regional variation
form1 <- y ~ s(year, m=2, by = region) + 
    s(site_code, bs = 're')

dat2$site_code <- factor(dat2$site_code)
dat2$region <- factor(dat2$region)


fit_gam <- function(spp_name, dat2, form1){
    # spp_name <- "Tosia australis"
    # spp_name <- "Trachinops caudimaculatus"
    # spp_name <- "Jasus edwardsii"
    dat_temp <- dat2 %>%
        mutate(y = !!sym(spp_name))
    
    fit1 <- gam(form1, 
                data = dat_temp,
                family = "poisson")

    datlm <- dat_temp %>%
        group_by(site_code, region, year) %>%
        summarize(y = mean(y)) %>% group_by(region, year) %>%
        summarize(y = mean(y))

    fitlm <- glm(y ~ year, data = datlm)
    sfitlm <- summary(fitlm)
    sfit1 <- summary(fit1)
    
    datpred <- expand.grid(year = unique(dat_temp$year),
                            region = unique(dat_temp$region)) %>%
                mutate(site_code = dat_temp$site_code[1])
    datpred$fit <- predict(fit1, type = "response", newdata = datpred, 
        )
    x <- datpred %>%
        group_by(region) %>%
        summarize(CV = var(fit)/mean(fit)) %>%
        ungroup() %>%
        summarize(CV = mean(CV))
    
    
    return(
        data.frame(
            spp_name = spp_name,
            CV = x$CV, 
            # trend = exp(sfit1$p.coeff[2]),
            trend = coef(fitlm)[2],
            trend_norm = coef(fitlm)[2]/
                mean(dat_temp$y, na.rm = TRUE),
            pval = round(sfitlm$coefficients[2,4],3)
        )
    )
}

dat3 <- filter(dat2) #%>%
    # filter(year>2005)

gamfits <- lapply(spp_names, fit_gam, dat3, form1)

bind_rows(gamfits)

#
# Species trends plots 
#

dat_sum <- dat2 %>%
    select(year, site_code, Region = region, 
        `Trachinops caudimaculatus`,
        `Tosia australis`,
        `Jasus edwardsii`,
        `Centrostephanus rodgersii`
        ) %>%
    pivot_longer(-c(year:Region), names_to = "Species", values_to = "N") %>%
    group_by(site_code, Region, year, Species) %>%
    summarize(N = mean(N)) %>%
    group_by(Region, year, Species) %>%
    summarize(N = mean(N))

g1 <- ggplot(dat_sum) + 
    aes(x = year, y = N, color = Region, group = Region) +
    geom_line(color = "grey") + 
    stat_smooth(se = FALSE) + 
    scale_y_sqrt() +
    facet_wrap(~Species, scales = "free") +
    labs(x = "Year", 
            y = "Mean abundance per transect \n (sqrt scale)") + 
    scale_colour_manual(values = c("hotpink", "black", "darkblue", "turquoise3"))

ggsave(g1, 
    filename = "Outputs/figures-methods-ms/species-mean-trends.png",
    width = 6,
    height = 4)
