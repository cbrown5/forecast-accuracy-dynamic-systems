#Function to fit and forecast using INLA

fit_models <- function(ispp, imodel, dat2, t1_train, tn_train, t0_test_vals, save_model_fits = FALSE, 
                       split_name = "none", save_plots = FALSE, central_stat_to_use = "0\\.5quant", 
                       savedate = "0000"){
  #ispp: species name
  #imodel: model name (ar1, rw1 or rw2)
  #dat2: data
  #t1_train: first year of training
  #tn_train: last year of training
  #t0_test_vals: forecast origin values
  #save_model_fits: logical, save model fits
  #set t0_test_vals == tn_train to do standard split test validation with rolling origin
  #set t0_test_vals > tn_train to do calibration then forecasting with fixed origin
  #split_name: name to use in filenames
  # save_plots: logical, save plots
  
  print(ispp)
  print(imodel)
  
  #check this species is in the data
  stopifnot(ispp %in% names(dat2))
  
  #check forecast origin is equal to or after last year of training
  stopifnot(min(t0_test_vals) >= tn_train)
  
  years_train <- t1_train:tn_train #years of training
  nt0_test_vals <- length(t0_test_vals)
  #Setup dataframes for fitting, calibration and forecasting
  dattemp <- dat2 %>%
    mutate(y = get(ispp))
  
  site_ase_scales <- dattemp %>%
    filter(!is.na(y)) %>%
    group_by(site_code, year) %>%
    arrange(year) %>%
    arrange(site_code) %>%
    summarize(y = mean(y)) %>%
    ungroup() %>%
    group_by(site_code) %>%
    reframe(ydiff = diff(y)) %>%
    group_by(site_code) %>%
    summarize(ASE_scale = mean(abs(ydiff)))
  
  dattemp <- dattemp %>%
    filter(year >= min(years_train))
  
  dattrain <- dattemp %>%
    filter(year %in% years_train)
  
  #
  # Fit the model to historical data 
  #
  form <- y ~ 1 +
    f(year, model = imodel, replicate = iregion) +
    f(isite, model = "iid") +
    f(PC1, model = "linear", mean.linear = 0, prec.linear = 1) +
    f(PC1_1yrlag, model = "linear", mean.linear = 0, prec.linear = 1) +
    f(PC1_2yrlag, model = "linear", mean.linear = 0, prec.linear = 1) 
  
  mtrain <- inla(form,
                 data = dattrain,
                 control.predictor = list(compute = TRUE, link = 1),
                 control.compute = list(return.marginals.predictor = TRUE,
                                        config = TRUE),
                 family = "poisson")
  
  #
  # Calibrate to calibration data and then forecast
  #
  dat_forecast_all <- NULL
  for (j in 1:nt0_test_vals){
    
    t0_test <- t0_test_vals[j]
    print(t0_test)
    #last year of trainign + 1 to forecast origin
    if(t0_test > tn_train){
      years_calibrate <- (tn_train+1):t0_test #last year of trainign + 1 to forecast origin
    } else {
      years_calibrate <- NULL
    }
    years_test <- (t0_test+1):max(years) #years of testing
    
    dat_forecast <- dattemp %>%
      #set years in forecast period to NA
      mutate(y = ifelse(year %in% years_test, NA, y))
    #
    #Get hyperparameters from mtrain 
    #
    
    #AR1 component
    hyper_names <- rownames(mtrain$summary.hyperpar)
    irow_prec_year <- grep("Precision for year", hyper_names)
    icol <- grep(central_stat_to_use, colnames(mtrain$summary.hyperpar))
    hyper_ar1 <- list(theta1 = #theta1 is precision, see ?inla.models
                        #no conversion neccessary. 
                        list(initial = log(mtrain$summary.hyperpar[irow_prec_year,icol]),
                             # mtrain$summary.hyperpar[irow_prec_year,icol],
                             fixed = TRUE))
    if (imodel == "ar1"){
      irow_rho <- grep("Rho", hyper_names)
      rho <- mtrain$summary.hyperpar[irow_rho,icol]
      hyper_ar1$theta2 <- #theta2 is related to rho, see ?inla.models for formula to convert between the two
        list(initial = log((1 + rho) / (1 - rho)),
             # rho, 
             fixed = TRUE)
    }
    # iid component
    irow_iid <- grep("isite", hyper_names)
    hyper_iid <- list(theta = list(initial = log(mtrain$summary.hyperpar[irow_iid,icol]),
                                   fixed = TRUE))
    
    #mean_temp component 
    row_fixed <- rownames(mtrain$summary.fixed)
    
    ifixed1 <- grep("^PC1$", row_fixed)
    ifixed2 <- grep("^PC1_1yrlag$", row_fixed)
    ifixed3 <- grep("^PC1_2yrlag$", row_fixed)
    
    form2 <- paste("y ~ 1 +", 
                   "f(year, model = imodel, replicate = iregion, hyper = hyper_ar1) +",
                   "f(isite, model = 'iid', hyper = hyper_iid) +",
                   "f(PC1, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed1, icol], prec.linear = 1E12) +",
                   "f(PC1_1yrlag, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed2, icol], prec.linear = 1E12) +",
                   "f(PC1_2yrlag, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed3, icol], prec.linear = 1E12)") %>%
      formula()
    
    # Forecast model 
    
    m_forecast <- inla(form2,family="poisson", data = dat_forecast, 
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute=list(config = TRUE),
                       control.fixed = list(mean.intercept = mtrain$summary.fixed[1,icol],
                                            prec.intercept = 1E12))
    
    #Plotting 
    #Make a plot of model forecasts
    dat_forecast$med <- m_forecast$summary.fitted.values$`0.5quant`
    dat_forecast$lwr <- m_forecast$summary.fitted.values$`0.025quant`
    dat_forecast$upr <- m_forecast$summary.fitted.values$`0.975quant`
    
    # Save results
    saveRDS(dat_forecast, 
            paste0("Outputs/inla-forecasts/",savedate,"_inla-forecast_",
                   split_name, "_", ispp, "_", imodel, "_", t0_test, ".rds"))
    
    if (save_model_fits){
      saveRDS(mtrain,
              paste0("Outputs/inla-forecasts/",savedate,"_inla-model-train_",
                     split_name, "_", ispp, "_", imodel, "_", t0_test,".rds"))
      saveRDS(m_forecast,
              paste0("Outputs/inla-forecasts/",savedate,"_inla-model-forecast_",
                     split_name, "_", ispp, "_", imodel, "_", t0_test, ".rds"))
    }
    if (save_plots){
      
      g1 <-  ggplot(dat_forecast) +
        geom_vline(xintercept = t0_test, color = "red", linetype = 2) +
        geom_vline(xintercept = tn_train, color = "black", linetype = 2) +
        geom_point(aes(x = year, y = y, color = site_code), alpha = 0.3) +
        geom_point(data = filter(dattemp, year > t0_test),
                   aes(x = year, y = y), color = "black",
                   alpha = 0.1) +
        geom_line(aes(x = year, y = med, color = site_code)) +
        geom_ribbon(aes(x = year, ymin = lwr, ymax = upr, fill = site_code),
                    alpha = 0.3, color = NA) +
        facet_wrap(~region, scales = "free") +
        scale_y_sqrt()
      ggsave(g1, 
             filename = paste0("Outputs/inla-forecasts/",savedate,"_inla-forecast-plot_",split_name, "_", ispp, "_", imodel, "_", t0_test, ".png"),
             width = 10, height = 6)  
    }
    dat_forecast <- select(dat_forecast, year, site_code, survey_id, 
    PC1, PC1_1yrlag, PC1_2yrlag,
                           ispp, 
                           y, med, lwr, upr, region) %>%
      mutate(y2 = .data[[ispp]]) %>%
      left_join(site_ase_scales, by = "site_code")
    
    dat_forecast$split_name <- split_name
    dat_forecast$ispp <- ispp
    dat_forecast$imodel <- imodel
    dat_forecast$t0_test <- t0_test
    dat_forecast_all <- rbind(dat_forecast_all, dat_forecast)
  }
  return(dat_forecast_all)
}



#Function to fit and forecast using INLA
#Fits to subset of sites
#This model has only a global year effect, not site specific
fit_models_sites <- function(ispp, imodel, dat2, sites_train, t1_train, 
                             tn_train, t0_test_val, save_model_fits = FALSE, 
                             split_name = "none", save_plots = FALSE, 
                             central_stat_to_use = "0\\.5quant"){
  #ispp: species name
  #imodel: model name (ar1, rw1 or rw2)
  #dat2: data
  # sites_train: sites to fit model to
  #t1_train: first year of training
  #tn_train: last year of training
  #t0_test_val: forecast origin value (only 1 allowed)
  #save_model_fits: logical, save model fits
  #set t0_test_val == tn_train to do standard split test validation with rolling origin
  #set t0_test_val > tn_train to do calibration then forecasting with fixed origin
  #split_name: name to use in filenames
  # save_plots: logical, save plots
  
  print(ispp)
  print(imodel)
  
  #check this species is in the data
  stopifnot(ispp %in% names(dat2))
  
  #check forecast origin is equal to or after last year of training
  stopifnot(min(t0_test_val) >= tn_train)
  
  years_train <- t1_train:tn_train #years of training
  #Setup dataframes for fitting, calibration and forecasting
  dattemp <- dat2 %>%
    mutate(y = get(ispp))
  
  site_ase_scales <- dattemp %>%
    filter(!is.na(y)) %>%
    group_by(site_code, year) %>%
    arrange(year) %>%
    arrange(site_code) %>%
    summarize(y = mean(y)) %>%
    ungroup() %>%
    group_by(site_code) %>%
    reframe(ydiff = diff(y)) %>%
    group_by(site_code) %>%
    summarize(ASE_scale = mean(abs(ydiff)))
  
  dattemp <- dattemp %>%
    filter(year >= min(years_train))
  
  dattrain <- dattemp %>%
    filter(year %in% years_train) %>%
    filter(isite %in% sites_train)
  
  #
  # Fit the model to historical data 
  #
  form <- y ~ 1 +
    f(year, model = imodel) +
    f(isite, model = "iid") +
    f(mean_temp, model = "linear", mean.linear = 0, prec.linear = 1) +
    f(mean_temp_1yrlag, model = "linear", mean.linear = 0, prec.linear = 1) +
    f(mean_temp_2yrlag, model = "linear", mean.linear = 0, prec.linear = 1) 
  
  mtrain <- inla(form,
                 data = dattrain,
                 control.predictor = list(compute = TRUE, link = 1),
                 control.compute = list(return.marginals.predictor = TRUE,
                                        config = TRUE),
                 family = "poisson")
  
  #
  # Calibrate to calibration data and then forecast
  #
  t0_test <- t0_test_val
  print(t0_test)
  #last year of trainign + 1 to forecast origin
  if(t0_test > tn_train){
    years_calibrate <- (tn_train+1):t0_test #last year of trainign + 1 to forecast origin
  } else {
    years_calibrate <- NULL
  }
  years_test <- (t0_test+1):max(years) #years of testing
  
  dat_forecast <- dattemp %>%
    #set years in forecast period to NA
    mutate(y = ifelse(year %in% years_test, NA, y))
  #
  #Get hyperparameters from mtrain 
  #
  
  #AR1 component
  hyper_names <- rownames(mtrain$summary.hyperpar)
  irow_prec_year <- grep("Precision for year", hyper_names)
  icol <- grep(central_stat_to_use, colnames(mtrain$summary.hyperpar))
  hyper_ar1 <- list(theta1 = #theta1 is precision, see ?inla.models
                      #no conversion neccessary. 
                      list(initial = log(mtrain$summary.hyperpar[irow_prec_year,icol]),
                           # mtrain$summary.hyperpar[irow_prec_year,icol],
                           fixed = TRUE))
  if (imodel == "ar1"){
    irow_rho <- grep("Rho", hyper_names)
    rho <- mtrain$summary.hyperpar[irow_rho,icol]
    hyper_ar1$theta2 <- #theta2 is related to rho, see ?inla.models for formula to convert between the two
      list(initial = log((1 + rho) / (1 - rho)),
           # rho, 
           fixed = TRUE)
  }
  # iid component
  irow_iid <- grep("isite", hyper_names)
  hyper_iid <- list(theta = list(initial = log(mtrain$summary.hyperpar[irow_iid,icol]),
                                 fixed = TRUE))
  
  #mean_temp component 
  row_fixed <- rownames(mtrain$summary.fixed)
  
  ifixed1 <- grep("^mean_temp$", row_fixed)
  ifixed2 <- grep("^mean_temp_1yrlag$", row_fixed)
  ifixed3 <- grep("^mean_temp_2yrlag$", row_fixed)
  
  form2 <- paste("y ~ 1 +", 
                 "f(year, model = imodel, hyper = hyper_ar1) +",
                 "f(isite, model = 'iid', hyper = hyper_iid) +",
                 "f(mean_temp, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed1, icol], prec.linear = 1E12) +",
                 "f(mean_temp_1yrlag, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed2, icol], prec.linear = 1E12) +",
                 "f(mean_temp_2yrlag, model = 'linear', mean.linear = mtrain$summary.fixed[ifixed3, icol], prec.linear = 1E12)") %>%
    formula()
  
  # Forecast model 
  
  m_forecast <- inla(form2,family="poisson", data = dat_forecast, 
                     control.predictor = list(compute = TRUE, link = 1),
                     control.compute=list(config = TRUE),
                     control.fixed = list(mean.intercept = mtrain$summary.fixed[1,icol],
                                          prec.intercept = 1E12))
  
  #Plotting 
  #Make a plot of model forecasts
  dat_forecast$med <- m_forecast$summary.fitted.values$`0.5quant`
  dat_forecast$lwr <- m_forecast$summary.fitted.values$`0.025quant`
  dat_forecast$upr <- m_forecast$summary.fitted.values$`0.975quant`
  
  dat_forecast <- select(dat_forecast, year, site_code, survey_id, 
                         mean_temp,mean_temp_1yrlag, mean_temp_2yrlag,
                         ispp, 
                         y, med, lwr, upr, region) %>%
    mutate(y2 = .data[[ispp]]) %>%
    left_join(site_ase_scales, by = "site_code")
  
  dat_forecast$split_name <- split_name
  dat_forecast$ispp <- ispp
  dat_forecast$imodel <- imodel
  dat_forecast$t0_test <- t0_test
  
  return(list(dat_forecast = dat_forecast, mtrain = mtrain, m_forecast = m_forecast))
  
}


fit_models_1site <- function(ispp, imodel, dat2, t1_train, tn_train, t0_test_vals, save_model_fits = FALSE, 
                       split_name = "none", save_plots = FALSE, central_stat_to_use = "0\\.5quant"){
  # Function to fit and forecast using INLA for a datsset with only one site
  #ispp: species name (used here to represent differnet levels of drop)
  #imodel: model name (ar1, rw1 or rw2)
  #dat2: data
  #t1_train: first year of training
  #tn_train: last year of training
  #t0_test_vals: forecast origin values
  #save_model_fits: logical, save model fits
  #set t0_test_vals == tn_train to do standard split test validation with rolling origin
  #set t0_test_vals > tn_train to do calibration then forecasting with fixed origin
  #split_name: name to use in filenames
  # save_plots: logical, save plots
  
  # print(ispp)
  # print(imodel)
  
  #check this species is in the data
  stopifnot(ispp %in% names(dat2))
  
  #check forecast origin is equal to or after last year of training
  stopifnot(min(t0_test_vals) >= tn_train)
  
  years_train <- t1_train:tn_train #years of training
  nt0_test_vals <- length(t0_test_vals)
  #Setup dataframes for fitting, calibration and forecasting
  dattemp <- dat2 %>%
    mutate(y = get(ispp))
  
  site_ase_scales <- dattemp %>%
    filter(!is.na(y)) %>%
    group_by(year) %>%
    arrange(year) %>%
    summarize(y = mean(y)) %>%
    ungroup() %>%
    reframe(ydiff = diff(y)) %>%
    summarize(ASE_scale = mean(abs(ydiff)))
  
  dattemp <- dattemp %>%
    filter(year >= min(years_train))
  
  dattrain <- dattemp %>%
    filter(year %in% years_train)
  
  #
  # Fit the model to historical data 
  #
  form <- y ~ 1 +
    f(year, model = imodel)
  
  mtrain <- inla(form,
                 data = dattrain,
                 control.predictor = list(compute = TRUE, link = 1),
                 control.compute = list(return.marginals.predictor = TRUE,
                                        config = TRUE),
                 family = "poisson")
  
  #
  # Calibrate to calibration data and then forecast
  #
  dat_forecast_all <- NULL
  for (j in 1:nt0_test_vals){
    
    t0_test <- t0_test_vals[j]
    # print(t0_test)
    #last year of trainign + 1 to forecast origin
    if(t0_test > tn_train){
      years_calibrate <- (tn_train+1):t0_test #last year of trainign + 1 to forecast origin
    } else {
      years_calibrate <- NULL
    }
    years_test <- (t0_test+1):max(years) #years of testing
    
    dat_forecast <- dattemp %>%
      #set years in forecast period to NA
      mutate(y = ifelse(year %in% years_test, NA, y))
    #
    #Get hyperparameters from mtrain 
    #
    
    #AR1 component
    hyper_names <- rownames(mtrain$summary.hyperpar)
    irow_prec_year <- grep("Precision for year", hyper_names)
    icol <- grep(central_stat_to_use, colnames(mtrain$summary.hyperpar))
    hyper_ar1 <- list(theta1 = #theta1 is precision, see ?inla.models
                        #no conversion neccessary. 
                        list(initial = log(mtrain$summary.hyperpar[irow_prec_year,icol]),
                             # mtrain$summary.hyperpar[irow_prec_year,icol],
                             fixed = TRUE))
    if (imodel == "ar1"){
      irow_rho <- grep("Rho", hyper_names)
      rho <- mtrain$summary.hyperpar[irow_rho,icol]
      hyper_ar1$theta2 <- #theta2 is related to rho, see ?inla.models for formula to convert between the two
        list(initial = log((1 + rho) / (1 - rho)),
             # rho, 
             fixed = TRUE)
    }
    
    form2 <- paste("y ~ 1 +", 
                   "f(year, model = imodel, hyper = hyper_ar1)") %>%
      formula()
    
    # Forecast model 
    
    m_forecast <- inla(form2,family="poisson", data = dat_forecast, 
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute=list(config = TRUE),
                       control.fixed = list(mean.intercept = mtrain$summary.fixed[1,icol],
                                            prec.intercept = 1E12))
    
    #Plotting 
    #Make a plot of model forecasts
    dat_forecast$med <- m_forecast$summary.fitted.values$`0.5quant`
    dat_forecast$lwr <- m_forecast$summary.fitted.values$`0.025quant`
    dat_forecast$upr <- m_forecast$summary.fitted.values$`0.975quant`
    
    # Save results
    saveRDS(dat_forecast, 
            paste0("Outputs/inla-sim-testing/2024-11-30_inla-forecast_",
                   split_name, "_", ispp, "_", imodel, "_", t0_test, ".rds"))
    
    if (save_model_fits){
      saveRDS(mtrain,
              paste0("Outputs/inla-sim-testing/2024-11-30_inla-model-train_",
                     split_name, "_", ispp, "_", imodel, "_", t0_test,".rds"))
      saveRDS(m_forecast,
              paste0("Outputs/inla-sim-testing/2024-11-30_inla-model-forecast_",
                     split_name, "_", ispp, "_", imodel, "_", t0_test, ".rds"))
    }
    if (save_plots){
      
      g1 <-  ggplot(dat_forecast) +
        geom_vline(xintercept = t0_test, color = "red", linetype = 2) +
        geom_vline(xintercept = tn_train, color = "black", linetype = 2) +
        geom_point(aes(x = year, y = y), alpha = 0.3) +
        geom_point(data = filter(dattemp, year > t0_test),
                   aes(x = year, y = y), color = "black",
                   alpha = 0.1) +
        geom_line(aes(x = year, y = med)) +
        geom_ribbon(aes(x = year, ymin = lwr, ymax = upr),
                    alpha = 0.3, color = NA) +
        scale_y_sqrt(limits = c(0, NA)) +
        labs(x = "Year", y = "Abundance")
      ggsave(g1, 
             filename = paste0("Outputs/inla-sim-testing/2024-11-30_inla-forecast-plot_",split_name, "_", ispp, "_", imodel, "_", t0_test, ".png"),
             width = 3, height = 2)  
    }
    dat_forecast <- select(dat_forecast, year, 
                           ispp, 
                           y, med, lwr, upr) %>%
      mutate(y2 = .data[[ispp]]) %>%
      mutate(ASE_scale = site_ase_scales$ASE_scale)
    
    dat_forecast$split_name <- split_name
    dat_forecast$ispp <- ispp
    dat_forecast$imodel <- imodel
    dat_forecast$t0_test <- t0_test
    dat_forecast_all <- rbind(dat_forecast_all, dat_forecast)
  }
  return(dat_forecast_all)
}