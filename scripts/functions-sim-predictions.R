#Functions to simulate predictions from a GAM
# 2023-03-24
# 
# based on code in ?gam.predict and  Wood, S.N. (2017) Generalized Additive Models: An
# Introduction with R (2nd edition). Chapman and
# Hall/CRC.

library(lazyeval)

rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}



simulate_gam_CIs <- function(model, 
                             newdata,
                             forms = list(~x), 
                             random_var = NULL,
                             offset = 0,
                             probs = c(0.025, 0.5, 0.975),
                             nsims = 1000,
                             func_to_apply = rep("quantile", length(forms))){
  
  #Function to simulate credible intervals for an mgcv GAM
  # 
  #Inputs:
  # model: gam object
  # newdata: data frame you want to predict too
  # forms: list of formulas to evaluate on the linear predictor
  # linear predictor is named 'x' internally. 
  # Can optionally be named list. 
  # Variables in the global environment can be
  # used in the formula
  # random_var: optional name variables you 
  # want to condition on, usually random effects
  # you want to ignore/set to zero
  # Can name multiple variables here as character vector. 
  # names must match variable names as they will appear
  # in the lpmatrix. Beware duplicate matches (using grepl here)
  # offset: numeric. offset to add back into predictions
  # probs: probability quantiles for CIs
  # func_to_apply: defaults to quantiles, other option is to sum
  # provide a character with same length of forms to 
  # specify summary function
  #
  #Output: returns a list of dataframes, were each dataframe
  # is prediction intervals for each formula. 
  
  nforms <- length(forms)
  
  #Predicting
  Xp <- predict(model, newdata = newdata, type="lpmatrix") 
  
  if (!is.null(random_var)){
    for (ivar in random_var){
      icol_rand <- grepl(ivar, colnames(Xp))
      Xp[,icol_rand] <- 0
    }
  }
  
  #simulate draws
  br <- rmvn(nsims, coef(model), model$Vp) 
  
  #Setup some vectors/matrices to store random draws of
  # predictions 
  
  effects <- lapply(1:nforms,
                    matrix,
                    ncol = nsims, 
                    nrow = nrow(newdata),
                    nsims)
  
  for (i in 1:1000){ 
    x <- Xp %*% br[i,] + offset ## replicate predictions
    #adding offset back in. This only matters for mean estimates
    #not important for the differences (sience the offset will cancel out)
    
    for (iforms in 1:nforms){
      effects[[iforms]][,i] <- f_eval(forms[[iforms]],
                                      data = list(x=x))
    }
    
  }
  
  # Compute CIs for each formula and cbind to original dataframe
  effects2 <- lapply(1:nforms, function(i){
    if (func_to_apply[i] == "quantile"){
      ztemp <- data.frame(t(apply(effects[[i]], 1, quantile, probs = probs)))  
    } else {
      
      ztemp <- data.frame(prob = apply(effects[[i]], 1, sum))/nsims  
    }
    
    ztemp <- setNames(ztemp, paste0(names(forms)[[i]], names(ztemp)))
    cbind(newdata, ztemp)
  })
  names(effects2) <- names(forms)
  return(effects2)
  
}


