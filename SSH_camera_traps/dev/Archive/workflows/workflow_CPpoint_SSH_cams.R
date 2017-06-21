    ###  Change point model for SSH camera trap photos in JAGS
    #  Calculate start and end for each spring in NH (later will expand to each 
    #   year site combination)
    #  Unmarked individuals, single season
    #  Marketa Zimova/Josh Nowak
    #  05/2017
################################################################################
    library(R2jags)
    library(readr)
    library(purrr)
    library(dplyr)
################################################################################
    #  Set working directory
    setwd("C:/Users/josh.nowak/Documents/GitHub/SSH-camera-trap-analysis")
    
    #  Path to data
    jjn <- "C:/Temp/NH Snowshoe Hare Photos_final.no.check.csv"

    #  Source functions
    source("code/utility_functions.R")

################################################################################
    #  Load data
    rawd <- read_csv(
      jjn,
      col_types = "icciicciccc"
    )
    
################################################################################
    #  Morph raw data
    hares <- morph_data(rawd)
    
################################################################################
    #  Call categorical model separately for each year by season combination
    fit_cat <- hares %>%
      group_by(Year, Season) %>%
      do(
        model = try(jags_call( 
          .,
          "Julian",
          model.file = "models/rand_camera_categorical.txt", 
          n.chains = 3,
          n.iter = 1000,
          n.burnin = 500,
          n.thin = 3 
        ))
      ) %>%
      ungroup
################################################################################
    #  Call continuous model separately for each year by season combination
    fit_cont <- hares %>%
      group_by(Year, Season) %>%
      do(
        model = try(jags_call( 
          .,
          "Julian",
          model.file = "models/rand_camera_continuous.txt", 
          n.chains = 3,
          n.iter = 1000,
          n.burnin = 500,
          n.thin = 3 
        ))
      ) %>%
      ungroup      
################################################################################
    #  Call all models in model folder over year and season, dumb idea given 
    #   the new models
    # mod2call <- list.files("./models", full.name = T)
    
    # fit <- map_df(mod2call, 
      # ~ hares %>%
          # group_by(Year, Season) %>%
          # do(
            # model = jags_call( 
              # .,
              # "Week",
              # model.file = .x, 
              # n.chains = 3,
              # n.iter = 1000,
              # n.burnin = 500,
              # n.thin = 3 
            # )
          # ) %>%
          # mutate(
            # model_file = .x
          # ) %>%
          # ungroup
      
    # )

################################################################################    
    #  Call models that aggregate years
    
    mod2call <- list.files("./models", pattern = "yr", full.name = T)
    
    fit_agg_yrs <- map_df(mod2call, 
      ~ hares %>%
          group_by(Season) %>%
          do(
            model = jags_call( 
              .,
              "Julian",
              model.file = .x, 
              n.chains = 3,
              n.iter = 1000,
              n.burnin = 500,
              n.thin = 3 
            )
          ) %>%
          mutate(
            model_file = .x
          ) %>%
          ungroup
      
    )

################################################################################

    #  Example print of a model
    fit_cat$model[[1]]

    #  Extracting mean of each model for muChange
    map_dbl(fit$model, ~.x$BUGS$mean$start)
    #  Extend the idea of extracting a mean to add a column
    mutate(fit, MeanStartDate = map_dbl(fit$model, ~.x$BUGS$mean$start))
################################################################################

    # Diagnostics
    # Simple plots
    #plot(fit$model[[1]]) # not sure what
    #traceplot(fit$model[[1]]) # history
    # convert model output into an MCMC object for more diagnostics
    fit.mcmc1 <- as.mcmc(fit$model[[1]]) 
    summary(fit.mcmc1)

    # library(coda)
    library(lattice)
    xyplot(fit.mcmc1, layout=c(3,3), aspect="fill") # chains history
    densityplot(fit.mcmc1, layout=c(3,3), aspect="fill") # posteriors
    #autocorr.plot(fit.mcmc1) # autocorrelation plot
    #gelman.plot(fit.mcmc1) 
    #geweke.diag(fit.mcmc1)

################################################################################
#  End
