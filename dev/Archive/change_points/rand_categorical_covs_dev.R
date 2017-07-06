    #  Debug script for hare models
    #  05/2017
    #  Josh Nowak
################################################################################
    library(R2jags)
    library(readr)
    library(purrr)
    library(dplyr)
################################################################################
    #  Set working directory
    #setwd("C:/Users/josh.nowak/Documents/GitHub/SSH-camera-trap-analysis")
    setwd("/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/GitHub/SSH-camera-trap-analysis")
    #  Path to data
    #jjn <- "C:/Temp/NH Snowshoe Hare Photos_final.no.check.csv"
    jjn <- "/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/NH_hare_data.csv"
    #jjn <- "/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/SSH CO 2014-2016.csv"
    
    #  Source functions
    source("code/utility_functions.R")

################################################################################
    #  Load data
    rawd <- read_csv(
      jjn,
      col_types = "ccciiccccciicicc"
    )
    
################################################################################
    #  Morph raw data
    hares <- morph_data(rawd) %>%
      filter(
        Season == "Spring",
        Year == 2015
      )
    
################################################################################
    #  Call a single model step by step - mimics jags_call
    #  Set time_scale for the analysis
    #  Options are in the column names of hares, Month, Week, Julian
    time_scale <- "Julian"
    
    load.module("glm")
    
    #  Subset to days - to reduce redundancy and ease inits and data create
    days <- as.integer(unlist(hares[,time_scale]))
    first_day <- min(days)
    last_day <- max(days)
    
    #  Inits
    inits <- function(){
      list(
        #beta = runif(1, -100, 100),
        tau = runif(3, 0, 1000),
        start = runif(1, first_day, last_day)
      ) 
    }
    
    #  Gather data
    dat <- list(
      nobs = nrow(hares),
      day = days, 
      white = hares$White,
      ncam = max(hares$CameraNum),
      cam = hares$CameraNum,
      first_day = first_day,
      last_day = last_day,
      snow = hares$Snow,
      elev = hares$Elevation
    )
    
    # Parameters to monitor
    parms <- c(
      "alpha",
      #"beta", 
      "tau", 
      "start_pop", 
      "end_pop", 
      "tau_start", 
      "tau_end",
      #"start_cam", 
      #"end_cam"
      "snow_eff",
      "elev_eff"
    )

    #  Call jags
    out <- jags(
      data = dat, 
      inits = inits,
      parameters.to.save = parms,
      model.file = "models/rand_camera_categorical_covs.txt", 
      n.chains = 3,
      n.iter = 3000,
      n.burnin = 1000,
      n.thin = 3 
    )
################################################################################
    #  Check
    alphas <- out$BUGS$mean$alpha
    #b <- out$BUGS$mean$beta
    strt <- as.numeric(out$BUGS$mean$start_pop)
    ends <- as.numeric(out$BUGS$mean$end_pop)

    tbl <- tibble(
      cam = hares$CameraNum,
      day = hares$Julian,
      white = hares$White
    ) %>%
    arrange(cam, day, white) %>%
    mutate(
      J = 1 + (day > strt) + (day > ends),
      #  continuous model formulation
      #white_est = alphas[J] + b * (day - strt) * (J == 2)
      #  categorical model formulation
      white_est = alphas[J]
    )

################################################################################
    #  Summarise number of cameras per time period (=day)
    ss_summ <- tbl %>%
      group_by(day) %>%
      summarise(
        n_cam = length(unique(cam))
      )
    
    #  Setup plot
    plot(0, 0, 
         type = "n", 
         col = "red",
         ylim = c(0, 100),
         xlim = c(0, 220),
         xlab = "Time",
         ylab = "% White",
         bty = "l"
    )
    
    #  Add sample size to plot as "histogram" at the bottom
    lines(ss_summ$day, ss_summ$n_cam, type = "h", lwd = 5, col = "gray90")
    
    #  Add all raw data points if desired
    points(tbl$day, tbl$white, pch = 19, cex = 1.5, col = "gray90")
    
    #  Summarise raw camera data for comparison
    cam_summ <- tbl %>%
      group_by(cam, day) %>%
      summarise(
        white = mean(white),
        ss = n()
      ) %>%
      do(plt = points(.$day, .$white, pch = 19, col = "gray30"))

    points(tbl$day, tbl$white_est, pch = 3, col = "blue")
    abline(v = strt, lty = 2, col = "red")
    abline(v = ends, lty = 2, col = "red")

    print(select(tbl, J, white_est) %>% distinct())

    #  print model results
    print(out)
    
    ################################################################################
    #Diagnostics plots
    out.mcmc <- as.mcmc(out) # Convert model output into an MCMC object

    library(lattice)
    xyplot(out.mcmc, layout=c(3,3), aspect="fill") # chains history
    densityplot(out.mcmc, layout=c(3,3), aspect="fill") # posteriors
    autocorr.plot(out.mcmc) # autocorrelation plot
    gelman.plot(out.mcmc) 
    
    library(ggplot2)
    ggplot(data=tbl, aes()) +
      geom_jitter(aes(day,white_est,colour=as.factor(tbl$cam)), width=1, height = 2, alpha = .7) +
      theme(legend.position = "none") 
    
    ################################################################################
    #End
    