    # Utility functions for change point model for SSH camera trap study
    #  Marketa Zimova/Josh Nowak
    #  05/2017
################################################################################
    morph_data <- function(x){
      #  A function to format the data for analysis
      #  Takes raw data as returned by read_csv 
      #  Returns a data.frame
      
      seasons <- c("Fall", "Spring")
      
      out <- x %>%
        mutate(
          Date = as.Date(Date, format = "%m/%d/%y"),
          Julian = as.integer(format(Date, "%j")),
          Year = as.integer(format(Date, "%Y")),
          Month = as.integer(format(Date, "%m")),
          Week = as.integer(format(Date, "%U")) + 1,
          CameraNum = as.integer(as.factor(Camera)),
          Area = as.factor(area),
          Morph = as.factor(Morph),
          Season = seasons[(Month %in% 1:7) + 1]
        )
        
      #  Check data
      check_data(out)
        
    return(out)
    }
################################################################################
    check_data <- function(x){
      #  A function to check a few columns in the data
      #  Takes a tibble/df with column names Date, White, Snow, Camera
      #  Returns a tibble
      #  Intended to be called within morph_data
      
      #  Date check
      dates <- x[,grepl("date", colnames(x), ignore.case = T)][[1]]
      dt_na <- sum(is.na(dates))
      
      stopifnot(class(dates) == "Date")
      
      if(dt_na > 0){
        cat(
          "\n\n",
          "Date conversion failed in row(s):",
          "\n",
          which(is.na(dates)),
          "\n\n"
        )
      }
      
      #  Print date range
      cat(
        "\n",
        "Date range:",
        as.character(min(dates)), "/", as.character(max(dates)),
        "\n\n"
      )

      #  White check
      wht <- unlist(x[,grepl("white", colnames(x), ignore.case = T)])
      
      stopifnot(class(wht) == "integer")
      
      if(any(wht < 0 || wht > 100)){
        cat(
          "\n\n",
          "Verify value of White in row(s):",
          "\n",
          which(wht < 0 || wht > 100),
          "\n\n"
        )
      }
      
      #  Snow check
      snw <- unlist(x[,grepl("snow$", colnames(x), ignore.case = T)])
      
      stopifnot(class(snw) == "integer")
      
      if(any(snw < 0 || snw > 100)){
        cat(
          "\n\n",
          "Verify value of Snow in row(s):",
          "\n",
          which(snw < 0 || snw > 100),
          "\n\n"
        )
      }
      
      #  Camera check - anticipates a numeric identifier for the camera
      cam <- unlist(x[,grepl("cameranum", colnames(x), ignore.case = T)])
      
      stopifnot(class(cam) == "integer")
      
      if(any(!is.finite(cam))){
        cat(
          "\n\n",
          "Verify value of Camera in row(s):",
          "\n",
          which(!is.finite(cam)),
          "\n\n"
        )
      }
    
    return()
    }

################################################################################
    #jags_call <- function(x, time_scale, ...){
      #  A wrapper to call jags on grouped data dplyr style
      #  Takes output of morph_data, time_scale (a column name defining grouping
      #   variable) and other arguments to pass to jags

      #load.module("glm")
      
      #  Subset to days - to reduce redundancy and ease inits and data create
      #days <- as.integer(unlist(x[,time_scale]))
      #first_day <- min(days)
      #last_day <- max(days)
      
      # #  Inits
      # inits <- function(){
      #   list(
      #     #beta = runif(1, -100, 100),
      #     tau = runif(3, 0, 10)
      #   ) 
      # }
      
      # #  Gather data
      # dat <- list(
      #   nobs = nrow(x),
      #   day = days, 
      #   white = x$White,
      #   ncam = max(x$CameraNum),
      #   cam = x$CameraNum,
      #   first_day = first_day,
      #   last_day = last_day,
      #   yr = x$Year - min(x$Year) + 1,
      #   nyr = length(unique(x$Year))
      # )
      # 
      # # Parameters to monitor
      # parms <- c(
      #   "alpha", 
      #   "beta", 
      #   "tau", 
      #   "start",
      #   "end",
      #   "start_pop", 
      #   "end_pop", 
      #   "tau_start", 
      #   "tau_end"#,
      #   #"start_cam", 
      #   #"end_cam"
      # )
      # 
      #  Call
    #   out <- jags(
    #     data = dat,
    #     inits = inits,
    #     parameters.to.save = parms,
    #     ...
    #   )
    #   
    # return(out)  
    # }
################################################################################
# End
