# Utility functions for change point model for SSH camera trap study
#  Marketa Zimova/Josh Nowak
#  05/2017
# 1. Data transform
# 2. JAGS call
################################################################################


# transform data
morph.data <- function(x){
  out <- x %>%
    mutate( 
      Date = as.Date(Date, format="%m/%d/%y"),
      Julian = yday(Date),
      Year = year(Date),
      Month = month(Date),
      Week = as.numeric(format(Date, "%U")),
      fYear = factor(Year)
    ) 
  return(out)
}
###########

# next function
# call to JAGS
################################################################################
jags_call <- function(x, time_scale, ...){
  load.module("glm")
  
  #  Inits
  inits <- function(){
    list(
      beta = c(0,-1), 
      beta1 = c(0,-1),
      tau = 5,
      start = 100
    ) 
  }
  
  #  Gather data
  dat <- list(
    nobs = nrow(x),
    day = as.numeric(unlist(x[,time_scale])), 
    white = x$White 
  )
  
  # Parameters to monitor
  parms <- c("alpha", "beta", "beta1", "tau", "start", "end")
  
  #  Call
  out <- jags(
    data = dat,
    inits = inits,
    parameters.to.save = parms,
    ...
  )
  
  return(out)  
}
################################################################################
# End
