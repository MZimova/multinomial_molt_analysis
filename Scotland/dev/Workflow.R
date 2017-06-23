     #  Set working directory
    setwd("/Users/marketzimova/Documents/WORK/DISSERTATION/2 Scotland/GitHub/ScottishHares") #MZ
    #setwd("/Users/...") #Josh

    #  Source functions
    source("helpers/Utility_funs.R")

################################################################################
    #  Load data and add temporal bits
    hare_dat <- readr::read_csv("/Users/marketzimova/Documents/WORK/DISSERTATION/2 Scotland/ANALYSIS and DATA/Scotland molt phenology data_averages.csv") 
    #hare_dat <- readr::read_csv("/Users/marketzimova/Documents/WORK/DISSERTATION/2 Scotland/ANALYSIS and DATA/Scotland molt phenology data_averages.csv") #Josh 
    morph_data(hare_dat)
################################################################################    
    
    # Consider using expand.grid and a join to accomplish this in long format.
    # Nested indexing may also be easier than a wide format.
    
    # The array is ragged at the moment and we want to predict the proportion of white:brown at each site whether it was visited or not. 
    # Like an encounter history for survival we want to make the data square.