    #  Hare coat color change data exploration
    #  Josh Nowak/Marketa Zimova
    #  02/2017
################################################################################
    #  Packages
    require(readr)
    require(tidyr)
    require(dplyr)
################################################################################
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
    
    