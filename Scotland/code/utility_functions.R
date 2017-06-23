    #  Hare analysis utility functions
    #  Josh Nowak/Marketa Zimova
    #  02/2017
################################################################################
    morph_data <- function(x){
      #  A function to perform basic cleaning of column entries and colnames
      #  Takes a tibble/data.frame, very specific
      #  Returns a tibble

      out <- x %>%
        transmute(
          Area = as.character(Area),
          Date = as.Date(Date, "%m/%d/%Y"),
          Julian = as.numeric(format(Date, "%j")),
          Year = as.numeric(format(Date, "%Y")),
          Count = Count,
          Molt = Molt.index
        ) %>%
        rowwise() %>%
        mutate(
          Year = replace(Year, nchar(Year) == 2, paste0("20", Year)),
          Date = as.Date(paste0(Julian, "/", Year), "%j/%Y")
        ) %>%
        ungroup() %>%
        group_by(Area) %>%
        filter(n() > 50) %>%
        ungroup() %>%
        replace_na(
          D = 0,
          DD = 0,
          L = 0,
          LD = 0,
          LL = 0
        ) %>%
        spread(Molt, Count) %>%
        mutate(
          Total = D + DD + L + LD + LL
        )

    return(out)
    }
################################################################################
