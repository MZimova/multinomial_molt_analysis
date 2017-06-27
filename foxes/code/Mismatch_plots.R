# Plots for fox mismatch
# 6/18/17
# Marketa Zimova

library(ggplot2)
library(gridExtra)
library(lubridate)
library(dplyr)

# load data
setwd('/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data')
fox <-read.csv("molts5.3.csv", header=T,sep=",")

# prepare data
fox$Date <- mdy(fox$Date)
fox$Julian <- yday(fox$Date)
fox$Year <- year(fox$Date)
fox$Week <- week(fox$Date)
fox$fYear <- factor(fox$Year)
fox$Mismatch <- fox$Color- fox$Snow

# Calculate weekly mean mismatch and SD
fox_df <- tbl_df(fox)
weekly <- group_by(fox_df,Year,Week, Morph, area)

means <- summarize(weekly, mm = mean(Mismatch), mmcount = n(),sd= sd(Mismatch))
means$upper <- means$mm + means$sd
means$lower <- means$mm - means$sd

# Weekly means and SDs
ggplot(data=means, aes(Week, mm)) +
  facet_grid(area ~ Year) +
  geom_line(data=means, aes(Week, mm, group=Morph:area, colour=Morph:area)) +
  geom_line(data=means, aes(Week, upper, group=Morph:area, colour=Morph:area),linetype=5, alpha=.7) +
  geom_line(data=means, aes(Week, lower, group=Morph:area, colour=Morph:area),linetype=5, alpha=.7) +
  geom_jitter(data=fox, aes(Week, Mismatch, colour=Morph:area), alpha=.3,size=.5) +
  labs(title="Average weekly mismatch")  
  
<<<<<<< HEAD:foxes/code/Mismatch_plots.R

=======
################## boxplots (but something's wrong)
H<- subset(fox,  area %in% c('helags'))
V<- subset(fox,  area %in% c('vindelfjallen'))

Helags <- ggplot(data=H, aes(Week, Mismatch)) +
  geom_boxplot(data=H, aes(Julian, Mismatch, group=Week,color=Morph)) +
         labs(title="Helags") 

Vindel <- ggplot(data=V, aes(Julian, Mismatch)) +
  geom_boxplot(data=V, aes(Julian, Mismatch, group=Week,color=Morph)) +
  labs(title="Vindel") 

grid.arrange(Helags, Vindel, nrow=2, ncol=1)

  
>>>>>>> origin/master:foxes/code/Mismatch_plots
