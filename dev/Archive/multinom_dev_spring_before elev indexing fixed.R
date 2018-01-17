#  Multinomial analysis to calculate molt start and end dates 
     # (and covariates effects) in SHH
# 3 categories of molt
# Random effects on camera traps
#  06/2017
#  Josh Nowak
################################################################################
library(R2jags)
library(readr)
library(purrr)
library(dplyr)
library(beepr)
library(mcmcplots)
library(ggplot2)
############################################################################
  # Load hare data  
setwd('//Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study')
# run to merge 3 sites hare data or run manually or load saved one
  # source('analysis and code/SSH_data prep.R')
# load from RData
load('data/hares_daymet_present_cl.RData')
hares_daymet <-hares_daymet_cl

# or to check whether better with clusters
#load('data/finalhares_3sites_WhiteNoNAs_1kmclusters.Rdata')

# all sites spring
hares <- subset(hares_daymet, Julian >50 & Julian <190 & Site=="Can")
# hares <- subset(hares_1kmcluster, Julian >50 & Julian <190 & Site=="Can")

ggplot(hares, aes(Julian, White,colour=as.factor(Site))) + geom_jitter(width = 0.05, height = 0.05)

################################################################################
#  Call a single model step by step - mimics jags_call
#  Set time_scale for the analysis (pptions are in the column names of hares, Month, Week, Julian)
setwd("/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis")
time_scale <- "Julian"
load.module("glm")

#  Subset to days - to reduce redundancy and ease inits and data create
days <- as.integer(unlist(hares[,time_scale]))
first_day <- min(days)
last_day <- max(days)

#  Create categorical response
response <- cut(hares$White3, 3, labels = 1:3)

# made it ordered
#response <- ordered(response)

#  Inits
inits <- function(){
  list(
    alpha = rnorm(3)
  )
}

#  Gather data
dat <- list(
  nobs = nrow(hares),
  day = days, 
  cam = as.numeric(as.factor(hares$CameraNum)), #change to camera
  y = response,
  nbins = 3,
  ndays = last_day,
  ncam = length(unique(hares$CameraNum)), #change to camera
  elev = as.numeric(hares$Elevation)
  #lat = as.numeric(hares$Lat)
  #season_tavg = as.numeric(hares$cl_Pres.spring_tmin)
)

# Parameters to monitor
parms <- c("beta", "alpha","pp","sigma_cam","tau_cam", "d_elev_eff" #"season_tavg_eff"
           #"elev_eff", "lat_eff"#,
           #"sigma", "rho", p_rand","cat_mu" 
)

#  Call jags
start.time <- Sys.time()
out <- jags.parallel(
  data = dat, 
  inits = NULL,
  parameters.to.save = parms,
  model.file = "models/multinom_covs_inter.txt", #multinom_covs
  n.chains = 3,
  n.iter = 1000,
  n.burnin = 100,
  n.thin = 3
)
end.time <- Sys.time();(time.taken <-end.time-start.time)
beep()

################################################################################
# Save results out
#save(out, file = "/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis/out_CO_lat_elev.RData")
#load("out_CO_yr_all.RData")

# Save results as csv
#writes csv with results
out.sum <- out$BUGS$summary 
#write.table(out.sum, file="/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis/results/test_1thrubin.csv",sep=",")

options(max.print=100) #extend maximum for print
print(out)
#out$BUGS$mean$lat_eff

#not avg, maybe min has negative effect on 2= p of changing as opposed to being brown decreases as temperature is warmer (in warmer temp youre more likely to be brown than changing)
################################################################################
# Plots
#  Find start dates
starts <- apply(out$BUGS$sims.list$pp[,3,], 1, function(x){ 
  min(which(x < 0.9)) 
})
hist(starts, xlab = "Day")
quantile(starts, c(0.025, 0.5, 0.975))

#  Find end dates
ends <- apply(out$BUGS$sims.list$pp[,1,], 1, function(x){ 
  min(which(x > 0.9)) 
})
hist(ends, xlab = "Day")
quantile(ends, c(0.025, 0.5, 0.975))

#  Find mid-points
mids <- apply(out$BUGS$sims.list$pp[,2,], 1, function(x){ 
  min(which(x == max(x)))
})
hist(mids, xlab = "Day")
quantile(mids, c(0.025, 0.5, 0.975))

#Plot start and end dates and mean pps
plot(0, 0, type = "n", col = "red", bty = "l",
     ylim = c(-.1, 1.1), xlim = c(0, 250),
     xlab = "Time", ylab = "Probability of being in bin 'x'")

day_seq <- 1:dim(out$BUGS$mean$pp)[2]
points(hares$Julian, jitter(hares$White), pch = 19, cex = 1, col = "gray70")

for(i in 1:3){
  lines(day_seq, out$BUGS$mean$pp[i,], col = i, type = "l")
}
abline(v=c(quantile(starts, 0.5)), col="black");abline(v=c(quantile(mids, 0.5)), col="red");abline(v=c(quantile(ends, 0.5)), col="green")
abline(v=c(quantile(starts, 0.025), quantile(starts, 0.975)), col = "black", lty = 3)
abline(v=c(quantile(mids, 0.025), quantile(mids, 0.975)), col = "red", lty = 3)
abline(v=c(quantile(ends, 0.025), quantile(ends, 0.975)), col = "green", lty = 3)
hist(starts, add = T, freq = F, col = "black", border = "black")
hist(ends, add = T, freq = F, col = "green", border = "green")  
hist(mids, add = T, freq = F, col = "red", border = "red")  
text(0, 0.2, paste("CO all, 5K/5H",
                   #"\nyear_eff15 =", quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.025),quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.5),quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.925),
                   #"\nyear_eff16 =", quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.025),quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.5),quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.925),
                   "\nStarts =", quantile(starts, 0.025),quantile(starts, 0.5),quantile(starts, 0.975),
                   "\nMids =", quantile(mids, 0.025),quantile(mids, 0.5),quantile(mids, 0.975),
                   "\nEnds =", quantile(ends, 0.025),quantile(ends, 0.5),quantile(ends, 0.975)), pos = 4, cex=0.9)

# end main spring plot

################################################################################
# FALL PLOT
################################################################################
#  Find start dates
starts <- apply(out$BUGS$sims.list$pp[,1,], 1, function(x){ 
  min(which(x < 0.9)) })
hist(starts, xlab = "Day");quantile(starts, c(0.025, 0.5, 0.975))
#starts.5 <- apply(out$BUGS$sims.list$pp[,1,], 1, function(x){min(which(x < 0.5)) }); hist(starts.5, xlab = "Day");quantile(starts.5, c(0.025, 0.5, 0.975))

#  Find end dates
ends <- apply(out$BUGS$sims.list$pp[,3,], 1, function(x){ 
  min(which(x > 0.9)) })
hist(ends, xlab = "Day");quantile(ends, c(0.025, 0.5, 0.975))
#ends.5 <- apply(out$BUGS$sims.list$pp[,3,], 1, function(x){min(which(x > 0.5)) }); hist(ends.5, xlab = "Day");quantile(ends.5, c(0.025, 0.5, 0.975))

#  Find mid-points
mids <- apply(out$BUGS$sims.list$pp[,2,], 1, function(x){ 
  min(which(x == max(x)))})
hist(mids, xlab = "Day");quantile(mids, c(0.025, 0.5, 0.975))

#Plot start and end dates and mean pps
plot(0, 0, type = "n", col = "red", bty = "l",
     ylim = c(-.1, 1.1), xlim = c(200, 365),
     xlab = "Time", ylab = "Probability of being in bin 'x'")

day_seq <- 1:dim(out$BUGS$mean$pp)[2]
points(hares$Julian, jitter(hares$White), pch = 19, cex = 1, col = "gray60")

for(i in 1:3){
  lines(day_seq, out$BUGS$mean$pp[i,], col = i, type = "l")
}
abline(v=c(quantile(starts, 0.5)), col="black");abline(v=c(quantile(ends, 0.5)), col="green");abline(v=c(quantile(mids, 0.5)), col="red")
abline(v=c(quantile(starts, 0.025), quantile(starts, 0.975)), col = "black", lty = 3)
abline(v=c(quantile(mids, 0.025), quantile(mids, 0.975)), col = "red", lty = 3)
abline(v=c(quantile(ends, 0.025), quantile(ends, 0.975)), col = "green", lty = 3)
hist(starts, add = T, freq = F, col = "black", border = "black");hist(ends, add = T, freq = F, col = "green", border = "green");hist(mids, add = T, freq = F, col = "red", border = "red")
# abline(v=c(quantile(starts.5, 0.5)), col="blue");abline(v=c(quantile(ends.5, 0.5)), col="blue")
# abline(v=c(quantile(starts.5, 0.025), quantile(starts.5, 0.975)), col = "green", lty = 3)
# abline(v=c(quantile(mids.5, 0.025), quantile(mids.5, 0.975)), col = "red", lty = 3)
# abline(v=c(quantile(ends.5, 0.025), quantile(ends.5, 0.975)), col = "black", lty = 3)
text(200, 0.3, paste("CAN all, 5K/500",
                     "\nyear_eff15 =", quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.025),quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.5),quantile(signif(out$BUGS$sims.list$year_eff[,1],digits=2),0.925),
                     "\nyear_eff16 =", quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.025),quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.5),quantile(signif(out$BUGS$sims.list$year_eff[,2],digits=2),0.925),
                     #"\nelev_eff1 =", quantile(signif(out$BUGS$sims.list$elev_eff[,1],digits=2),0.025),quantile(signif(out$BUGS$sims.list$elev_eff[,1],digits=2),0.5),quantile(signif(out$BUGS$sims.list$elev_eff[,1],digits=2),0.925),
                     #"\nelev_eff2 =", quantile(signif(out$BUGS$sims.list$elev_eff[,2],digits=2),0.025),quantile(signif(out$BUGS$sims.list$elev_eff[,2],digits=2),0.5),quantile(signif(out$BUGS$sims.list$elev_eff[,2],digits=2),0.925),
                     "\nStarts=", quantile(starts, 0.025),quantile(starts, 0.5),quantile(starts, 0.975),
                     "\nMids =", quantile(mids, 0.025),quantile(mids, 0.5),quantile(mids, 0.975),
                     "\nEnds=", quantile(ends, 0.025),quantile(ends, 0.5),quantile(ends, 0.975)),
     #"\nStarts.5 =", quantile(starts.5, 0.025),quantile(starts.5, 0.5),quantile(starts.5, 0.975),
     #"\nEnds.5 =", quantile(ends.5, 0.025),quantile(ends.5, 0.5),quantile(ends.5, 0.975)), 
     pos = 4, cex=0.9)

# the only thing copied from fall dev is the plot above
########################################################################################################################


########################################################################################################################
#  Plot with random effects
plot(0, 0, type = "n", col = "red", bty = "l",
     ylim = c(-.1, 1.1), xlim = c(0, 200),
     xlab = "Time",ylab = "Probability of being in bin 'x'")

day_seq <- 1:dim(out$BUGS$mean$pp)[2]

pr_dim <- dim(out$BUGS$mean$p_rand)
ncamera <- pr_dim[2]
ncategories <- pr_dim[1]
nday <- pr_dim[3]

#  Save values of per camera bin probabilities for export
mat <- expand.grid( 
  cats = 1:ncategories,
  cam = 1:ncamera,
  day = day_seq
)

# Create df with all cameras
for(i in 1:nrow(mat)){
  mat$bin_prob[i] <- out$BUGS$mean$p_rand[
    mat$cats[i],
    mat$cam[i],
    mat$day[i]
    ]
}
#mat

#  Add lines to plot for each camera
points(hares$Julian, jitter(hares$White3/100), pch = 19, cex = 1, col = "gray80")
for(i in 1:ncategories){
  for(j in 1:ncamera){
    lines(day_seq, out$BUGS$mean$p_rand[i,j,], col = "gray70", type = "l") #or col =i
  }
}

for(i in 1:3){
  lines(day_seq, out$BUGS$mean$pp[i,], col = i, type = "l")
}
abline(v=c(quantile(starts, 0.5)), col="black");abline(v=c(quantile(mids, 0.5)), col="red");abline(v=c(quantile(ends, 0.5)), col="green")
abline(v=c(quantile(starts, 0.025), quantile(starts, 0.975)), col = "black", lty = 3)
abline(v=c(quantile(mids, 0.025), quantile(mids, 0.975)), col = "red", lty = 3)
abline(v=c(quantile(ends, 0.025), quantile(ends, 0.975)), col = "green", lty = 3)
hist(starts, add = T, freq = F, col = "black", border = "black")
hist(ends, add = T, freq = F, col = "green", border = "green")  
hist(mids, add = T, freq = F, col = "red", border = "red")  
text(0, 0.2, paste("Starts =", quantile(starts, 0.025),quantile(starts, 0.5),quantile(starts, 0.975),
                   "\nMids =", quantile(mids, 0.025),quantile(mids, 0.5),quantile(mids, 0.975),
                   "\nEnds =", quantile(ends, 0.025),quantile(ends, 0.5),quantile(ends, 0.975)), pos = 4, cex=0.9)
# legend("topright",legend = paste(c(0,50,100),"% white"),lty = 1,col = 1:3)



################################################################################ 
#Diagnostics plots
mcmcplot(out, parms = c("pp", "beta", "alpha", "sigma", "rho","elev_eff"#
                        #, "p_rand", "cat_mu")
))
################################################################################ 
hist(out$BUGS$sims.list$elev_eff[,1], breaks = 200);hist(out$BUGS$sims.list$elev_eff[,2], breaks = 200)
hist(out$BUGS$sims.list$sigma[,1], breaks = 200);hist(out$BUGS$sims.list$sigma[,2], breaks = 200)
hist(out$BUGS$sims.list$rho, breaks = 200)
hist(out$BUGS$sims.list$alpha[,1], breaks=200);hist(out$BUGS$sims.list$alpha[,2], breaks=200)
hist(out$BUGS$sims.list$beta[,1], breaks=200);hist(out$BUGS$sims.list$beta[,2], breaks=200)

out.mcmc <- as.mcmc(out) # Convert model output into an MCMC object
str(out.mcmc)
library(coda)
plot(out.mcmc)
#out.mtx <- as.matrix(out.mcmc)
#out.df <- as.data.frame(out.mcmc) # all itterations
#mymodel.p <- out.df[, grep("p[", colnames(out.df), fixed=T)] #only p's
#write.csv(mymodel.p, file = "model.csv") # all itterations for p

#print(out$BUGS$sd) # or instead of mean: sd, median
#a <-print(out$BUGSoutput$sims.array)
#str(out)

hist(out$BUGS$sims.list$beta[,1])
plot(density(out$sims.matrix[,"pp"]))

jpeg(file =/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/data/SSH/myplot.jpeg)
#or
library(lattice)
xyplot(out.mcmc, layout=c(10,10), aspect="fill") # chains history
dev.off()


densityplot(out.mcmc, layout=c(10,10), aspect="fill") # posteriors
#autocorr.plot(out.mcmc) # autocorrelation plot
#gelman.plot(out.mcmc) 

#####