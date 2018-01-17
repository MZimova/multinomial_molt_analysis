#  Multinomial analysis to calculate molt start and end dates 
     # (and covariates effects) in SHH
# 3 categories of molt
# Random effects on camera traps
#  06/2017
#  Josh Nowak and Marketa Zimova
################################################################################
library(R2jags)
library(readr)
library(purrr)
library(dplyr)
library(beepr)
library(mcmcplots)
library(ggplot2)
############################################################################
# Load data  

# 1. Load hare data  
load("/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/hares_cl.RData")
hares <- select(hares_cl, Site, White3, White, Julian, Cluster, Year, fCluster, Camera) %>%
  #subset(Julian >50 & Julian <200)# & Site=="NH")# & Year==2015)
  subset(Julian >220 & Site=="CO")# & Year==2015)
# if running with temp remove CO_68: it is missing daymet! 
#hares <-filter(hares,Cluster != 68)
(fClustercount<-length(unique(hares$fCluster)))

# 2. Load camera data (for selected hares data above)
load("/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/cameras.RData")
#camera_names <- unique(hares$Camera)
#cameras <- filter(cameras, Camera %in% camera_names)
# for CO or CLUSTERS:
camera_names <- unique(hares$fCluster)
cameras <- filter(cameras, fCluster %in% camera_names)
(fClustercount <-length(unique(cameras$fCluster)))

# summary <- hares %>%
#   group_by(Cluster,White3) %>%
#   summarize(counts=n())
# ggplot(summary, aes(as.factor(Cluster),White3,size=counts,colour=as.factor(White3))) + 
#   geom_jitter(width = 0, height = .3)

ggplot(hares, aes(Julian, White,colour=as.factor(Cluster))) + geom_jitter(width = 1, height = 1)
################################################################################
#  Call a single model step by step - mimics jags_call
#  Set time_scale for the analysis (pptions are in the column names of hares, Month, Week, Julian)
time_scale <- "Julian" #MUTE and write "Julian" into days directly?
load.module("glm")

#  Subset to days - to reduce redundancy and ease inits and data create
# unlist makes anything into a vector here subsetting for one column: [all rows,column 'Julian']
days <- as.integer(unlist(hares[,time_scale]))
(first_day <- min(days))
(last_day <- max(days))
(my_year <-  as.integer(unique(hares$Year)))

#  Create categorical response
hares$response <- cut(hares$White3, 3, labels = 1:3) # 1 is brown
# response <- cut(hares$White3, 3, labels = 1:3) # used to have (or: response <- ordered(response) to make it ordered but same result)

#  Inits
inits <- function(){
  list(
    alpha = rnorm(3)
  )
}
# 
#  Prepare data - all covariates on camera basis so run following to make camera cluster basis:
GeoCL <- cameras %>%
  #filter(Site =="NH")  %>% # filter for each site if not working with all at once
  group_by(fCluster) %>% # hide if interested in running models on camera basis
  summarize(ElevCl= mean(Elevation),
            LatCl=mean(Lat),
            LonCl=mean(Lon)) %>%
  mutate(ElevSC=scale(ElevCl),
         LatSC=scale(LatCl),
        LonSC=scale(LonCl)) %>%  #centers and scales (both default true)
  arrange(fCluster)

# #  Snow livneh
# load('/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/past_livneh_snowseas.RData')
# livneh <- past_livneh %>%
#   filter(fCluster %in% camera_names) %>%
#   mutate(FirstSC =scale(First),
#        LastSC =scale(Last)) %>%  #centers and scales (both default true)
#   arrange(Cluster)
# length(unique(livneh$fCluster))
# 
# #  Temperature (seasonal 30-yr normal 1980-2009)
# load('/Users/marketzimova/Documents/WORK/DISSERTATION/3 Camera Traps Study/data/past_daymet_seas.RData')
# past_daymet <- past_daymet_seas %>%
#   filter(fCluster %in% camera_names)
# 
# past_daymet<- past_daymet %>%
#   filter(season==1)  %>% # 1 is spring
#   group_by(fCluster) %>% # hide if interested in running models on camera basis
#   summarize(tavg30 =mean(season_30tavg),
#             tmin30 =mean(season_30tmin),
#             tmax30 =mean(season_30tmax),
#             Cluster =mean(Cluster)) %>%
#   mutate(tavg30SC =scale(tavg30),
#          tmin30SC =scale(tmin30),
#          tmax30SC =scale(tmax30)) %>%  #centers and scales (both default true)
#   arrange(Cluster)
# length(unique(past_daymet$fCluster))

# Gather data 
dat <- list(
  nobs = nrow(hares),
  day = days, 
  cam = as.numeric(as.factor(hares$fCluster)), 
  y = hares$response,
  nbins = 3,
  #first_day = first_day, #MZ added
  ndays = last_day,
  ncam = length(unique(hares$fCluster)),
  #tavg30 = as.numeric(GeoCL$ElevSC) #or elev instead of tavg30
  elev = as.numeric(GeoCL$ElevSC)
  #tavg30 = as.numeric(past_daymet$tmax30SC)
  #tavg30 = as.numeric(livneh$LastSC)
)

# Parameters to monitor
parms <- c("beta", "alpha","pp","sigma_cam","tau_cam",
           "elev_eff","tavg30_eff"#, "lat_eff"#,
           #"sigma", "rho", p_rand","cat_mu" 
)

#  Call jags
setwd("/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis") #for the model location
start.time <- Sys.time()
out <- jags.parallel(
  data = dat, 
  inits = NULL,
  parameters.to.save = parms,
  model.file = "models/multinom_covs.txt", #multinom_covs or multinom
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 2000,
  n.thin = 3
)
end.time <- Sys.time();(time.taken <-end.time-start.time)
beep()

#### 
out.sum <- out$BUGS$summary 
#write.table(out.sum, file="/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis/results/test_1thrubin.csv",sep=",")

options(max.print=10000) #extend maximum for print
print(out)

#out$BUGS$mean$lat_eff
recompile(out)
out <- update(out, n.iter=5000)

# cool plot, doesn't work now
# ggplot(hares) +
#   #facet_grid(Year ~ .) +
#   geom_jitter(aes(Julian, response, fill=as.numeric(elev)), width = 1, height = .3,
#               pch=21, colour="Black", size=2.5, alpha=0.65) +
#   ggtitle("New England") +
#   scale_fill_gradientn(colours = rev(rainbow(2))) # or rev(rainbow(2)) or terrain.colors(10)
################################################################################
# Save results out
#save(out, file = "/Users/marketzimova/Documents/WORK/DISSERTATION/GitHub/multinomial_molt_analysis/out_CO_lat_elev.RData")
#load("out_CO_yr_all.RData")
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
points(hares$Julian, jitter(hares$White/100), pch = 19, cex = 1, col = "gray70")

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
text(0, 0.2, paste("Can all w elev, 5K/1H",
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