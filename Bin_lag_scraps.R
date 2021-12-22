### Performing a logistic regression which checks if bin lag is influencing novelty classificaiton ###

# Require packages
library(data.table)

# Create a giant data frame for binary values
com_frame_B <- rbind(rbindlist(Fish_Communities_B$BioRealm_Novelty_B$palearctic_novelty_B), 
                   rbindlist(Fish_Communities_B$BioRealm_Novelty_B$nearctic_novelty_B),
                   rbindlist(Fish_Communities_B$BioRealm_Novelty_B$afrotropics_novelty_B),
                   rbindlist(Fish_Communities_B$BioRealm_Novelty_B$neotropics_novelty_B),
                   rbindlist(Fish_Communities_B$BioRealm_Novelty_B$australasia_novelty_B))

com_frame_B_2 <- rbind(rbindlist(Fish_Communities_B$BioRealm_Novelty_B_2$palearctic_novelty_B), 
                     rbindlist(Fish_Communities_B$BioRealm_Novelty_B_2$nearctic_novelty_B),
                     rbindlist(Fish_Communities_B$BioRealm_Novelty_B_2$afrotropics_novelty_B),
                     rbindlist(Fish_Communities_B$BioRealm_Novelty_B_2$neotropics_novelty_B))

# Do the same for abundance values
com_frame_A <- rbind(rbindlist(Fish_Communities_A$BioRealm_Novelty_A$palearctic_novelty_A), 
                       rbindlist(Fish_Communities_A$BioRealm_Novelty_A$nearctic_novelty_A),
                       rbindlist(Fish_Communities_A$BioRealm_Novelty_A$afrotropics_novelty_A),
                       rbindlist(Fish_Communities_A$BioRealm_Novelty_A$neotropics_novelty_A),
                       rbindlist(Fish_Communities_A$BioRealm_Novelty_A$australasia_novelty_A))

com_frame_A_2 <- rbind(rbindlist(Fish_Communities_A$BioRealm_Novelty_A_2$palearctic_novelty_A), 
                          rbindlist(Fish_Communities_A$BioRealm_Novelty_A_2$nearctic_novelty_A),
                          rbindlist(Fish_Communities_A$BioRealm_Novelty_A_2$afrotropics_novelty_A),
                          rbindlist(Fish_Communities_A$BioRealm_Novelty_A_2$neotropics_novelty_A))



# Add a biorealm variable
com_frame_A_2$BioRealm <- 0
com_frame_A$BioRealm <- 0
com_frame_B_2$BioRealm <- 0
com_frame_B$BioRealm <- 0


for (i in 1:nrow(com_frame_B)) {
  print(i)
  for (j in 1:nrow(time_series_data)){
    if (com_frame_B$site[i] == time_series_data$TimeSeriesID[j]){
      com_frame_B$BioRealm[i] <- time_series_data$BioRealm[j]
    }
  }
}

for (i in 1:nrow(com_frame_A)) {
  print(i)
  for (j in 1:nrow(time_series_data)){
    if (com_frame_A$site[i] == time_series_data$TimeSeriesID[j]){
      com_frame_A$BioRealm[i] <- time_series_data$BioRealm[j]
    }
  }
}

for (i in 1:nrow(com_frame_B_2)) {
  print(i)
  for (j in 1:nrow(time_series_data)){
    if (com_frame_B_2$site[i] == time_series_data$TimeSeriesID[j]){
      com_frame_B_2$BioRealm[i] <- time_series_data$BioRealm[j]
    }
  }
}

for (i in 1:nrow(com_frame_A_2)) {
  print(i)
  for (j in 1:nrow(time_series_data)){
    if (com_frame_A_2$site[i] == time_series_data$TimeSeriesID[j]){
      com_frame_A_2$BioRealm[i] <- time_series_data$BioRealm[j]
    }
  }
}



com_frame_B$BioRealm <- as.factor(com_frame_B$BioRealm)
com_frame_B$TimeSeries_Length <- com_frame_B$n
com_frame_B$bins <- as.numeric(com_frame_B$bins)

com_frame_A$BioRealm <- as.factor(com_frame_A$BioRealm)
com_frame_A$TimeSeries_Length <- com_frame_A$n
com_frame_A$bins <- as.numeric(com_frame_A$bins)

com_frame_A_2$BioRealm <- as.factor(com_frame_A_2$BioRealm)
com_frame_A_2$TimeSeries_Length <- com_frame_A_2$n
com_frame_A_2$bin_lag_corrected <- com_frame_A_2$bin.lag/2
com_frame_A_2$bins <- as.numeric(com_frame_A_2$bins)

com_frame_B_2$BioRealm <- as.factor(com_frame_B_2$BioRealm)
com_frame_B_2$TimeSeries_Length <- com_frame_B_2$n
com_frame_B_2$bin_lag_corrected <- com_frame_B_2$bin.lag/2
com_frame_B_2$bins <- as.numeric(com_frame_B_2$bins)


# Convert Boolean to binary values
test_frame <- communities_A_2 [, c("instant", "cumul", "novel") ] + 0
communities_A_2[, c("instant", "cumul", "novel")] <- test_frame[,c("instant", "cumul", "novel")]


# Generalized linear model 

bin_lag_model <- glmer(novel ~ bin.lag + TimeSeries_Length + (1|site),
                     data=com_frame, family=binomial)


bin_lag_model_dropped <- glm(novel ~ bin.lag + TimeSeries_Length,
                              data=com_frame, family=binomial)

# Get coefficients
sim.res <- simulateResiduals(bin_lag_model_A)
disp.test <- testDispersion(sim.res)
print(ifelse(disp.test$p.value <=0.05, paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),"Dispersal okay :)"))
  

  
pred.df <- as.data.frame(summary(bin_lag_model_A)$coefficients)
pred.df$taxa.rand <- summary(bin_lag_model_A)$varcor$taxa[1,1]
pred.df$length.rand <- summary(bin_lag_model_A)$varcor$length[1,1]

write.csv(pred.df,
          date.wrap(paste0("./outputs/", novelty_type,
                           " GLMM summary"),".csv")) 



# Get coefficients
summary(bin_lag_model)

summary(bin_lag_model_A)

# Plot some data

spline_plotter <- function(community){
  
  a <- (by(data = community$novel, INDICES = community$bins, FUN = mean))
  w <- (by(data = community$novel, INDICES = community$bins, FUN = length))
  b <- as.numeric(names(a))
  w <- as.vector(w)
  a <- as.vector(a)
  A <- as.data.frame(cbind(b, a, w))
  
  a_c <- (by(data = community$cumul, INDICES = community$bins, FUN = mean))
  w_c <- (by(data = community$cumul, INDICES = community$bins, FUN = length))
  b_c <- as.numeric(names(a_c))
  w_c <- as.vector(w_c)
  a_c <- as.vector(a_c)
  A_c <- as.data.frame(cbind(b_c, a_c, w_c))
  
  a_i <- (by(data = community$instant, INDICES = community$bins, FUN = mean))
  w_i <- (by(data = community$instant, INDICES = community$bins, FUN = length))
  b_i <- as.numeric(names(a_i))
  w_i <- as.vector(w_i)
  a_i <- as.vector(a_i)
  A_i <- as.data.frame(cbind(b_i, a_i, w_i))

# Smoothing spline

spline_plot <- smooth.spline(A$b, A$a, w = w, df = 10) %>%
    broom::augment() %>%
    ggplot(aes(x = b)) + ylim(c(0,0.15)) +
    geom_point(aes(y = a), colour = 'gold') +
    geom_line(aes(y=.fitted), colour = 'gold') + 
  xlab("Bin Age") + ylab('Probability of Novelty (%)') +
  
  geom_point(aes(y = a_c), colour = 'blue') +
  
  geom_point(aes(y = a_i), colour = 'red')


  
  

  show(spline_plot)
  return(spline_plot)
  
}



a <- spline_plotter(communities_B)


test_stuff$length <- rowSums(test_stuff[, c("back", "instant", "cumul", "novel")])

test_stuff$richness <- NA
for (i in 1:nrow(test_stuff)) {
  print(i)
  sub <- subset(Survey_Data, TimeSeriesID == test_stuff$TimeSeries_ID[i])
  test_stuff$richness[i] <- length(unique(sub$Species))
}

test_stuff$average_lag <- NA

means <- by(communities_A_2$bin.lag, INDICES = communities_A_2$site, FUN = mean)

for (i in 1:nrow(test_stuff)) {
  print(i)
  for (j in 1:nrow(means)) {
    if (test_stuff$TimeSeries_ID[i] == names(means[j])){
      test_stuff$average_lag[i] <- means[j] 
    
    }
  }
}

taxa.prob.m <- glm(cbind(novel, non.novel) ~ 1 + log(length) + log(richness) + log(average_lag),
                   data=test_stuff, family=binomial)

summary(taxa.prob.m)

a <- predict(taxa.prob.m, type = "response")
x <- c(1:1331)
par(mfrow = c(1,1))
plot(x, a, main = "DF = 1")
lines(x, exp(mu), lwd = 2, col = "grey")
lines(x, a, col = "orange", lwd = 2)









