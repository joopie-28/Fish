### Creating a complete model using covariates such as bin lag and timeseries richness ####

# Require packages
library(data.table)

# Create a giant data frame for individual communities

frame_builder_function <- function(fish_community, 
                                   bin_width, 
                                   data_type){
  
  if (bin_width == 2){
    
    if(data_type == "A"){
      
      matrices <- lapply(c(1:5), function(x){ 
        
        lapply(names(fish_community$BioRealm_Novelty_A_2[[x]]), 
               function(y){
                 
                 temp <- fish_community$BioRealm_Novelty_A_2[[x]][[y]]
                 temp$position <- 0
                 
                 for (i in 1:nrow(temp)) {
                   
                   temp$position[i] <- (i + 5)
                   
                 }
                 return(temp)
                 
               })
      })

    }
    
    else {
      
      matrices <- lapply(c(1:5), function(x){ 
        
        lapply(names(fish_community$BioRealm_Novelty_B_2[[x]]), 
               function(y){
                 
                 temp <- fish_community$BioRealm_Novelty_B_2[[x]][[y]]
                 temp$position <- 0
                 
                 for (i in 1:nrow(temp)) {
                   
                   temp$position[i] <- (i + 5)
                   
                 }
                 return(temp)
                 
               })
      })
      
    }
  }
  
  if (bin_width == 1){
    
    if(data_type == "A"){
    
      matrices <- lapply(c(1:5), function(x){ 
        
        lapply(names(fish_community$BioRealm_Novelty_A[[x]]), 
               function(y){
                 
                 temp <- fish_community$BioRealm_Novelty_A[[x]][[y]]
                 temp$position <- 0
                 
                 for (i in 1:nrow(temp)) {
                   
                   temp$position[i] <- (i + 5)
                   
                 }
                 return(temp)
                 
               })
      })
    
    }
    else {
      
      matrices <- lapply(c(1:5), function(x){ 
        
        lapply(names(fish_community$BioRealm_Novelty_B[[x]]), 
               function(y){
                 
                 temp <- fish_community$BioRealm_Novelty_B[[x]][[y]]
                 temp$position <- 0
                 
                 for (i in 1:nrow(temp)) {
                   
                   temp$position[i] <- (i + 5)
                   
                 }
                 return(temp)
                 
               })
      })
      
    }
  }
  
  # Bind everything together into one large frame which we will use 
  # in our models
  
  temp_frame <- rbind(rbindlist(matrices[[1]]), 
                      rbindlist(matrices[[2]]),
                      rbindlist(matrices[[3]]),
                      rbindlist(matrices[[4]]))
  
  
  
  # Add the BioRealm variable in case we want to use it
  
  temp_frame$BioRealm <- 0
  
  # Use a loop to fill it (might take some time)
  
  for (i in 1:nrow(temp_frame)) {
    print(i)
    for (j in 1:nrow(time_series_data)){
      if (temp_frame$site[i] == time_series_data$TimeSeriesID[j]){
        temp_frame$BioRealm[i] <- time_series_data$BioRealm[j]
      }
    }
  }
  
  
  # Ensure the other variables are all in the right format and tidy
  
  temp_frame$BioRealm <- as.factor(temp_frame$BioRealm)
  temp_frame$TimeSeries_Length <- as.numeric(temp_frame$n)
  temp_frame$bins <- as.numeric(temp_frame$bins)
  
  # Convert Boolean to binary values, this ensures a binary response variable
  
  test_frame <- temp_frame [, c("instant", "cumul", "novel")] + 0
  temp_frame[, c("instant", "cumul", "novel")] <- test_frame[,c("instant", "cumul", "novel")]
  
  # We now also want to add timeseries richness 
  
  temp_frame$richness <- NA
  
  # A simple loop gives us richness for the entire timeseries
  for (i in 1:nrow(temp_frame)) {
    print(i)
    
    sub <- subset(Survey_Data, 
                  TimeSeriesID == temp_frame$site[i])
    
    temp_frame$richness[i] <- length(unique(sub$Species))
  }
  
  temp_frame <- temp_frame[, c("instant", "cumul", "novel", "site", "bins", "bin.lag", "position", "BioRealm", "TimeSeries_Length", "richness")]
  
  return(temp_frame)

}

model_input <- frame_builder_function(Fish_Communities_B, 1, "B")

# This is a GLMM with site or timeseriesID treated as a random intercept

covariate_full_glmm <- glmer(novel ~  bin.lag + TimeSeries_Length + (1|site) + position + richness,
                       data=model_input, family=binomial)

summary(covariate_full_glmm)

model_input_B_1 <- model_input

saveRDS(model_input_B_1, "./outputs/GLM_input_B_1.rds")
