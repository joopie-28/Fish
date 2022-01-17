#### This function takes a list of novelty matrices as input and returns a dataframe with 
#### all communities and their respective attributes, which can be used in the GLM's.

# Require packages
library(data.table)

# Create a giant data frame for individual communities
# This allows us to use a binary response variable for the GLM's 
# we use later.

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
  
 
  temp_frame$length <- 0
  temp_frame$years_before <- 0
  temp_frame$richness <- NA
  temp_frame$bins <- as.numeric(temp_frame$bins)
  temp_frame$TimeSeries_Length <- as.numeric(temp_frame$n)
  
  # Use a loop to fill 

  for (i in 1:nrow(temp_frame)){
    print(paste("length, years before and richness", i))
    
    sub <- subset(Survey_Data, 
                  TimeSeriesID == temp_frame$site[i])
    
    temp_frame$length[i] <- (max(unique(sub$Year)-min((unique(sub$Year)))))
    temp_frame$years_before[i] <- ((2021-temp_frame$bins[i]) - min((unique(sub$Year))))
    
    # Richness per community
    if (bin_width == 2) {
      
      sub <- subset(sub,
                    Year == (2021-temp_frame$bins[i]) | Year == (2022-temp_frame$bins[i]))
    }
    
    else {
      sub <- subset(sub,
                    Year == (2021-temp_frame$bins[i])) 
    }
    
    temp_frame$richness[i] <- length(unique(sub$Species))
    
  }
  
    
  
  
  
  # Ensure the other variables are all in the right format and tidy

  
  # Convert Boolean to binary values, this ensures a binary response variable
  
  test_frame <- temp_frame [, c("instant", "cumul", "novel")] + 0
  temp_frame[, c("instant", "cumul", "novel")] <- test_frame[,c("instant", "cumul", "novel")]

  temp_frame <- temp_frame[, c("instant", "cumul", "novel", "site", "bins", "bin.lag", "position", "TimeSeries_Length", "richness", "years_before", "length")]
  
  return(temp_frame)
  
}




