# creates the return data needed for modelling

inv.frame.builder <- function(converted_matrix_list, country){
  
   
  
  list_mat <- lapply(names(converted_matrix_list[]), function(ID){
    print(ID)  
    
    matrix_1 <- converted_matrix_list[[ID]]
    
    if(is.na(matrix_1)){
      return(NA)
    }
    
    
    
    
    # Lets label!
    if(country == "USA"){
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A_2$nearctic_mat_A[[ID]], 
                                        alpha = 0.05,
                                        metric = "jaccard",
                                        plot = TRUE, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    
    else{
    label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A[[ID]], 
                                      alpha = 0.05,
                                      metric = "jaccard",
                                      plot = TRUE, 
                                      site = ID,
                                      plot.data = FALSE,
                                      gam.max.k = -1)
    }
    
    
    matrix_1$cat <- label_frame$cat
    matrix_1$site <- ID
    matrix_1$position <- 0
    matrix_1$bin_lag <- label_frame$bin.lag
    matrix_1$country <- 0
    matrix_1$novel <- label_frame$novel + 0
    matrix_1$instant <- label_frame$instant + 0
    matrix_1$cumul <- label_frame$cumul + 0
    
    # HydroBasin
    for (i in 1:nrow(time_series_data)){
      if (ID == time_series_data$TimeSeriesID[i]){
        matrix_1$basin <- time_series_data$HydroBasin[i]
      }
    }
    
    # Country 
    for (i in 1:nrow(time_series_data)){
      if (ID == time_series_data$TimeSeriesID[i]){
        matrix_1$country <- time_series_data$Country[i]
      }
    }
    
    for(i in 1:nrow(matrix_1)){
      matrix_1$position[i] <- i
    }
    
    # Now that we have added all variables, it is a good idea to scale them for modelling
    #matrix_1$NNC <- scale(matrix_1$NNC, center = T, scale = T)
    #matrix_1$NNC_increase <- scale(matrix_1$NNC_increase, center = T, scale = T)
    #matrix_1$NAC <- scale(matrix_1$NAC, center = T, scale = T)
    #matrix_1$NAC_increase <- scale(matrix_1$NAC_increase, center = T, scale = T)
    #matrix_1$INC <- scale(matrix_1$INC, center = T, scale = T)
    #matrix_1$INC_increase <- scale(matrix_1$INC_increase, center = T, scale = T)
    #matrix_1$bin_lag <- scale(as.numeric(matrix_1$bin_lag, center = T, scale = T ))
    #matrix_1$position <- scale(matrix_1$position, center = T, scale = T)
    #matrix_1$basin <- as.factor(matrix_1$basin)
    #matrix_1[is.na(matrix_1)] <- 0
    
    
    # Remove first 5 as these are biased
    matrix_1 <- matrix_1[-(1:5),] 
    
    # Remove species
    matrix_1 <- matrix_1[, c("cat", "site", "position", 
                             "bin_lag", "basin", "NNC", 
                             "NNC_increase", "NAC", "NAC_increase", "INC", "INC_increase", 
                             "bins", "novel", "instant", "cumul", "country")]
  
    return(matrix_1)
    
  
    
    
    
  })
  
  list_mat_2 <- list_mat[!is.na(list_mat)]
  
  list_mat_2 <- rbindlist(list_mat_2)
  
  list_mat_2$NNC <- scale(list_mat_2$NNC, center = T, scale = T)
  list_mat_2$NNC_increase <- scale(list_mat_2$NNC_increase, center = T, scale = T)
  list_mat_2$NAC <- scale(list_mat_2$NAC, center = T, scale = T)
  list_mat_2$NAC_increase <- scale(list_mat_2$NAC_increase, center = T, scale = T)
  list_mat_2$INC <- scale(list_mat_2$INC, center = T, scale = T)
  list_mat_2$INC_increase <- scale(list_mat_2$INC_increase, center = T, scale = T)
  list_mat_2$bin_lag <- scale(as.numeric(list_mat_2$bin_lag, center = T, scale = T ))
  list_mat_2$position <- scale(list_mat_2$position, center = T, scale = T)
  list_mat_2$basin <- as.factor(list_mat_2$basin)
  list_mat_2$site <- as.factor(list_mat_2$site)
  list_mat_2[is.na(list_mat_2)] <- 0
  
  return(list_mat_2)
}
