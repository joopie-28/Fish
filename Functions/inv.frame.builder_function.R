# creates the return data needed for modelling

inv.frame.builder <- function(converted_matrix_list){
  
  list_mat <- lapply(names(converted_matrix_list[]), function(ID){

    print(ID)
    
    matrix_1 <- converted_matrix_list[[ID]]
    
    if(any(is.na(matrix_1))){
      return(NA)
    }
    
    c.ID <- strsplit(ID, split = "-")[[1]][2]
    n.ID <- strsplit(ID, split = "-")[[1]][1]
    print(n.ID)
    
    # Lets label!
    
    if(c.ID == "USA" | c.ID == "CAN"){
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[n.ID]]
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[n.ID]], 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = F, 
                                        site = n.ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    
    if(c.ID == "BRA" | c.ID == "COL"){
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[n.ID]]
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[n.ID]], 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = T, 
                                        site = n.ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    
    if(c.ID == "BWA" | c.ID == "CIV"){
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[n.ID]]
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[n.ID]], 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = F, 
                                        site = n.ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    } 
    
    if(c.ID == "FRA" | c.ID == "GBR"| c.ID == "ESP" | c.ID == "FIN"| c.ID == "BEL" | c.ID == "SWE"|c.ID == "HUN" | c.ID == "JPN") {
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[n.ID]]
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[n.ID]], 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = T, 
                                        site = n.ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    
    if(c.ID == "AUS"){
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[n.ID]]
      label_frame <- identify.novel.gam(Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[n.ID]], 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = T, 
                                        site = n.ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    
    
    matrix_1$cat <- label_frame$cat
    matrix_1$site <- n.ID
    matrix_1$position <- 0
    matrix_1$bin_lag <- label_frame$bin.lag
    matrix_1$country <- 0
    matrix_1$Latitude <- 0
    matrix_1$Longitude <- 0
    matrix_1$novel <- label_frame$novel + 0
    matrix_1$instant <- label_frame$instant + 0
    matrix_1$cumul <- label_frame$cumul + 0
    matrix_1$total.n <- label_frame$n
    

    # evenness
    
    even.vector <- NULL
    shannon.vector <- NULL
    diversity.vector <- NULL
    diversity.prev.vector <- NA
    diversity.next.vector <- NULL
    for (i in 1:nrow(matrix)){
      # Get number of species
      n.sp <- ncol(matrix[,matrix[i,] > 0, drop = F])
      add <- diversity(matrix[i,], index = "shannon")/log(n.sp)
      add.s <- diversity(matrix[i,], index = "shannon")
      even.vector <- c(even.vector, add)
      shannon.vector <- c(shannon.vector, add.s)
      diversity.vector<- c(diversity.vector, n.sp)
      diversity.prev.vector <- c(diversity.prev.vector, n.sp)
      diversity.next.vector <- c(diversity.next.vector, n.sp)
    }
    
    
    diversity.next.vector <- c(diversity.next.vector, NA)
    
    matrix_1$evenness <- even.vector
    matrix_1$shannon.d <- shannon.vector
    matrix_1$diversity <- diversity.vector
    matrix_1$diversity.previous <- diversity.prev.vector[-length(diversity.prev.vector)]
    matrix_1$diversity.next <- diversity.next.vector[-1]

    
    # HydroBasin
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$basin <- time_series_data$HydroBasin[i]
      }
    }
    
    # Lat and Long
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$Latitude <- time_series_data$Latitude[i]
        matrix_1$Longitude <- time_series_data$Longitude[i]
      }
    }
    
    
    # Country 
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$country <- time_series_data$Country[i]
      }
    }
    
    for(i in 1:nrow(matrix_1)){
      matrix_1$position[i] <- i
    }
    
    ### Local origination
    
    
    sub_series <- subset(Survey_Data, TimeSeriesID == n.ID)
    
    # Create a vector with years
    years <- sort(unique(sub_series$Year))
    
    # Apply over all years to find the species list for each year
    species_checklist <- lapply(years, FUN = function(x){
      
      checklist <- subset(sub_series, sub_series$Year == x)
      species_checklist <- matrix(unique(checklist[, c("Species")]))
      return(species_checklist)
      
    })
    
    loc.orig.vector <- NA
    checker <- species_checklist
    
    for (i in 2:length(species_checklist)){
      
      new.sp <- length(checker[[i]][checker[[i]] %!in% checker[[i-1]]])
      
      loc.orig.vector <- c(loc.orig.vector, ifelse(length(new.sp) == 0, 0, new.sp))
      
      checker[[i]] <- unique(c(checker[[i]], checker[[i-1]]))
      
    }
    
    matrix_1$orig <- loc.orig.vector
    
    ### Local extinction
    
    loc.ext.vector <- NA
    checker.e <- species_checklist
    
    for (i in (length(species_checklist)-1):1){
      
      ext.sp <- length(checker.e[[i]][checker.e[[i]] %!in% checker.e[[i+1]]])
      
      loc.ext.vector <- c(loc.ext.vector, ifelse(length(ext.sp) == 0, 0, ext.sp))
      
      checker.e[[i]] <- unique(c(checker.e[[i]], checker.e[[i+1]]))
      
    }
    
    matrix_1$ext <- rev(loc.ext.vector)
    
    
    
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
                             "bin_lag", "basin", "Longitude", "Latitude", "NNC", 
                             "NNC_increase", "NAC", "NAC_increase", "INC", "INC_increase", 
                             "bins", "novel", "instant", "cumul", "country", "evenness", "orig", "ext", "shannon.d", 'diversity', 'diversity.previous', 'diversity.next',
                             'total.n')]
  
    return(matrix_1)
    
  
    
    
    
  })
  
  list_mat_2 <- list_mat[!is.na(list_mat)]
  
  list_mat_2 <- rbindlist(list_mat_2)
  
  list_mat_2$BioRealm <- NA
  country_biorealm <- unique(data.frame("country" = time_series_data$Country, "BioRealm" = time_series_data$BioRealm))
  for (i in 1:nrow(list_mat_2)){
    print("Biorealms")
    for (j in 1:nrow(country_biorealm)){
      
      if(list_mat_2$country[i] == country_biorealm$country[j]){
        
        list_mat_2$BioRealm[i] <- country_biorealm$BioRealm[j]
        next
      }
      
    }
    
  }
  
  # center all features
  
  #list_mat_2$NNC <- scale(list_mat_2$NNC, center = T, scale = T)
  #list_mat_2$NNC_increase <- scale(list_mat_2$NNC_increase, center = T, scale = T)
  #list_mat_2$NAC <- scale(list_mat_2$NAC, center = T, scale = T)
  #list_mat_2$NAC_increase <- scale(list_mat_2$NAC_increase, center = T, scale = T)
  #list_mat_2$INC <- scale(list_mat_2$INC, center = T, scale = T)
  #list_mat_2$INC_increase <- scale(list_mat_2$INC_increase, center = T, scale = T)
  #list_mat_2$bin_lag <- scale(as.numeric(list_mat_2$bin_lag, center = T, scale = T ))
  #list_mat_2$position <- scale(list_mat_2$position, center = T, scale = T)
  #list_mat_2$evenness <- scale(list_mat_2$evenness, center = T, scale = T)
  #list_mat_2$orig <- scale(list_mat_2$orig, center = T, scale = T)
  #list_mat_2$ext <- scale(list_mat_2$ext, center = T, scale = T)
  #list_mat_2$shannon.d <- scale(list_mat_2$shannon.d, center = T, scale = T)
  list_mat_2$basin <- as.factor(list_mat_2$basin)
  list_mat_2$site <- as.factor(list_mat_2$site)
  list_mat_2[is.na(list_mat_2)] <- 0
  
  return(list_mat_2)
}

# extra little trick
`%!in%` = Negate(`%in%`)





# Incorporates seasonality
inv.frame.builder.V2 <- function(converted_matrix_list){
  
  list_mat <- lapply(names(converted_matrix_list), function(ID){
  
    matrix_1 <- converted_matrix_list[[ID]]
  
    n.ID <- strsplit(as.character(ID), ".",
                     fixed = TRUE)[[1]][1]
    quarter.ID <-strsplit(as.character(ID), ".",
                          fixed = TRUE)[[1]][2]
    
    print(ID)
    
    # Lets label!
    matrix <- matrix_list_seasonality[[ID]]
    label_frame <- identify.novel.gam(matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    
    
    
    matrix_1$cat <- label_frame$cat
    matrix_1$site <- ID
    matrix_1$position <- 0
    matrix_1$bin_lag <- label_frame$bin.lag
    matrix_1$country <- 0
    matrix_1$Latitude <- 0
    matrix_1$Longitude <- 0
    matrix_1$novel <- label_frame$novel + 0
    matrix_1$instant <- label_frame$instant + 0
    matrix_1$cumul <- label_frame$cumul + 0
    matrix_1$total.n <- label_frame$n
    matrix_1$BioRealm <- NA
    
    
    # evenness
    
    even.vector <- NULL
    shannon.vector <- NULL
    diversity.vector <- NULL
    diversity.prev.vector <- NA
    diversity.next.vector <- NULL
    for (i in 1:nrow(matrix)){
      # Get number of species
      n.sp <- ncol(matrix[,matrix[i,] > 0, drop = F])
      add <- diversity(matrix[i,], index = "shannon")/log(n.sp)
      add.s <- diversity(matrix[i,], index = "shannon")
      even.vector <- c(even.vector, add)
      shannon.vector <- c(shannon.vector, add.s)
    }
    
    
    diversity.next.vector <- c(diversity.next.vector, NA)
    
    matrix_1$evenness <- even.vector
    matrix_1$shannon.d <- shannon.vector
 
   
    # Country 
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$country <- time_series_data$Country[i]
        matrix_1$BioRealm <- time_series_data$BioRealm[i]
      }
    }
    
    for(i in 1:nrow(matrix_1)){
      matrix_1$position[i] <- i
    }
    
    # Lat and Long
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$Latitude <- time_series_data$Latitude[i]
        matrix_1$Longitude <- time_series_data$Longitude[i]
      }
    }
    
  
    
    # HydroBasin
    for (i in 1:nrow(time_series_data)){
      if (n.ID == time_series_data$TimeSeriesID[i]){
        matrix_1$basin <- time_series_data$HydroBasin[i]
      }
    }
    
    
  
    # Remove first 5 as these are biased
    matrix_1 <- matrix_1[-(1:5),] 
    
    # Remove species
    matrix_1 <- matrix_1[, c("cat", "site", "position", 
                             "bin_lag","BioRealm","basin","Longitude", "Latitude", "NNC", 
                             "NNC_increase","NNC_spec","NNC_spec_increase", "NAC", "NAC_increase","NAC_spec","NAC_spec_increase", "INC", "INC_increase", "INC_spec", "INC_spec_increase", 
                             "bins", "novel", "instant", "cumul", "country", "evenness", "shannon.d",
                             'total.n')]
    
    return(matrix_1)
    
  })
  
  list_mat_2 <- list_mat[!is.na(list_mat)]
  
  list_mat_2 <- rbindlist(list_mat_2)
  
 
  
  # center all features
  
  #list_mat_2$NNC <- scale(list_mat_2$NNC, center = T, scale = T)
  #list_mat_2$NNC_increase <- scale(list_mat_2$NNC_increase, center = T, scale = T)
  #list_mat_2$NAC <- scale(list_mat_2$NAC, center = T, scale = T)
  #list_mat_2$NAC_increase <- scale(list_mat_2$NAC_increase, center = T, scale = T)
  #list_mat_2$INC <- scale(list_mat_2$INC, center = T, scale = T)
  #list_mat_2$INC_increase <- scale(list_mat_2$INC_increase, center = T, scale = T)
  #list_mat_2$bin_lag <- scale(as.numeric(list_mat_2$bin_lag, center = T, scale = T ))
  #list_mat_2$position <- scale(list_mat_2$position, center = T, scale = T)
  #list_mat_2$evenness <- scale(list_mat_2$evenness, center = T, scale = T)
  #list_mat_2$orig <- scale(list_mat_2$orig, center = T, scale = T)
  #list_mat_2$ext <- scale(list_mat_2$ext, center = T, scale = T)
  #list_mat_2$shannon.d <- scale(list_mat_2$shannon.d, center = T, scale = T)

  list_mat_2$site <- as.factor(list_mat_2$site)
  list_mat_2[is.na(list_mat_2)] <- 0
  
  return(list_mat_2)
}
