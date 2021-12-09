### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the presence/absence data pathway. Remember to load in functions, which are stored in separate files.


# Assign individual time series to BioRealm groups

palearctic_ID <- as.list(subset(time_series_data, BioRealm == "Palearctic")[,3])

nearctic_ID <- as.list(subset(time_series_data, BioRealm == "Nearctic")[,3])

afrotropics_ID <- as.list(subset(time_series_data, BioRealm == "Afrotropics")[,3])

neotropics_ID <- as.list(subset(time_series_data, BioRealm == "Neotropics")[,3])

australasia_ID <- as.list(subset(time_series_data, BioRealm == "Australasia")[,3])

# Store them in a list

check <- list(palearctic_ID, nearctic_ID, afrotropics_ID, neotropics_ID, australasia_ID)
names(check) <- c("palearctic_ID", "nearctic_ID", "afrotropics_ID", "neotropics_ID", "australasia_ID")

# Store the presence/absence matrices in the list according to BioRealm

palearctic_matrices <- lapply(check$palearctic_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- binary_converter_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(palearctic_matrices) <- check$palearctic_ID

nearctic_matrices <- lapply(check$nearctic_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- binary_converter_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(nearctic_matrices) <- check$nearctic_ID

afrotropics_matrices <- lapply(check$afrotropics_ID, function(TimeSeries_ID){
 
  print(TimeSeries_ID)
  temp <- binary_converter_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(afrotropics_matrices) <- check$afrotropics_ID

neotropics_matrices <- lapply(check$neotropics_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- binary_converter_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(neotropics_matrices) <- check$neotropics_ID

australasia_matrices <- lapply(check$australasia_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- binary_converter_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(australasia_matrices) <- check$australasia_ID

# This is just a placeholder, the final data list will have a different configuration
matrix_lists <- list(palearctic_matrices, nearctic_matrices, afrotropics_matrices, neotropics_matrices, australasia_matrices)
names(matrix_lists) <- c("palearctic_mat","nearctic_mat","afrotropics_mat","neotropics_mat", "australasia_mat")


# Compute the novelty and add them to the master list
palearctic_novelty <- lapply(names(matrix_lists$palearctic_mat[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists$palearctic_mat[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  # This line had to be added because there was one timeseries with 0 change over 20 years....
  if (ID == "G7555"){
    return(NA)
  }
  else{
    if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
    else{
      if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
      
        temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                             gam.max.k = -1)
        # Remove first 5 bins
        temp <- temp[-c(1:5),]
        return(temp)
      }
      else {
        return(NA)
      }
    }
  }  
})
names(palearctic_novelty) <- check$palearctic_ID
# Remove NA's
palearctic_novelty <- palearctic_novelty[!sapply(palearctic_novelty, function(x) all(is.na(x)))]

nearctic_novelty <- lapply(names(matrix_lists$nearctic_mat[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists$nearctic_mat[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
    
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                               gam.max.k = -1)
      # Remove first 5 bins
      temp <- temp[-c(1:5),]
      return(temp)
    }
    else {
      return(NA)
    }
  }  
})
names(nearctic_novelty) <- check$nearctic_ID
# Remove NA's
nearctic_novelty <- nearctic_novelty[!sapply(nearctic_novelty, function(x) all(is.na(x)))]

afrotropics_novelty <- lapply(names(matrix_lists$afrotropics_mat[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists$afrotropics_mat[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
    
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                             gam.max.k = -1)
      # Remove first 5 bins
      temp <- temp[-c(1:5),]
      return(temp)
    }
    else {
      return(NA)
    }
  }
})
names(afrotropics_novelty) <- check$afrotropics_ID
# Remove NA's
afrotropics_novelty <- afrotropics_novelty[!sapply(afrotropics_novelty, function(x) all(is.na(x)))]

neotropics_novelty <- lapply(names(matrix_lists$neotropics_mat[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists$neotropics_mat[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
      
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                             gam.max.k = -1)
      # Remove first 5 bins
      temp <- temp[-c(1:5),]
      return(temp)
    }
    else{
      return(NA)
    }
  }
})
names(neotropics_novelty) <- check$neotropics_ID
# Remove NA's
neotropics_novelty <- neotropics_novelty[!sapply(neotropics_novelty, function(x) all(is.na(x)))]

australasia_novelty <- lapply(names(matrix_lists$australasia_mat[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists$australasia_mat[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >= 5) {
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                             gam.max.k = -1)
      # Remove first 5 bins
      temp <- temp[-c(1:5),]
      return(temp)
    }
    else {
      return(NA)
    }
  }
})
names(australasia_novelty) <- check$australasia_ID
# Remove NA's
australasia_novelty <- australasia_novelty[!sapply(australasia_novelty, function(x) all(is.na(x)))]




# Creating the master list, B is for binary (presence/absence)

palearctic_data <- list(palearctic_ID, palearctic_matrices, palearctic_novelty)
names(palearctic_data) <- c("palearctic_ID", "palearctic_matrices_B", 
                            "palearctic_novelty_B")

nearctic_data <- list(nearctic_ID, nearctic_matrices, nearctic_novelty)
names(nearctic_data) <- c("nearctic_ID", "nearctic_matrices_B", 
                            "nearctic_novelty_B")

afrotropics_data <- list(afrotropics_ID, afrotropics_matrices, afrotropics_novelty)
names(afrotropics_data) <- c("afrotropics_ID", "afrotropics_matrices_B", 
                          "afrotropics_novelty_B")

neotropics_data <- list(neotropics_ID, neotropics_matrices, neotropics_novelty)
names(neotropics_data) <- c("neotropics_ID", "neotropics_matrices_B", 
                             "neotropics_novelty_B")

australasia_data <- list(australasia_ID, australasia_matrices, australasia_novelty)
names(australasia_data) <- c("australasia_ID", "australasia_matrices_B", 
                            "australasia_novelty_B")

# Here it is
Fish_Communities <- list(palearctic_data, nearctic_data, afrotropics_data, neotropics_data, australasia_data)
names(Fish_Communities) <- c("palearctic_data", "nearctic_data", "afrotropics_data", "neotropics_data", "australasia_data")



