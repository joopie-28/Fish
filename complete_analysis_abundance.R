### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the unaltered abundance data pathway. Remember to load in functions, which are stored in separate files.

# Store the presence/absence matrices in the list according to BioRealm

palearctic_matrices_A <- lapply(check$palearctic_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- abundance_identifier_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(palearctic_matrices_A) <- check$palearctic_ID

nearctic_matrices_A <- lapply(check$nearctic_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- abundance_identifier_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(nearctic_matrices_A) <- check$nearctic_ID

afrotropics_matrices_A <- lapply(check$afrotropics_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- abundance_identifier_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(afrotropics_matrices_A) <- check$afrotropics_ID

neotropics_matrices_A <- lapply(check$neotropics_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- abundance_identifier_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(neotropics_matrices_A) <- check$neotropics_ID

australasia_matrices_A <- lapply(check$australasia_ID, function(TimeSeries_ID){
  
  print(TimeSeries_ID)
  temp <- abundance_identifier_function(TimeSeries_ID)
  if(class(temp)=="matrix"){
    # remove bins within sites with fewer species than cut off
    temp <- temp[rowSums(temp) > rich.cutoff,]
  }
  return(temp)
})
names(australasia_matrices_A) <- check$australasia_ID

# This is just a placeholder, the final data list will have a different configuration
matrix_lists_A <- list(palearctic_matrices_A, nearctic_matrices_A, afrotropics_matrices_A, neotropics_matrices_A, australasia_matrices_A)
names(matrix_lists_A) <- c("palearctic_mat_A","nearctic_mat_A","afrotropics_mat_A","neotropics_mat_A", "australasia_mat_A")


# Compute the novelty and add them to the master list
palearctic_novelty_A <- lapply(names(matrix_lists_A$palearctic_mat_A[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists_A$palearctic_mat_A[ID])[]
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
        
        temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
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
names(palearctic_novelty_A) <- check$palearctic_ID
# Remove NA's
palearctic_novelty_A <- palearctic_novelty_A[!sapply(palearctic_novelty_A, function(x) all(is.na(x)))]

nearctic_novelty_A <- lapply(names(matrix_lists_A$nearctic_mat_A[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists_A$nearctic_mat_A[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
      
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
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
names(nearctic_novelty_A) <- check$nearctic_ID
# Remove NA's
nearctic_novelty_A <- nearctic_novelty_A[!sapply(nearctic_novelty_A, function(x) all(is.na(x)))]

afrotropics_novelty_A <- lapply(names(matrix_lists_A$afrotropics_mat_A[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists_A$afrotropics_mat_A[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
      
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
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
names(afrotropics_novelty_A) <- check$afrotropics_ID
# Remove NA's
afrotropics_novelty_A <- afrotropics_novelty_A[!sapply(afrotropics_novelty_A, function(x) all(is.na(x)))]

neotropics_novelty_A <- lapply(names(matrix_lists_A$neotropics_mat_A[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists_A$neotropics_mat_A[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
      
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
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
names(neotropics_novelty_A) <- check$neotropics_ID
# Remove NA's
neotropics_novelty_A <- neotropics_novelty_A[!sapply(neotropics_novelty_A, function(x) all(is.na(x)))]

australasia_novelty_A <- lapply(names(matrix_lists_A$australasia_mat_A[]), function(ID){
  print(ID)
  site.sp.mat <- (matrix_lists_A$australasia_mat_A[ID])[]
  site.sp.mat <- site.sp.mat[[ID]]
  if(typeof(site.sp.mat) == "character"){
    return(NA)
  }
  else{
    if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >= 5) {
      temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
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
names(australasia_novelty_A) <- check$australasia_ID
# Remove NA's
australasia_novelty_A <- australasia_novelty_A[!sapply(australasia_novelty_A, function(x) all(is.na(x)))]




# Creating the master list, B is for binary (presence/absence)

palearctic_data_A <- list(palearctic_ID, palearctic_matrices_A, palearctic_novelty_A)
names(palearctic_data_A) <- c("palearctic_ID", "palearctic_matrices_A", 
                            "palearctic_novelty_A")

nearctic_data_A <- list(nearctic_ID, nearctic_matrices_A, nearctic_novelty_A)
names(nearctic_data_A) <- c("nearctic_ID", "nearctic_matrices_A", 
                          "nearctic_novelty_A")

afrotropics_data_A <- list(afrotropics_ID, afrotropics_matrices_A, afrotropics_novelty_A)
names(afrotropics_data_A) <- c("afrotropics_ID", "afrotropics_matrices_A", 
                             "afrotropics_novelty_A")

neotropics_data_A <- list(neotropics_ID, neotropics_matrices_A, neotropics_novelty_A)
names(neotropics_data_A) <- c("neotropics_ID", "neotropics_matrices_A", 
                            "neotropics_novelty_A")

australasia_data_A <- list(australasia_ID, australasia_matrices_A, australasia_novelty_A)
names(australasia_data_A) <- c("australasia_ID", "australasia_matrices_A", 
                             "australasia_novelty_A")

# Here it is
Fish_Communities_A <- list(palearctic_data_A, nearctic_data_A, afrotropics_data_A, neotropics_data_A, australasia_data_A)
names(Fish_Communities_A) <- c("palearctic_data_A", "nearctic_data_A", "afrotropics_data_A", "neotropics_data_A", "australasia_data_A")



