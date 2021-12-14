### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the presence/absence data pathway. Remember to load in functions, which are stored in separate files.

# Assign individual time series to BioRealm groups

palearctic_ID <- as.list(subset(time_series_data, BioRealm == "Palearctic")[,3])

nearctic_ID <- as.list(subset(time_series_data, BioRealm == "Nearctic")[,3])

afrotropics_ID <- as.list(subset(time_series_data, BioRealm == "Afrotropics")[,3])

neotropics_ID <- as.list(subset(time_series_data, BioRealm == "Neotropics")[,3])

australasia_ID <- as.list(subset(time_series_data, BioRealm == "Australasia")[,3])

# Create a list to hold these, serves as input for the matrix creator
ID_list <- list(palearctic_ID, nearctic_ID, afrotropics_ID, neotropics_ID, australasia_ID )

names(ID_list) <- c("palearctic_ID", "nearctic_ID", "afrotropics_ID", "neotropics_ID", "australasia_ID")
rm(palearctic_ID, nearctic_ID, nearctic_ID, afrotropics_ID, neotropics_ID, australasia_ID)

# Create matrix list using this function
list_matrix_B_function <- function(check_list){

  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
             function(TimeSeries_ID){
               print(TimeSeries_ID)
               temp <- binary_converter_function(TimeSeries_ID)
               if(class(temp)=="matrix"){
                 # remove bins within sites with fewer species than cut off
                 temp <- temp[rowSums(temp) > rich.cutoff,]
                 }
               return(temp)
      }) 
  
    if (names(check_list[i]) == "palearctic_ID"){
      palearctic_mat_B <- nam
      names(palearctic_mat_B) <- check_list$palearctic_ID
    }
    if (names(check_list[i]) == "nearctic_ID"){
      nearctic_mat_B <- nam
      names(nearctic_mat_B) <- check_list$nearctic_ID
    }
    if (names(check_list[i]) == "afrotropics_ID"){
      afrotropics_mat_B <- nam
      names(afrotropics_mat_B) <- check_list$afrotropics_ID
    }
    if (names(check_list[i]) == "neotropics_ID"){
      neotropics_mat_B <- nam
      names(neotropics_mat_B) <- check_list$neotropics_ID
    }
    if (names(check_list[i]) == "australasia_ID"){
      australasia_mat_B <- nam
      names(australasia_mat_B) <- check_list$australasia_ID
    }
  
  }
  list_matrix <- list(palearctic_mat_B, nearctic_mat_B, afrotropics_mat_B, neotropics_mat_B, australasia_mat_B)
  names(list_matrix) <- c("palearctic_mat_B","nearctic_mat_B","afrotropics_mat_B","neotropics_mat_B", "australasia_mat_B")
  return(list_matrix)
}

matrix_list_B <- list_matrix_B_function(ID_list)

# Calculate novelty and return output in a list
list_novelty_function <- function(matrix_list){
  
  for (i in 1:length(matrix_list)) {
    
    nam <- lapply(names(matrix_list[[i]][]), 
             function(ID){
               print(ID)
               site.sp.mat <- (matrix_list[[i]][ID])
               site.sp.mat <- site.sp.mat[[ID]]
               #site.sp.mat <- site.sp.mat[[ID]]
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
    
    # Select correct data, delete NA's
  
    if (names(matrix_list[i]) == "palearctic_mat_B"){
      palearctic_novelty_B <- nam
      names(palearctic_novelty_B) <- ID_list$palearctic_ID
      palearctic_novelty_B <- palearctic_novelty_B[!sapply(palearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "nearctic_mat_B"){
      nearctic_novelty_B <- nam
      names(nearctic_novelty_B) <- ID_list$nearctic_ID
      nearctic_novelty_B <- nearctic_novelty_B[!sapply(nearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "afrotropics_mat_B"){
      afrotropics_novelty_B <- nam
      names(afrotropics_novelty_B) <- ID_list$afrotropics_ID
      afrotropics_novelty_B <- afrotropics_novelty_B[!sapply(afrotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "neotropics_mat_B"){
      neotropics_novelty_B <- nam
      names(neotropics_novelty_B) <- ID_list$neotropics_ID
      neotropics_novelty_B <- neotropics_novelty_B[!sapply(neotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "australasia_mat_B"){
      australasia_novelty_B <- nam
      names(australasia_novelty_B) <- ID_list$australasia_ID
      australasia_novelty_B <- australasia_novelty_B[!sapply(australasia_novelty_B, function(x) all(is.na(x)))]
    }
    
    # Add a list of TimeSeries which will come in handy when evaluating the results.
    
    
  }
  list_novelty <- list(palearctic_novelty_B, nearctic_novelty_B, afrotropics_novelty_B, neotropics_novelty_B, australasia_novelty_B)
  names(list_novelty) <- c("palearctic_novelty_B","nearctic_novelty_B","afrotropics_novelty_B","neotropics_novelty_B", "australasia_novelty_B")
  return(list_novelty)
  
}

novelty_list_B <- list_novelty_function(matrix_list_B)

# Create a final master list for Binary results
Fish_Communities_B <- list(ID_list, matrix_list_B, novelty_list_B)
names(Fish_Communities_B) <- c("BioRealm_ID", "BioRealm_Matrices_B", "BioRealm_Novelty_B")
rm(matrix_list_B, novelty_list_B)



