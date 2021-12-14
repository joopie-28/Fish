### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the "unaltered abundance" data pathway. Remember to load in functions, which are stored in separate files.

# Assign individual time series to BioRealm groups (if you have attempted to do the analysis using the presence/absence pathway
# earlier, this data will be loaded in already).

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

list_matrix_A_function <- function(check_list){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- abundance_identifier_function(TimeSeries_ID)
                    if(class(temp)=="matrix"){
                      # remove bins within sites with fewer species than cut off
                      temp <- temp[rowSums(temp) > rich.cutoff,]
                    }
                    return(temp)
                  }) 
    
    if (names(check_list[i]) == "palearctic_ID"){
      palearctic_mat_A <- nam
      names(palearctic_mat_A) <- check_list$palearctic_ID
    }
    if (names(check_list[i]) == "nearctic_ID"){
      nearctic_mat_A <- nam
      names(nearctic_mat_A) <- check_list$nearctic_ID
    }
    if (names(check_list[i]) == "afrotropics_ID"){
      afrotropics_mat_A <- nam
      names(afrotropics_mat_A) <- check_list$afrotropics_ID
    }
    if (names(check_list[i]) == "neotropics_ID"){
      neotropics_mat_A <- nam
      names(neotropics_mat_A) <- check_list$neotropics_ID
    }
    if (names(check_list[i]) == "australasia_ID"){
      australasia_mat_A <- nam
      names(australasia_mat_A) <- check_list$australasia_ID
    }
    
  }
  list_matrix <- list(palearctic_mat_A, nearctic_mat_A, afrotropics_mat_A, neotropics_mat_A, australasia_mat_A)
  names(list_matrix) <- c("palearctic_mat_A","nearctic_mat_A","afrotropics_mat_A","neotropics_mat_A", "australasia_mat_A")
  return(list_matrix)
}

matrix_list_A <- list_matrix_A_function(ID_list)

# Calculate novelty and return output in a list. No difference between Binary and Abundance data here except the similarity index used
# (Jaccard for binary, Bray-Curtis for abundance). Just remember to use the correct matrices as input (list_A for abundance, list_B for binary).

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
    
    # Select correct data, delete NA's
    
    if (names(matrix_list[i]) == "palearctic_mat_A"){
      palearctic_novelty_A <- nam
      names(palearctic_novelty_A) <- ID_list$palearctic_ID
      palearctic_novelty_A <- palearctic_novelty_A[!sapply(palearctic_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "nearctic_mat_A"){
      nearctic_novelty_A <- nam
      names(nearctic_novelty_A) <- ID_list$nearctic_ID
      nearctic_novelty_A <- nearctic_novelty_A[!sapply(nearctic_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "afrotropics_mat_A"){
      afrotropics_novelty_A <- nam
      names(afrotropics_novelty_A) <- ID_list$afrotropics_ID
      afrotropics_novelty_A <- afrotropics_novelty_A[!sapply(afrotropics_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "neotropics_mat_A"){
      neotropics_novelty_A <- nam
      names(neotropics_novelty_A) <- ID_list$neotropics_ID
      neotropics_novelty_A <- neotropics_novelty_A[!sapply(neotropics_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "australasia_mat_A"){
      australasia_novelty_A <- nam
      names(australasia_novelty_A) <- ID_list$australasia_ID
      australasia_novelty_A <- australasia_novelty_A[!sapply(australasia_novelty_A, function(x) all(is.na(x)))]
    }
    
    # Add a list of Time Series which will come in handy when evaluating the results.
    
    
  }
  list_novelty <- list(palearctic_novelty_A, nearctic_novelty_A, afrotropics_novelty_A, neotropics_novelty_A, australasia_novelty_A)
  names(list_novelty) <- c("palearctic_novelty_A","nearctic_novelty_A","afrotropics_novelty_A","neotropics_novelty_A", "australasia_novelty_A")
  return(list_novelty)
  
}

novelty_list_A <- list_novelty_function(matrix_list_A)

# Create a final master list for abundance results
Fish_Communities_A <- list(ID_list, matrix_list_A, novelty_list_A)
names(Fish_Communities_A) <- c("BioRealm_ID", "BioRealm_Matrices_A", "BioRealm_Novelty_A")
rm(matrix_list_A, novelty_list_A)
