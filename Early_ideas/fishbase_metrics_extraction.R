### Extracting interesting metrics for fish species in the database ###

fish_species <- as.list(unique(Survey_Data$Species))

# Function for extracting metrics from FishBase, for
# a list of species.

fishbase_metrics_function <- function(species_list){
  
  fish_diet_list <- lapply(species_list[], 
                           function(fish_name){
                             print(fish_name)
                             
                             value <- estimate(fish_name) %>% 
                               select(dplyr::matches("Troph"))
                             
                             troph <- round(value$Troph, digits = 3)
                             
                             return(troph)
                             
                           })
  
  # Give names to list
  for (i in 1:length(fish_diet_list)) {
    
    names(fish_diet_list)[i] <- fish_species[i]
    
  }
  
  # Turn into dataframe
  fish_metrics <- (as.data.frame(fish_diet_list))
  
  fish_features <- transpose(fish_metrics)
  
  # Get row and column names in order
  colnames(fish_features) <- "Trophic_level"
  rownames(fish_features) <- colnames(fish_metrics)
  
  rownames(fish_features) <- gsub(rownames(fish_features), 
                                  pattern ="\\.", 
                                  replacement = " ")
  
  return(fish_features)
  
}



# Experimenting with implementing these metrics into the timeseries database
# so that we can do some analysis.

# Select one matrix to work on

matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A$G106

matrix$average_trophic <- NA
matrix["Trophic_level",] <- NA

# Insert trophic values into matrix.

for (i in 1:(ncol(matrix)-1)){
  
  for (j in 1:nrow(fish_features)) {
    
    if (colnames(matrix)[i] == rownames(fish_features)[j]){
      
      matrix["Trophic_level",i] <- fish_features$Trophic_level[j]
      
    }
  }
}

# Compute average trophic level for each community.

for (i in 1:(nrow(matrix)-1)) {
  
  temp <- matrix[i,]*matrix["Trophic_level",]
  trophic_sum <- rowSums(temp[1, 1:(ncol(temp)-1)])
  total <- rowSums(matrix[i, 1:(ncol(temp)-1)])
  
  matrix[i,"average_trophic"] <- trophic_sum/total
  
}