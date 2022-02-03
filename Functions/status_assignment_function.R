# This function assigns invasive status to fish based on both basin
# and country level data.

status_assignment_function <- function(survey_ID, country){
  
  if(country == "USA"){
    test_matrix <- Fish_Communities_A$BioRealm_Matrices_A_2$nearctic_mat_A[[survey_ID]]
  }
  else{
  test_matrix <- Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A[[survey_ID]]
  }
  
  if (typeof(test_matrix) == "character"){
    return(NA)
  }
  
  if (is.na(test_matrix)){
    return(NA)
  }
  
  # Find the basin name we are looking for
  
  index <- which(time_series_data$TimeSeriesID == survey_ID)
  
  basin_name <- time_series_data$Basin_name[index]
  
  invasive_status <- subset(invasives_data, 
                            Basin %in% subset(invasives_data, 
                                              Basin == basin_name)$Basin)
  
  if(country == "FRA"){
    country_level <- france_country_level
    invaders <- france_invaders
  }
  
  if(country == "GBR"){
    country_level <- gbr_country_level
    invaders <- gbr_invaders
  }
  
  if(country == "SWE"){
    country_level <- swe_country_level
    
  }
  

  ### Matching basin names and ID using Geolocation, so that we can determine invasives per basin ####
  
  # Match the basin to its ID.
  
  if(!is.na(basin_name)){
    
    #### Link the basin names to their ID's ###
    invasive_status$Hydro_ID <- NA
    
    for (i in 1:nrow(invasive_status)) {
      
      for (j in 1:nrow(basin_name_code)) {
        
        if (invasive_status$Basin[i] == basin_name_code$Basin_name[j]){
          
          invasive_status$Hydro_ID[i] <- basin_name_code$HydroBasin[j]
        }
      }
    }
  }
  
  #  Extract species names per basin
  species_vector <- as.data.frame(as.matrix(colnames(test_matrix)))
  
  # We will start with the "unestablished invaders" as this category takes
  # priority, this is always at a country level.
  
  if(country != "SWE"){
  
    for (i in 1:nrow(species_vector)) {
      print("Allocating invaders")
      for (j in 1:nrow(invaders)) {
        if (species_vector$V1[i] == invaders$Species[j]) {
          species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", invaders$Status[j], " C")
        }
      }  
    }  
  }
  
  
  
  # If there is no basin ID go straight to country level
  if (is.na(basin_name)){
    for (i in 1:nrow(species_vector)) {
      print("Allocating at country level")
      for (j in 1:nrow(country_level)) {
        if (species_vector$V1[i] == country_level$Species[j]) {
          species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", country_level$Status[j], " C")
        }
      }  
    }  
  }
  
  
  # Now that the ones that could not be done on basin level have been filled in
  # fill in the others on basin level data.
  else {
    print("Allocating at basin level")
    for (i in 1:nrow(species_vector)) {
      for (j in 1:nrow(invasive_status)){
        
        if ((species_vector$V1[i] %in% invasive_status$Species)){
          
          if(species_vector$V1[i] == invasive_status$Species[j]) {
            species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", invasive_status$Status[j], " B")
          }
        }
      }
    }
    
    
    # Now we have to clean up the scraps; there will be some fish that still do not have a status
    # because tedesco did not include them in his database. We thus run the country level again.
    for (i in 1:nrow(species_vector)) {
      print("Reached basin level but not included by Tedesco")
      
      if (species_vector$V1[i] %in% country_level$Species) {
        for (j in 1:nrow(country_level)) {
          if (species_vector$V1[i] == country_level$Species[j]) {
            species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", country_level$Status[j], " C")
          }
        }  
      }
    }
  }
  
  colnames(test_matrix) <- species_vector$V1
  
  return(test_matrix)
  
}
