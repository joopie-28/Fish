# Local origination script

# J.M. Sassen 30-03-2022

# Purpose of this script is to investigate the association between 
# primary origination species and novel community emergence; 
# the two functions here are made to do just that.

# Local origination function

loc.orig.function <- function(survey_identifier){
  
  test_matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[survey_identifier]]

  # We want to create a variable "new.species" 

  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)

  # Create a vector with years
  years <- sort(unique(sub_series$Year))

  # Apply over all years to find the species list for each year
  species_checklist <- lapply(years, FUN = function(x){
  
   checklist <- subset(sub_series, sub_series$Year == x)
   species_checklist <- matrix(unique(checklist[, c("Species")]), nrow = 1)
   return(species_checklist)
  
  })

  loc.orig.vector <- NA
  checker <- species_checklist

  for (i in 2:length(species_checklist)){
  
    new.sp <- length(checker[[i]][checker[[i]] %!in% checker[[i-1]]])
  
    loc.orig.vector <- c(loc.orig.vector, ifelse(length(new.sp) == 0, 0, new.sp))
  
    checker[[i]] <- unique(c(checker[[i]], checker[[i-1]]))
  
  }

  test_matrix$orig <- loc.orig.vector

  return(test_matrix)

}

# Local extinction function 

loc.ext.function <- function(survey_identifier){
  
  test_matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[survey_identifier]]
  
  # We want to create a variable "new.species" 
  
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  
  # Create a vector with years
  years <- sort(unique(sub_series$Year))
  
  # Apply over all years to find the species list for each year
  species_checklist <- lapply(years, FUN = function(x){
    
    checklist <- subset(sub_series, sub_series$Year == x)
    species_checklist <- matrix(unique(checklist[, c("Species")]))
    return(species_checklist)
    
  })
  
  loc.ext.vector <- NA
  checker <- species_checklist
  
  for (i in (length(species_checklist)-1):1){
    
    ext.sp <- length(checker[[i]][checker[[i]] %!in% checker[[i+1]]])
    
    loc.ext.vector <- c(loc.ext.vector, ifelse(length(ext.sp) == 0, 0, ext.sp))
    
    checker[[i]] <- unique(c(checker[[i]], checker[[i+1]]))
    
  }
  
  test_matrix$ext <- rev(loc.ext.vector)
  
  return(test_matrix)
  
}


















