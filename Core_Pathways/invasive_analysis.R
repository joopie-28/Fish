#### Analysis of Invasives ####
#### DOODLES AND IDEAS ########

library(rgdal)
library(raster)
library(rgeos)
library(sf)
library(tidyverse)
library(rfishbase)
library(stringi)

# Quick negation for use in loops
`%notin%` <- Negate(`%in%`)


# Import the invasive by basin data set, the column names are not configured well
# so we'll rename them.

# These contain native/nonnative data for species PER BASIN
invasives_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Occurrence_Table.csv")
colnames(invasives_data) <- c("Basin", "Species", "Status", "TSN_ITIS_Code", "FishBase_code", "Valid_name", "Sighting")

# These contain positional information for NAMED basins
basin_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Drainage_Basins_Table.csv")
colnames(basin_data) <- c("Basin", "Country", "BioRealm", "Endorheic", "Lon", "Lat", "Med_Lon", "Med_Lat", "Surface_Area")

# Load in the shapefile which contains basin names AND coordinate-based shapefiles (Tedesco et Al.)
occurence_shapefile <- read_sf(dsn = "/Users/sassen/Desktop/datatoFigshare/Basin042017_3119.shp")



# Our species names are formatted awkwardly so we need to fix that too.

invasives_data$Species <- gsub(invasives_data$Species, 
                               pattern ="\\.", 
                               replacement = " ")

invasives_data$Valid_name <- gsub(invasives_data$Valid_name, 
                                  pattern ="\\.", 
                                  replacement = " ")



# We are going to look at invasives per country, so let's make a subset of Hydro Basins for a 
# particular country. The invasive status for each species in French Basins

survey_ID <- "G8146"

index <- which(time_series_data$TimeSeriesID == survey_ID)

basin_name <- time_series_data$Basin_name[index]

invasive_status <- subset(invasives_data, 
                          Basin %in% subset(invasives_data, 
                                            Basin == basin_name)$Basin)

#### Link the basin names to their ID's ###

invasive_status$Hydro_ID <- NA


### Matching basin names and ID using Geolocation, so that we can determine invasives per basin ####


# We now need to allocate a basin name for every TimeSeries_ID in time_series_data,
# so that we can reliably tag our fish species.

# This function matches available names to HydroBasin ID's.
basin_name_code <- basin_name_match_function(occurence_shapefile)

# Match the basins to their ID.

for (i in 1:nrow(invasive_status)) {
  
  for (j in 1:nrow(basin_name_code)) {
  
    if (invasive_status$Basin[i] == basin_name_code$Basin_name[j]){
      
      invasive_status$Hydro_ID[i] <- basin_name_code$HydroBasin[j]
    }
  }
}

# Now we can determine which fish is invasive in each HydroBasin
# Unfortunately we can not use basins that have no ID...

invasive_status <- na.omit(invasive_status)

# As mentioned before, in the cases where there are species missing from the
# basin-specific database OR if the name for a basin is missing all together, 
# we will descend to the country level and use that database to assign 
# native/nonnative status.

# Let's play around with the matrices.

survey_ID <- "G7977"

# This function finds the basin of the survey, and first tries to assign
# native/nonnative on a basin scale. If there is no basin level data it
# plugs in to country level data. Only working for France right now.

france_country_level <- france_species_status

status_assignment_function <- function(survey_ID){
  
  test_matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[survey_ID]]
  
  if (typeof(test_matrix) == "character"){
    return("Not enough Data")
  }
  
  if (is.na(test_matrix)){
    return("Not enough Data")
  }

  # Find the basin name we are looking for
  
  index <- which(time_series_data$TimeSeriesID == survey_ID)
  
  basin_name <- time_series_data$Basin_name[index]
  
  invasive_status <- subset(invasives_data, 
                            Basin %in% subset(invasives_data, 
                                              Basin == basin_name)$Basin)
  
 
  
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
  
  # If there is no basin ID go straight to country level
  if (is.na(basin_name)){
    for (i in 1:nrow(species_vector)) {
      print("Allocating at country level")
      for (j in 1:nrow(france_country_level)) {
        if (species_vector$V1[i] == france_country_level$Species[j]) {
          species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", france_country_level$Status[j], " C")
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
    
      if (species_vector$V1[i] %in% france_country_level$Species) {
        for (j in 1:nrow(france_country_level)) {
          if (species_vector$V1[i] == france_country_level$Species[j]) {
            species_vector$V1[i] <- paste0(species_vector$V1[i], ", ", france_country_level$Status[j], " C")
          }
        }  
      }
    }
  }
  
  colnames(test_matrix) <- species_vector$V1
 
  return(test_matrix)

}

test <- status_assignment_function(survey_ID)


big_list <- as.list(subset(time_series_data, Country == "FRA")$TimeSeriesID)

# Run for all french timeseries, just a test
output <- lapply(big_list[], 
              function(TimeSeries_ID){
                print((TimeSeries_ID))
                temp <- status_assignment_function(TimeSeries_ID)
                  
                return(temp)
              })

for (i in 1:length(output)) {
  print(i)
  names(output)[i] <- big_list[[i]]
}

check <- identify.novel.gam(matrix,
                            metric = "bray",
                            site = "G8168",
                            alpha = 0.05)


# Compute contribution of natives
matrix_invasive <- matrix
matrix_invasive$NNC <- 0



for (i in 1:nrow(matrix_invasive)) {
  NNC <- 0
  for (j in 1:ncol(matrix_invasive)) {
    if(stri_detect_fixed(colnames(matrix_invasive)[j], "exotic")){
      NNC <- NNC + matrix_invasive[i,j]
    }
  
  }
  matrix_invasive$NNC[i] <- NNC/rowSums(matrix_invasive[i,])
  
}


gam(min.dist.tr ~ s(bins, bs="cr", k= set.k), 
    family=betar(),
    method="REML")



# Practice run with just one Hydrobasin
{ID <- 2080021030

newmap <- getMap(resolution = "high")

svg("map_Europe.svg", height=4, width=6)
par(mar=c(1, 1, 1, 1))
waterColor <- 'cadetblue1'
plot(newmap, xlim = c(-4.5, 8), ylim = c(44, 49), 
     asp = 1,lty=1, lwd=1,
     bg= "white", col="#ffffff") 

plot(flow_data, col= "grey" , add=T) # plot basins

# 2080021030

current <- subset(time_series_data, HydroBasin == 2080021030)
space <- current[,c("Longitude", "Latitude")]
points(x=space$Longitude, y=space$Latitude, pch= 16, col="red", cex=0.3)

# 2080022540
current_1 <- subset(time_series_data, HydroBasin == 2080022540)
space_1 <- current_1[,c("Longitude", "Latitude")]
points(x=space_1$Longitude, y=space_1$Latitude, pch= 16, col="blue", cex=0.3)

# 2080016510
current_2 <- subset(time_series_data, HydroBasin == 2080016510)
space_2 <- current_2[,c("Longitude", "Latitude")]
points(x=space_2$Longitude, y=space_2$Latitude, pch= 16, col="green", cex=0.3)

# 2080016680
current_3 <- subset(time_series_data, HydroBasin == 2080016680)
space_3 <- current_3[,c("Longitude", "Latitude")]
points(x=space_3$Longitude, y=space_3$Latitude, pch= 16, col="gold", cex=0.3)

# 2080020330
current_4 <- subset(time_series_data, HydroBasin == 2080020330)
space_4 <- current_4[,c("Longitude", "Latitude")]
points(x=space_4$Longitude, y=space_4$Latitude, pch= 16, col="gold", cex=0.3)
}










