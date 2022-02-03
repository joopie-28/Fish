# GBR country level species list

# Let's extract all the French species and label them using FishBase.

gbr_timeseries <- subset(time_series_data, 
                            Country == "GBR")$TimeSeriesID

gbr_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% gbr_timeseries)$Species))

colnames(gbr_species) <- "Species"
gbr_species$Status <- 0

# Connect to fishbase

gbr_species$Species[34] <- "Salvelinus alpinus" # subspecies
gbr_species$Species[37] <- "Chelon auratus" # liza aurata synonym

for (i in 1:nrow(gbr_species)) {
  
  print(paste0(i, " out of ", nrow(gbr_species)))
  
  if (gbr_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(gbr_species$Species[i],
                 server = "fishbase") %>%
      
      select(matches(c("country", "Status"))) %>%
      
      filter(country == "UK") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    gbr_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (gbr_species$Status[i] %in% c("native", "endemic")){
    gbr_species$Status[i] <- "Native"
  }
  
  else {
    gbr_species$Status[i] <- "Exotic"
  }
}

gbr_species$Species[34] <- "Salvelinus alpinus alpinus" # subspecies
gbr_species$Species[37] <- "Liza aurata" # liza aurata synonym

gbr_country_level <- gbr_species
rm(gbr_species)

# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

gbr_invaders <- read.csv("/Users/sassen/Desktop/gbr_invaders.csv")

saveRDS(gbr_invaders, "./Exotics_databases_countries/gbr_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(gbr_country_level)){
  if(gbr_country_level$Species[i] %in% gbr_invaders$Species){
    gbr_country_level[i,] <- NA
  }
}
gbr_country_level <- na.omit(gbr_country_level)

saveRDS(gbr_country_level, "./Exotics_databases_countries/gbr_country_level.rds")









