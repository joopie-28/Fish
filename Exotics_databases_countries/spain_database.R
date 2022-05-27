## Spain Database ###

# Let's extract all the French species and label them using FishBase.

spain_timeseries <- subset(time_series_data, 
                            Country == "ESP")$TimeSeriesID

spain_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% spain_timeseries)$Species))

colnames(spain_species) <- "Species"
spain_species$Status <- 0

# Fisbase connection 


spain_species$Status[21] <- "Exotic" # not in records yet
spain_species$Species[28] <- "Chelon ramada" # not in records yet
spain_species$Species[29] <- "Gasterosteus aculeatus" # synonym
spain_species$Status[30] <- "Native" # New species

for (i in 1:nrow(spain_species)) {
  
  print(paste0(i, " out of ", nrow(spain_species)))
  
  if (spain_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(spain_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Spain") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    spain_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (spain_species$Status[i] %in% c("native", "endemic")){
    spain_species$Status[i] <- "Native"
  }
  
  else {
    spain_species$Status[i] <- "Exotic"
  }
}

spain_species$Species[28] <- "Liza ramada" # not in records yet
spain_species$Species[29] <- "Gasterosteus gymnurus" # synonym

spain_country_level <- spain_species
rm(spain_species)

# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

spain_invaders <- read.csv("/Users/sassen/Desktop/spain_invaders.csv")

saveRDS(spain_invaders, "./Exotics_databases_countries/spain_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(spain_country_level)){
  if(spain_country_level$Species[i] %in% spain_invaders$Species){
    spain_country_level[i,] <- NA
  }
}

spain_country_level <- na.omit(spain_country_level)

# Create a list that treats all established species as native
spain_country_level_nn <- spain_country_level
spain_country_level_nn["Status"] <- "Native"



saveRDS(spain_country_level, "./Exotics_databases_countries/spain_country_level.rds")
saveRDS(spain_country_level_nn, "./Exotics_databases_countries/spain_country_level_nn.rds")















