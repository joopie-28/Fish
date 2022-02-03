### France country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

france_timeseries <- subset(time_series_data, 
                            Country == "FRA")$TimeSeriesID

france_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% france_timeseries)$Species))

colnames(france_species) <- "Species"
france_species$Status <- 0

# Fix up one name that is not used anymore (Liza ramada), some species have also 
# only recently become apparent in France so I've filled those manually. This is
# tedious but frankly I can't think of any other way to do it!


france_species$Species[55] <- "Chelon ramada" # synonym
france_species$Species[58] <- "Alosa fallax" # synonym
france_species$Species[72] <- "Chelon auratus" # synonym
france_species$Status[80] <- "Exotic" # Not recorded in Fishbase
france_species$Status[62] <- "Exotic" # goby invasion
france_species$Status[63] <- "Exotic" # goby invasion
france_species$Status[64] <- "Exotic" # goby invasion
france_species$Status[69] <- "Exotic" # goby invasion
france_species$Species[81] <- "Salvelinus alpinus" # subspecies

for (i in 1:nrow(france_species)) {
  
  print(paste0(i, " out of ", nrow(france_species)))
  
  if (france_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(france_species$Species[i],
                 server = "fishbase") %>%
      
      select(matches(c("country", "Status"))) %>%
      
      filter(country == "France") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    france_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (france_species$Status[i] %in% c("native", "endemic")){
    france_species$Status[i] <- "Native"
  }
  
  else {
    france_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

france_species$Species[55] <- "Liza ramada" 
france_species$Species[58] <- "Alosa agone" 
france_species$Species[72] <- "Liza aurata"
france_species$Species[81] <- "Salvelinus alpinus alpinus" # subspecies

france_country_level <- france_species
rm(france_species)

# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

france_invaders <- read.csv("/Users/sassen/Desktop/france_invaders.csv")

saveRDS(france_invaders, "./Exotics_databases_countries/france_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(france_country_level)){
  if(france_country_level$Species[i] %in% france_invaders$Species){
    france_country_level[i,] <- NA
  }
}
france_country_level <- na.omit(france_country_level)

saveRDS(france_country_level, "./Exotics_databases_countries/france_country_level.rds")














