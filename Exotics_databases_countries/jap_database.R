### Japan country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

jap_timeseries <- subset(time_series_data, 
                         Country == "JPN")$TimeSeriesID

jap_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% jap_timeseries)$Species))

colnames(jap_species) <- "Species"
jap_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly

jap_species$Species[2] <- "Oncorhynchus masou"
jap_species$Species[7] <- "Plecoglossus altivelis"
jap_species$Species[9] <- "Salvelinus leucomaenis"




for (i in 1:nrow(jap_species)) {
  
  print(paste0(i, " out of ", nrow(jap_species)))
  
  if (jap_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(jap_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Japan") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    jap_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (jap_species$Status[i] %in% c("native", "endemic")){
    jap_species$Status[i] <- "Native"
  }
  
  else {
    jap_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

jap_species$Species[2] <- "Oncorhynchus masou masou"
jap_species$Species[7] <- "Plecoglossus altivelis altivelis"
jap_species$Species[9] <- "Salvelinus leucomaenis leucomaenis"

jap_country_level <- jap_species
rm(jap_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

jap_invaders <- read.csv("/Users/sassen/Desktop/jap_invaders.csv")

saveRDS(jap_invaders, "./Exotics_databases_countries/jap_invaders.rds")



# Create a list that treats all established species as native
jap_country_level_nn <- jap_country_level
jap_country_level_nn["Status"] <- "Native"

saveRDS(jap_country_level, "./Exotics_databases_countries/jap_country_level.rds")
saveRDS(jap_country_level_nn, "./Exotics_databases_countries/jap_country_level_nn.rds")

