### Canada country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

canada_timeseries <- subset(time_series_data, 
                            Country == "CAN")$TimeSeriesID

canada_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% canada_timeseries)$Species))

colnames(canada_species) <- "Species"
canada_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly




for (i in 1:nrow(canada_species)) {
  
  print(paste0(i, " out of ", nrow(canada_species)))
  
  if (canada_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(canada_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Canada") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    canada_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (canada_species$Status[i] %in% c("native", "endemic")){
    canada_species$Status[i] <- "Native"
  }
  
  else {
    canada_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries


canada_country_level <- canada_species
rm(canada_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

canada_invaders <- read.csv("/Users/sassen/Desktop/canada_invaders.csv")

saveRDS(canada_invaders, "./Exotics_databases_countries/canada_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(canada_country_level)){
  if(canada_country_level$Species[i] %in% canada_invaders$Species){
    canada_country_level[i,] <- NA
  }
}
canada_country_level <- na.omit(canada_country_level)

# Create a list that treats all established species as native
canada_country_level_nn <- canada_country_level[-40,]
canada_country_level_nn["Status"] <- "Native"

saveRDS(canada_country_level, "./Exotics_databases_countries/canada_country_level.rds")
saveRDS(canada_country_level_nn, "./Exotics_databases_countries/canada_country_level_nn.rds")













