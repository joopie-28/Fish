# Sweden database


# Let's extract all the French species and label them using FishBase.

swe_timeseries <- subset(time_series_data, 
                         Country == "SWE")$TimeSeriesID

swe_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% swe_timeseries)$Species))

colnames(swe_species) <- "Species"
swe_species$Status <- 0

# Connect to fishbase

swe_species$Species[22] <- "Salvelinus alpinus" # subspecies

for (i in 1:nrow(swe_species)) {
  
  print(paste0(i, " out of ", nrow(swe_species)))
  
  if (swe_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(swe_species$Species[i],
                 server = "fishbase") %>%
      
      select(matches(c("country", "Status"))) %>%
      
      filter(country == "Sweden") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    swe_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (swe_species$Status[i] %in% c("native", "endemic")){
    swe_species$Status[i] <- "Native"
  }
  
  else {
    swe_species$Status[i] <- "Exotic"
  }
}

swe_species$Species[22] <- "Salvelinus alpinus alpinus"

# No recent invaders in Sweden Time Series
saveRDS(swe_species, "./Exotics_databases_countries/swe_country_level.rds")

# Create a list that treats all established species as native
swe_country_level_nn <- swe_country_level
swe_country_level_nn["Status"] <- "Native"

saveRDS(swe_country_level_nn, "./Exotics_databases_countries/swe_country_level_nn.rds")


rm(swe_timeseries, swe_species, species_vector)


