### Botswana country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

bwa_timeseries <- subset(time_series_data, 
                         Country == "BWA")$TimeSeriesID

bwa_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% bwa_timeseries)$Species))

colnames(bwa_species) <- "Species"
bwa_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly

bwa_species$Status[1] <- "Native"
bwa_species$Status[11] <- "Native"
bwa_species$Status[12] <- "Native"
bwa_species$Species[13] <- "Enteromius paludinosus"
bwa_species$Status[24] <- "Native"
bwa_species$Status[25] <- "Native"
bwa_species$Species[31] <- "Enteromius unitaeniatus"
bwa_species$Species[32] <- "Enteromius radiatus"
bwa_species$Status[36] <- "Native"
bwa_species$Species[41] <- "Enteromius bifrenatus"
bwa_species$Species[43] <- "Enteromius trimaculatus"
bwa_species$Status[44] <- "Native"
bwa_species$Species[47] <- "Enteromius barnardi"

for (i in 1:nrow(bwa_species)) {
  
  print(paste0(i, " out of ", nrow(bwa_species)))
  
  if (bwa_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(bwa_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Botswana") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    bwa_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (bwa_species$Status[i] %in% c("native", "endemic")){
    bwa_species$Status[i] <- "Native"
  }
  
  else {
    bwa_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

bwa_species$Species[13] <- "Barbus paludinosus"
bwa_species$Species[31] <- "Barbus unitaeniatus"
bwa_species$Species[32] <- "Barbus radiatus"
bwa_species$Species[41] <- "Barbus bifrenatus"
bwa_species$Species[43] <- "Barbus trimaculatus"
bwa_species$Species[47] <- "Barbus barnardi"

bwa_country_level <- bwa_species
rm(bwa_species)







# Create a list that treats all established species as native
bwa_country_level_nn <- bwa_country_level
bwa_country_level_nn["Status"] <- "Native"

saveRDS(bwa_country_level, "./Exotics_databases_countries/bwa_country_level.rds")
saveRDS(bwa_country_level_nn, "./Exotics_databases_countries/bwa_country_level_nn.rds")

