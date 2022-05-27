### Colombia country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

col_timeseries <- subset(time_series_data, 
                         Country == "COL")$TimeSeriesID

col_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% col_timeseries)$Species))

colnames(col_species) <- "Species"
col_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly





for (i in 1:nrow(col_species)) {
  
  print(paste0(i, " out of ", nrow(col_species)))
  
  if (col_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(col_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Colombia") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    col_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (col_species$Status[i] %in% c("native", "endemic")){
    col_species$Status[i] <- "Native"
  }
  
  else {
    col_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

col_country_level <- col_species
rm(col_species)

# Create a list that treats all established species as native
col_country_level_nn <- col_country_level
col_country_level_nn["Status"] <- "Native"

saveRDS(col_country_level, "./Exotics_databases_countries/col_country_level.rds")
saveRDS(col_country_level_nn, "./Exotics_databases_countries/col_country_level_nn.rds")

