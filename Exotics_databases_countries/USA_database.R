## USA invader database ###

# Let's extract all the French species and label them using FishBase.

usa_timeseries <- subset(time_series_data, 
                            Country == "USA")$TimeSeriesID

usa_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% usa_timeseries)$Species))

colnames(usa_species) <- "Species"
usa_species$Status <- 0

# 445 Species is quite a lot; lucky we can use fishbase for the majority.

for (i in 1:nrow(usa_species)) {
  
  print(paste0(i, " out of ", nrow(usa_species)))
  
  if (usa_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(usa_species$Species[i],
                 server = "fishbase") %>%
      select(matches(c("country", "Status"))) %>% filter(country == "USA") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    usa_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (usa_species$Status[i] %in% c("native", "endemic")){
    usa_species$Status[i] <- "Native"
  }
  
  else {
    usa_species$Status[i] <- "Exotic"
  }
}
