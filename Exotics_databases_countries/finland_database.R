# Finland species database #

# Let's extract all the French species and label them using FishBase.

finland_timeseries <- subset(time_series_data, 
                            Country == "FIN")$TimeSeriesID

finland_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% finland_timeseries)$Species))

colnames(finland_species) <- "Species"
finland_species$Status <- 0


# Allocate using Fishbase

for (i in 1:nrow(finland_species)) {
  
  print(paste0(i, " out of ", nrow(finland_species)))
  
  if (finland_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(finland_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Finland") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    finland_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (finland_species$Status[i] %in% c("native", "endemic")){
    finland_species$Status[i] <- "Native"
  }
  
  else {
    finland_species$Status[i] <- "Exotic"
  }
}

finland_country_level <- finland_species
rm(finland_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

finland_invaders <- read.csv("/Users/sassen/Desktop/finland_invaders.csv")

saveRDS(finland_invaders, "./Exotics_databases_countries/finland_invaders.rds")




# Create a list that treats all established species as native
finland_country_level_nn <- finland_country_level
# All "exotics" have actually been in finland for a long time
finland_country_level_nn["Status"] <- "Native"


saveRDS(finland_country_level_nn, "./Exotics_databases_countries/finland_country_level_nn.rds")
saveRDS(finland_country_level, "./Exotics_databases_countries/finland_country_level.rds")
