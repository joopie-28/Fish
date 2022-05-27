### Hungary country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

hun_timeseries <- subset(time_series_data, 
                         Country == "HUN")$TimeSeriesID

hun_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% hun_timeseries)$Species))

colnames(hun_species) <- "Species"
hun_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly

hun_species$Status[2] <- "Native"
hun_species$Status[31] <- "Native"
hun_species$Status[39] <- "Native"


for (i in 1:nrow(hun_species)) {
  
  print(paste0(i, " out of ", nrow(hun_species)))
  
  if (hun_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(hun_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Hungary") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    hun_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (hun_species$Status[i] %in% c("native", "endemic")){
    hun_species$Status[i] <- "Native"
  }
  
  else {
    hun_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries



hun_country_level <- hun_species
rm(hun_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

hun_invaders <- read.csv("/Users/sassen/Desktop/hun_invaders.csv")

saveRDS(hun_invaders, "./Exotics_databases_countries/hun_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(hun_country_level)){
  if(hun_country_level$Species[i] %in% hun_invaders$Species){
    hun_country_level[i,] <- NA
  }
}
hun_country_level <- na.omit(hun_country_level)

# Create a list that treats all established species as native
hun_country_level_nn <- hun_country_level
hun_country_level_nn["Status"] <- "Native"

saveRDS(hun_country_level, "./Exotics_databases_countries/hun_country_level.rds")
saveRDS(hun_country_level_nn, "./Exotics_databases_countries/hun_country_level_nn.rds")

