### Canada country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

aus_timeseries <- subset(time_series_data, 
                            Country == "AUS")$TimeSeriesID

aus_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% aus_timeseries)$Species))

colnames(aus_species) <- "Species"
aus_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly

aus_species$Species[11] <- "Percalates novemaculeatus"
aus_species$Species[16] <- "Anguilla australis"
aus_species$Status[45] <- "Exotic"
aus_species$Species[46] <- "Arrhamphus sclerolepis"
aus_species$Species[62] <- "Melanotaenia splendida"
aus_species$Status[69] <- "Exotic"


for (i in 1:nrow(aus_species)) {
  
  print(paste0(i, " out of ", nrow(aus_species)))
  
  if (aus_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(aus_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Australia") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    aus_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (aus_species$Status[i] %in% c("native", "endemic")){
    aus_species$Status[i] <- "Native"
  }
  
  else {
    aus_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

aus_species$Species[11] <- "Macquaria novemaculeata"
aus_species$Species[16] <- "Anguilla australis australis"
aus_species$Species[46] <- "Arrhamphus sclerolepis sclerolepis"
aus_species$Species[62] <- "Melanotaenia splendida splendida"


aus_country_level <- aus_species
rm(aus_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

aus_invaders <- read.csv("/Users/sassen/Desktop/aus_invaders.csv")

saveRDS(aus_invaders, "./Exotics_databases_countries/aus_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(aus_country_level)){
  if(aus_country_level$Species[i] %in% aus_invaders$Species){
    aus_country_level[i,] <- NA
  }
}
aus_country_level <- na.omit(aus_country_level)

# Create a list that treats all established species as native
aus_country_level_nn <- aus_country_level
aus_country_level_nn["Status"] <- "Native"

saveRDS(aus_country_level, "./Exotics_databases_countries/aus_country_level.rds")
saveRDS(aus_country_level_nn, "./Exotics_databases_countries/aus_country_level_nn.rds")













