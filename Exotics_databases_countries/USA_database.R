## USA invader database ###

# Let's extract all the French species and label them using FishBase.

usa_timeseries <- subset(time_series_data, 
                            Country == "USA")$TimeSeriesID

usa_species <- as.data.frame(unique(subset(Survey_Data, 
                                              TimeSeriesID %in% usa_timeseries)$Species))

colnames(usa_species) <- "Species"
usa_species$Status <- 0

# 445 Species is quite a lot; lucky we can use fishbase for the majority.

usa_species$Species[6] <- "Notropis buccata" # synonym
usa_species$Species[10] <- "Moxostoma duquesnei" # spelling
usa_species$Species[102] <- "Esox americanus" # subspecies
usa_species$Species[138] <- "Opsopoeodus emiliae" # subspecies
usa_species$Species[171] <- "Esox americanus" # subspecies, both are native so OK, American pickerels
usa_species$Status[237] <- "Exotic" # Marine coral trout????
usa_species$Species[256] <- "Hysterocarpus traskii" # Subspecies
usa_species$Status[441] <- "Exotic" # Not recorded in USA per Fishbase 2007 but hails from Europe/Asia.


for (i in 1:nrow(usa_species)) {
  
  print(paste0(i, " out of ", nrow(usa_species)))
  
  if (usa_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(usa_species$Species[i],
                 server = "fishbase") %>%
      dplyr::select(matches(c("country", "Status"))) %>% filter(country == "USA") 
    
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

# Return to RIVfishtime names
usa_species$Species[6] <- "Ericymba buccata" # synonym
usa_species$Species[10] <- "Moxostoma duquesnii" # spelling
usa_species$Species[102] <- "Esox americanus americanus" # subspecies
usa_species$Species[138] <- "Opsopoeodus emiliae emiliae" # subspecies
usa_species$Species[171] <- "Esox americanus vermiculatus" # subspecies
usa_species$Species[256] <- "Hysterocarpus traskii traskii" # Subspecies

# Save list

usa_country_level <- usa_species
rm(usa_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

usa_invaders <- read.csv("/Users/sassen/Desktop/usa_invaders.csv")
usa_invaders[1,1] <- "Ctenopharyngodon idella"
usa_invaders[2,1] <- "Hypophthalmichthys nobilis"
usa_invaders[3,1] <- "Hypophthalmichthys molitrix"
usa_invaders[4,1] <- "Neogobius melanostomus"
usa_invaders[5,1] <- "Mylopharyngodon piceus"
usa_invaders[6,1] <- "Rutilus rutilus"


saveRDS(usa_invaders, "./Exotics_databases_countries/usa_invaders.rds")



# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.


saveRDS(usa_country_level, "./Exotics_databases_countries/usa_country_level.rds")
# Create a list that treats all established species as native
usa_country_level_nn <- usa_country_level
usa_country_level_nn["Status"] <- "Native"


saveRDS(usa_country_level_nn, "./Exotics_databases_countries/usa_country_level_nn.rds")



