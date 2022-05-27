### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

bra_timeseries <- subset(time_series_data, 
                         Country == "BRA")$TimeSeriesID

bra_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% bra_timeseries)$Species))

colnames(bra_species) <- "Species"
bra_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly
bra_species$Species[2] <- "Psalidodon paranae"
bra_species$Species[15] <- "Psalidodon fasciatus"
bra_species$Species[43] <- "Psalidodon anisitsi"
bra_species$Species[50] <- "Psalidodon schubarti"
bra_species$Species[53] <- "Coptodon rendalli"
bra_species$Species[99] <- "Megaleporinus obtusidens"
bra_species$Species[105] <- "Megaleporinus obtusidens"
bra_species$Species[114] <- "Megaleporinus macrocephalus"
bra_species$Species[120] <- "Cyphocharax naegelii"
bra_species$Species[125] <- "Megaleporinus piavussu"
bra_species$Species[130] <- "Myloplus tiete"
bra_species$Status[171] <- "Native"
bra_species$Species[173] <- "Diapoma itaimbe"
bra_species$Species[182] <- "Deuterodon luetkenii"
bra_species$Status[183] <- "Native"



for (i in 1:nrow(bra_species)) {
  
  print(paste0(i, " out of ", nrow(bra_species)))
  
  if (bra_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(bra_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Brazil") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    bra_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (bra_species$Status[i] %in% c("native", "endemic")){
    bra_species$Status[i] <- "Native"
  }
  
  else {
    bra_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

bra_species$Species[2] <- "Astyanax paranae"
bra_species$Species[15] <- "Astyanax fasciatus"
bra_species$Species[43] <- "Hyphessobrycon anisitsi"
bra_species$Species[50] <- "Astyanax schubarti"
bra_species$Species[53] <- "Tilapia rendalli"
bra_species$Species[99] <- "Leporinus obtusidens"
bra_species$Species[105] <- "Leporinus obtusidens"
bra_species$Species[114] <- "Leporinus macrocephalus"
bra_species$Species[120] <- "Cyphocharax nagelii"
bra_species$Species[125] <- "Leporinus piavussu"
bra_species$Species[130] <- "Myleus tiete"
bra_species$Species[173] <- "Cyanocharax itaimbe"
bra_species$Species[182] <- "Hyphessobrycon luetkenii"

bra_country_level <- bra_species
rm(bra_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 




bra_invaders <- read.csv("/Users/sassen/Desktop/bra_invaders.csv")

saveRDS(bra_invaders, "./Exotics_databases_countries/bra_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

for (i in 1:nrow(bra_country_level)){
  if(bra_country_level$Species[i] %in% bra_invaders$Species){
    bra_country_level[i,] <- NA
  }
}
bra_country_level <- na.omit(bra_country_level)

# Create a list that treats all established species as native
bra_country_level_nn <- bra_country_level
bra_country_level_nn["Status"] <- "Native"

saveRDS(bra_country_level, "./Exotics_databases_countries/bra_country_level.rds")
saveRDS(bra_country_level_nn, "./Exotics_databases_countries/bra_country_level_nn.rds")













