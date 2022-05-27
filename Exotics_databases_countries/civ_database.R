### Ivory coast country level species list ############
### Based on FishBase ############################
### Manual assignment of names where necessary ###


# Let's extract all the French species and label them using FishBase.

civ_timeseries <- subset(time_series_data, 
                         Country == "CIV")$TimeSeriesID

civ_species <- as.data.frame(unique(subset(Survey_Data, 
                                           TimeSeriesID %in% civ_timeseries)$Species))

colnames(civ_species) <- "Species"
civ_species$Status <- 0

# Fix up names that have changed or do not mathc fishbase exactly

civ_species$Species[32] <- "Coptodon zillii"
civ_species$Species[39] <- "Labeobarbus bynni"
civ_species$Species[41] <- "Enteromius trispilos"
civ_species$Status[42] <- "Native"
civ_species$Species[45] <- "Enteromius macrops"
civ_species$Status[47] <- "Native"
civ_species$Status[48] <- "Native"
civ_species$Status[53] <- "Native"
civ_species$Status[56] <- "Native"
civ_species$Status[57] <- "Native"
civ_species$Status[59] <- "Native"
civ_species$Species[61] <- "Coptodon dageti"
civ_species$Status[65] <- "Native"
civ_species$Status[70] <- "Native"
civ_species$Status[73] <- "Native"

for (i in 1:nrow(civ_species)) {
  
  print(paste0(i, " out of ", nrow(civ_species)))
  
  if (civ_species$Status[i] == 0){
    
    # Extract metrics from FishBase
    
    a <- country(civ_species$Species[i],
                 server = "fishbase") %>%
      
      dplyr::select(matches(c("country", "Status"))) %>%
      
      filter(country == "Cote d'Ivoire") 
    
    # Tidy up
    
    a <- unique(a[, c("Status", "country")])
    
    # Add to the species mainframe
    
    civ_species$Status[i] <- a$Status
    
  }
  else{
    next
  }
  
  # Now let's clean up the definitions, FishBase has a few categories
  # which are in essence the same for our purposes. 
  
  if (civ_species$Status[i] %in% c("native", "endemic")){
    civ_species$Status[i] <- "Native"
  }
  
  else {
    civ_species$Status[i] <- "Exotic"
  }
}

# Return these to synonyms used in timeseries

civ_species$Species[32] <- "Tilapia zillii"
civ_species$Species[39] <- "Barbus bynni"
civ_species$Species[41] <- "Barbus trispilos"
civ_species$Species[45] <- "Barbus macrops"
civ_species$Species[61] <- "Tilapia dageti"




civ_country_level <- civ_species
rm(civ_species)


# Because there is a difference between established non-natives and new invaders
# I am creating a separate status list for species that not established prior to 
# 1970. 

civ_invaders <- data.frame("Species" = "Oreochromis aureus", "Status" = "Invader")

saveRDS(civ_invaders, "./Exotics_databases_countries/civ_invaders.rds")


# I will now remove the species in the invader list from the complete list, 
# as they are invaders country wide and thus in every basin.

# Create a list that treats all established species as native
civ_country_level_nn <- civ_country_level
civ_country_level_nn["Status"] <- "Native"

saveRDS(civ_country_level, "./Exotics_databases_countries/civ_country_level.rds")
saveRDS(civ_country_level_nn, "./Exotics_databases_countries/civ_country_level_nn.rds")

