# This is a shelved idea (11/03/2022).

# J.M. Sassen

# Novel Communities of Freshwater Fish

# Investigate spreading of novelty through geographical space
# using our data on Hydrobasins and novelty emergence.

# This might be a good area for future work, the code below 
# represents some early ideas and doodles on the subject.


###########################################
###### Step X. HydroBasin Aggregation #####
###########################################


hydrosheds <- read_sf(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")

hydrobasins_12 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev12_v1c")

hydrobasins_8 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev08_v1c")

hydrobasins_0 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev00_v1c")



# Attempt to match geometry

loc_points <- time_series_data[, c("Longitude", "Latitude", "TimeSeriesID")]

# Converts lat-long coordinates to geometry
loc_points_sf <- st_as_sf(loc_points, 
                          coords = c("Longitude", "Latitude"), 
                          crs = st_crs(hydrobasins_12))

# Set ID to character instead of integer
hydrobasins_12$HYBAS_ID <- as.character(hydrobasins_12$HYBAS_ID)


# Essential line, functions will not run without this setting 
sf::sf_use_s2(FALSE)

# Create a matrix containing HydroBasin ID and Tedesco Basin name
# This piece of code will find the intersection of the coordinates
# and the named basin.

loc_points <- loc_points_sf %>% mutate(
  
  intersection = as.integer(st_intersects(geometry, 
                                          hydrobasins_12)),
  
  hybas_id = if_else(is.na(intersection), '', 
                     hydrobasins_12$HYBAS_ID[intersection])
)

# Store timeseriesID - subbasin combo's in a 'by' list
loc_by_subbasin <- by(loc_points$TimeSeriesID, INDICES = loc_points$hybas_id, FUN = list)

# Filter out the French ones

france_ID <- subset(time_series_data, Country == "FRA")$TimeSeriesID

for (i in 1:length(loc_by_subbasin)){
  print(i)
  for(j in 1:length(loc_by_subbasin[[i]])){
    if(loc_by_subbasin[[i]][j] %notin% france_ID){
      loc_by_subbasin[[i]] <- NA
      
    }
  }
}

# Clean up which leaves the French ones

loc_by_subbasin <- loc_by_subbasin[!is.na(loc_by_subbasin)]


totals_inv <- read_csv("/Users/sassen/Desktop/Invaders_country.csv")
totals_inv$totals <- 0
totals_inv$totals[1] <- nrow(france_country_level_nn)
totals_inv$totals[2] <- nrow(usa_country_level_nn)
totals_inv$totals[3] <- nrow(gbr_country_level_nn)
totals_inv$totals[4] <- nrow(spain_country_level_nn)
totals_inv$totals[5] <- nrow(swe_country_level_nn)
totals_inv$totals[6] <- nrow(finland_country_level_nn)


ggplot(totals_inv) +
  geom_bar(aes(x = Country, y = totals), stat = "identity", colour = "green", fill = "green") +
  scale_color_grey() + theme_classic() +
  geom_bar(aes(x = Country, y = Invaders), stat = "identity", colour = "red", fill = "red") + 
  ylab("Unique species per country")




