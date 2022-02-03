

# Here, we have used the sf package to match the lat/long coordinates of
# Our known basins, which is data from RIVFishTime, to the basin names 
# Pedesco et Al. devised. It is not 100% perfect, as some small basins do not have names,
# but it works well for our study!

basin_name_match_function <- function(occurence_shapefile){
  
  loc_points <- time_series_data[, c("Longitude", "Latitude", "HydroBasin")]
  
  # Converts lat-long coordiantes to geometry
  loc_points_sf <- st_as_sf(loc_points, 
                            coords = c("Longitude", "Latitude"), 
                            crs = st_crs(occurence_shapefile))
  
  # Essential line, functions will not run without this setting 
  sf::sf_use_s2(FALSE)
  
  # Create a matrix containing HydroBasin ID and Tedesco Basin name
  # This piece of code will find the intersection of the coordinates
  # and the named basin.
  
  loc_points <- loc_points_sf %>% mutate(
    
    intersection = as.integer(st_intersects(geometry, 
                                            occurence_shapefile)),
    
    area = if_else(is.na(intersection), '', 
                   occurence_shapefile$BasinName[intersection])
  )
  
  # Tidy up the frame so we just have basin codes and names. Due to geometry
  # being a special object type we have to do this in a slightly convoluted way.
  # This returns a dataframe with HydroBasin ID and Tedesco hydrobasin name!
  
  basin_name_code <- as.data.frame(as.matrix(loc_points$HydroBasin))
  colnames(basin_name_code) <- "HydroBasin"
  basin_name_code$Basin_name <- loc_points$area
  
  # Some more tidying
  
  for (i in 1:nrow(basin_name_code)) {
    print(i)
    if(basin_name_code$Basin_name[i] == ""){
      
      basin_name_code$Basin_name[i] <- NA
    }
  }
  
  basin_name_code <- na.omit(basin_name_code)
  
  # This leaves us with 194 basins, out of 402 found in the RIVFishTime database..
  
  basin_name_code <- unique(basin_name_code)
  
  return(basin_name_code)
}