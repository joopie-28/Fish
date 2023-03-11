# Function for creating a data frame with variables extracted from the HydroAtlas

create_ENV_frame <- function(df, HYBAS_Level, HYBAS_scheme, variables_list){
  
  lay = ogrListLayers("/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb")[HYBAS_Level]
  
  # test
  HydroAtlas_lev <- st_read(dsn = "/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb", 
                             layer = lay)
  
  HydroAtlas_lev$HYBAS_ID <- as.factor(HydroAtlas_lev$HYBAS_ID)
  
  
  HYBAS_sub <- HYBAS_scheme[, c(1, which(colnames(HYBAS_scheme) == paste0('HYBAS_ID_', HYBAS_Level)))]
  HYBAS_sub[, 2] <- as.factor(HYBAS_sub[,2])
  class.frame <- df|> 
    rename('HYBAS_ID' = 'Basin') |> 
    left_join(HYBAS_sub, by = c("Site" = "TimeSeriesID")) |>
    dplyr::select(!"HYBAS_ID") |>
    rename("HYBAS_ID" = paste0('HYBAS_ID_', HYBAS_Level)) |>
    
    left_join(as.data.frame(HydroAtlas_lev[,c('HYBAS_ID',variables_list)]), by = 'HYBAS_ID')

  return(class.frame)
}


# Extract all sites that had enough data to apply novelty detection
sites <- as.character(unique(fullNovFrame_complete$site_ID))
sites_quartered <- as.character(unique(fullNovFrame_complete$site))

environ.df <- data.frame('ID' = sites_quartered, 
                         "Lat" = NA,
                         "Long" = NA)

# Cross-reference with survey data to obtain coordinates
index <- which(time_series_data$TimeSeriesID %in% sites)

# Find relevant data and convert to SF
geo.timeseries <- time_series_data[index, c("TimeSeriesID", "Longitude", "Latitude")]

# Account for seasonality

for (i in 1:length(sites_quartered)){
  temp_name <- strsplit(as.character(sites_quartered[i]), ".",
                        fixed = TRUE)[[1]][1]
  geo_index<-which(geo.timeseries == temp_name)
  environ.df$Lat[i] <- geo.timeseries$Latitude[geo_index]
  environ.df$Long[i] <- geo.timeseries$Longitude[geo_index]
  environ.df$Site[i] <- geo.timeseries$TimeSeriesID[geo_index]
  
}

# Spatial dataframe for modelling
geo.timeseries.sf <- st_as_sf(environ.df, coords = c("Long", "Lat"), crs = WG84) 

# Combine novelty data with environmental data in a spatial data frame.
geo.timeseries.full <- cbind(geo.timeseries.sf, 
                             rbindlist(lapply(geo.timeseries.sf$ID, 
                                              function(site){
                                                print(site)
                                                # Add up all the novelty metrics for binomial model
                                                back <- length(which(fullNovFrame_complete$site == site & 
                                                                       fullNovFrame_complete$cat == "back"))
                                                
                                                instant <- length(which(fullNovFrame_complete$site == site & 
                                                                          fullNovFrame_complete$cat == "instant"))
                                                
                                                cumul <- length(which(fullNovFrame_complete$site == site & 
                                                                        fullNovFrame_complete$cat == "cumul"))
                                                
                                                novel <- length(which(fullNovFrame_complete$site == site & 
                                                                        fullNovFrame_complete$cat == "novel"))
                                                
                                                # Include novelty classes based on persistence length for
                                                # model variation.
                                                if(novel > 0){
                                                  index <- (which(fullNovFrame_complete$site == site & 
                                                                    fullNovFrame_complete$cat == "novel"))[1]
                                                  print(index)
                                                  
                                                  novelty.class <- fullNovFrame_complete[index, "novel.class"]
                                                }
                                                else{
                                                  novelty.class <- "NONE"
                                                }
                                                
                                                # Add some more variables
                                                indices <- which(fullNovFrame_complete$site== site)
                                                
                                                
                                                Country <- fullNovFrame_complete$country.x[indices][1]
                                                BioRealm <- fullNovFrame_complete$BioRealm[indices][1]
                                                Basin <- fullNovFrame_complete$basin[indices][1]
                                                
                                                # Add in total length of timeseries as a covariate
                                                length <- back + instant + cumul +novel
                                                
                                                # Return a clean df for modelling
                                                df <- data.frame("back" =back, "instant" = instant, 
                                                                 "cumul" = cumul, "novel"= novel, "length" = length, "Class" = novelty.class,
                                                                 "Country" = Country, 'BioRealm' = BioRealm, 
                                                                 'Basin' = Basin)
                                                return(df)
                                              })))



# Create a frame that does things on the basin level ?
unique_basins<-unique(EnviroByTS_L12$HYBAS_ID)
for (i in 1:length(unique_basins)){
  sub_bas <- subset(EnviroByTS_L12, HYBAS_ID == as.character(unique_basins[i]))
  plot(sub_bas$Shape, add = T)
  plot(sub_bas$geometry, add = T, col = 'green',
       pch = 19, cex =0.1)
}

plot(ne_countries(scale = "Large", geounit = c('France'), type = 'map_unit',
                              returnclass = c("sp", "sf")))


# Create a frame that does things on the basin level ?
unique_basins<-unique(EnviroByTS_L12$HYBAS_ID)
basinLevel_df <- data.frame('HYBAS_ID' = unique_basins)

# loop through each basin and extract data
basinLevel_df$count_novels <- NA
for (i in 1:nrow(basinLevel_df)){
  indices <- which(EnviroByTS_L12$HYBAS_ID == basinLevel_df$HYBAS_ID[i])
  basinLevel_df$count_novels[]
  
  
  
  
}



