#### Characterizing the environmental drivers of novelty emergence #######

############################################################################
### This script includes all main analyses linking geographic factors 
### such as land-use type, human population, elevation and climate to
### the emergence of novelty.

# J.M. Sassen 06-06-2022 

# Load required packages
package.loader(c("rgdal", "raster", "rgeos", "sf", 
                 "tidyverse", "rfishbase","mgcv",
                 "vegan", "lme4", "nlme", 
                 "DHARMa", "merTools", "shape",
                 "multcomp", "maptools", "sp", 
                 "divDyn", "plotrix", "raster",
                 "rgeos", "fun", "analogue",
                 "brms", "data.table", "data.table", 
                 "stringi", "tinytex", "knitr",
                 "sjPlot", "rworldmap", "ggeffects", 
                 "gridExtra"))

# Load CRS
WG84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  

### Step 0. Pre-process novel  ommunity data ###

# We have a couple of different goals here; we want to look at things
# like elevation, which do not change over time. However, we also wish
# to look at temperature, land-use, population, all og which have a temporal
# aspect. Thus, for some variables all we need is geographical data for the
# survey location, for others we want time-location information. Let's start
# with a simple object containing spatial information only.

# Extract all sites that had enough data to apply novelty detection
sites <- as.character(unique(full.novel.mat$site))

# Cross-reference with survey data to obtain coordinates

index <- which(time_series_data$TimeSeriesID %in% sites)

# Find relevant data and convert to SF
geo.timeseries <- time_series_data[index, c("TimeSeriesID", "Longitude", "Latitude")]

geo.timeseries.sf <- st_as_sf(geo.timeseries, coords = c("Longitude", "Latitude"), crs = WG84) # convert to sf


### Step 1. Extracting geographic data from multiple sources ###

# Elevation, from the WORLDCLIM database
elevation.data <- raster("./inputs/wc2.1_30s_elev.tif")

# Built-up score at 250 M (Mollweide)
builtup.data <- raster("./inputs/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0.tif")

# Land-use Globcover 2009
land.use.data <- raster("./inputs/GLOBCOVER_L4_200901_200912_V2.3.tif")

# Population GHSL 2015 at 250 M (Mollweide)

population.data <- raster("./inputs/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0.tif")

# Set the correct CRS where necessary

builtup.WG84 <- projectRaster(builtup.data, crs = WG84)

population.WG84 <- projectRaster(population.data, crs =WG84)
 
# Extract elevation data for our time series

geo.timeseries.sf$Elevation <- raster::extract(elevation.data, geo.timeseries.sf)
geo.timeseries.sf$Builtup <- raster::extract(builtup.WG84, geo.timeseries.sf)
geo.timeseries.sf$Land_use <- raster::extract(land.use.data, geo.timeseries.sf)
geo.timeseries.sf$Pop <- raster::extract(population.WG84, geo.timeseries.sf)

# Extract novelty information data and add to geodata
geo.timeseries.full <- cbind(geo.timeseries.sf, 
                             rbindlist(lapply(geo.timeseries.sf$TimeSeriesID, 
                                              function(site){

  back <- length(which(full.novel.mat$site == site & 
                         full.novel.mat$cat == "back"))
  instant <- length(which(full.novel.mat$site == site & 
                            full.novel.mat$cat == "instant"))
  cumul <- length(which(full.novel.mat$site == site & 
                          full.novel.mat$cat == "cumul"))
  novel <- length(which(full.novel.mat$site == site & 
                          full.novel.mat$cat == "novel"))
  # Add in total length
  length <- back + instant + cumul +novel
  
  df <- data.frame("back" =back, "instant" = instant, 
                   "cumul" = cumul, "novel"= novel, "length" = length)
  return(df)
})))

summary(glm(cbind(geo.timeseries.full$novel, (geo.timeseries.full$length-geo.timeseries.full$novel))~(length) + Land_use+ Builtup + Elevation, 
    data = geo.timeseries.full, family = "binomial"))



### Step 2. Extracting climate data and spatiotemporal manipulation ###














