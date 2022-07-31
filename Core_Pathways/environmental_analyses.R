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
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"

#### Step 0. Pre-process novel  community data ####

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

geo.timeseries.sf <- st_as_sf(geo.timeseries, coords = c("Longitude", "Latitude"), crs = WG84) 


#### Step 1. Extracting geographic data from multiple sources ####


# Elevation, from the WORLDCLIM database
elevation.data <- raster("/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_elev.tif")

# Built-up GHSL
builtup.data.2 <- raster("/Users/sassen/Desktop/Fish_GeoData/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-1975-30ss-v1-geotiff/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-1975-30ss-v1/ghsl-population-built-up-estimates-degree-urban-smod_built-30ss-1975.tif")
builtup.data.3 <- raster("/Users/sassen/Desktop/Fish_GeoData/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2014-2015-30ss-v1-geotiff/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2014-2015-30ss-v1/ghsl-population-built-up-estimates-degree-urban-smod_built-30ss-2014.tif")

# Land-use Globcover 2009
land.use.data <- raster("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER_L4_200901_200912_V2.3.tif")

# Also import codes legend for land use
land.use.legend <- read.csv("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER2009_legend.csv")

# Population GHSL 2015 at 250 M (Mollweide)
population.data.2 <- raster("/Users/sassen/Desktop/Fish_GeoData/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-1975-30ss-v1-geotiff/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-1975-30ss-v1/ghsl-population-built-up-estimates-degree-urban-smod_pop-30ss-1975.tif")
population.data.3 <- raster("/Users/sassen/Desktop/Fish_GeoData/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2014-2015-30ss-v1-geotiff/ghsl-population-built-up-estimates-degree-urban-smod-ghsl-2014-2015-30ss-v1/ghsl-population-built-up-estimates-degree-urban-smod_pop-30ss-2015.tif")

# Human modification
anthromod.data <- raster("/Users/sassen/Desktop/Fish_GeoData/lulc-human-modification-terrestrial-systems-geographic-geotiff/lulc-human-modification-terrestrial-systems_geographic.tif")

# Temperature from worldclim

path.69 <- "/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_tmax_1960-1969/wc2.1_2.5m_tmax_1969-"

path.18 <- "/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_tmax_2010-2018/wc2.1_2.5m_tmax_2018-"

path.69.min <- "/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_tmin_1960-1969/wc2.1_2.5m_tmin_1969-"

path.18.min # need to download

file.list <- list("01.tif", "02.tif","03.tif", "04.tif",
                  "05.tif", "06.tif","07.tif", "08.tif",
                  "09.tif", "10.tif","11.tif", "12.tif")

temp.loader <- function(path, file.list){
  
  test <- lapply(file.list, function(x){
  
   file.name <- paste0(path, x)
   
    raster.data <- raster(file.name)
  
    raster::extract(raster.data, geo.timeseries.sf)
  
  })
  k<-as.data.frame(do.call(cbind, test))
  
  k$mean<-(rowSums(k)/length(test))
  
  return(k$mean)
}

# Extract environmental data and store in a new data frame.
geo.timeseries.sf$Builtup.1975 <- raster::extract(builtup.data.2, geo.timeseries.sf)
geo.timeseries.sf$Builtup.2015 <- raster::extract(builtup.data.3, geo.timeseries.sf)
geo.timeseries.sf$Builtup.delta <- geo.timeseries.sf$Builtup.2015 - geo.timeseries.sf$Builtup.1975

# Population
geo.timeseries.sf$Pop.1975 <- raster::extract(population.data.2, geo.timeseries.sf)
geo.timeseries.sf$Pop.2015 <- raster::extract(population.data.3, geo.timeseries.sf)
geo.timeseries.sf$Pop.delta <- geo.timeseries.sf$Pop.2015 - geo.timeseries.sf$Pop.1975

# Geographical Data
geo.timeseries.sf$Lat <- geo.timeseries$Latitude
geo.timeseries.sf$Long <- geo.timeseries$Longitude
geo.timeseries.sf$Elevation <- raster::extract(elevation.data, geo.timeseries.sf)

# Anthropological impact index
geo.timeseries.sf$Anthromod <- raster::extract(anthromod.data, geo.timeseries.sf)

# Temperatures
geo.timeseries.sf$tmax.mean.1969 <- temp.loader(path.69, file.list)
geo.timeseries.sf$tmax.mean.2019 <- temp.loader(path.18, file.list)
geo.timeseries.sf$tmin.mean.1969 <- temp.loader(path.69.min, file.list)
geo.timeseries.sf$mean.delta <- geo.timeseries.sf$tmax.mean.2019 - geo.timeseries.sf$tmax.mean.1969

# These are classifications, not continuous data. Convert to factor.
geo.timeseries.sf$Land_use <- raster::extract(land.use.data, geo.timeseries.sf)

# Some categories are more interesting than others, i.e. we are not that concerned
# with the difference between an open forest or closed forest for our questions.
# Therefore, we will merge some groups which will allow for a more meaningful comparison.
# Here, I create a new df which holds info on original classes and new classes.

land.use.legend$New_Class <- c(rep("Cropland", 4), rep("Natural Vegetation", 14), "Urban and Artificial", "Bare area", "Water body", "Snow and Ice", "No data")

# Iterate to replace codes with category names
for(i in 1:nrow(geo.timeseries.sf)){
  
  index <- which(land.use.legend$Value == geo.timeseries.sf$Land_use[i])
  
  geo.timeseries.sf$Land_use[i] <- land.use.legend$New_Class[index]
  
}

# Set as factor
geo.timeseries.sf$Land_use <- as.factor(geo.timeseries.sf$Land_use)


# Combine novelty data with environmental data in a spatial data frame.
geo.timeseries.full <- cbind(geo.timeseries.sf, 
                             rbindlist(lapply(geo.timeseries.sf$TimeSeriesID, 
                                              function(site){
  # Add up all the novelty metrics for binomial model
  back <- length(which(full.novel.mat$site == site & 
                         full.novel.mat$cat == "back"))
  instant <- length(which(full.novel.mat$site == site & 
                            full.novel.mat$cat == "instant"))
  cumul <- length(which(full.novel.mat$site == site & 
                          full.novel.mat$cat == "cumul"))
  novel <- length(which(full.novel.mat$site == site & 
                          full.novel.mat$cat == "novel"))
  
  # Include novelty classes based on persistence length for
  # model variation.
  if(novel > 0){
    index <- (which(full.novel.mat$site == site & 
                    full.novel.mat$cat == "novel"))[1]
    print(index)
  
    novelty.class <- full.novel.mat[[index, "novel.class"]]
  }
  else{
    novelty.class <- "NONE"
  }
  
  # Add Mean invader increase in each time series
  indices <- which(full.novel.mat$site == site)
  
  mean.INV.inc<-sum(full.novel.mat[indices, "INC_increase"])/length(indices)
  
  # Add in total length of timeseries as a covariate
  length <- back + instant + cumul +novel
  
  # Return a clean df for modelling
  df <- data.frame("back" =back, "instant" = instant, 
                   "cumul" = cumul, "novel"= novel, "length" = length,
                   "Mean_Invader_Increase" =mean.INV.inc, "Class" = novelty.class)
  return(df)
})))




#### Step 2. Modelling novelty emergence as an effect of environment ####

geo.timeseries.full.copy <- geo.timeseries.full

geo.timeseries.full.copy$Mean_Invader_Increase <- scale(geo.timeseries.full.copy$Mean_Invader_Increase, center = T, scale = T)
geo.timeseries.full.copy$length <- geo.timeseries.full$length
geo.timeseries.full.copy$mean.delta <- scale(geo.timeseries.full.copy$mean.delta, center = T, scale = T)
geo.timeseries.full.copy$tmin.mean.1969 <- scale(geo.timeseries.full.copy$tmin.mean.1969, center = T, scale = T)
geo.timeseries.full.copy$Anthromod <- scale(geo.timeseries.full.copy$Anthromod, center = T, scale = T)
geo.timeseries.full.copy$Elevation <- scale(geo.timeseries.full.copy$Elevation, center = T, scale = T)

fit <- glm(cbind(novel, (length-novel))~scale(length, center = T, scale = T)+Anthromod+mean.delta+Land_use+Elevation+Mean_Invader_Increase, 
           data = geo.timeseries.full.copy, family = "binomial")

summary(fit)


fit.no.blip <- glm(cbind(novel, (length-novel))~scale(length, center = T, scale = T)+Elevation+
                     Anthromod+mean.delta+
                     Mean_Invader_Increase+Land_use, 
                   data = subset(geo.timeseries.full.copy ,Class !="BLIP" & Class != "END"), family = "binomial")


summary(fit.no.blip)






#### Step 3. Plot Model Results ####

plot_model(fit.no.blip, type = "pred",  terms=c("Anthromod [all]"), 
           title = "Probability of Novelty Emergence explained by Invader dynamics",
           axis.title = c("Net change in Invader Relative Abundance","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green", show.data = T)



#### End of Analyses ####











