# Environmental drivers analysis #

# Joop Sassen 14-12-2022 #

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
                 "gridExtra", "ncdf4"))

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
sites <- as.character(unique(full.novel.mat.season$site_ID))
sites_quartered <- as.character(unique(full.novel.mat.season$site))

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


#### Step 1. Extracting geographic data from multiple sources ####

# Elevation, from the WORLDCLIM database
elevation.data <- raster("/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_elev.tif")

# Land-use Globcover 2009
land.use.data <- raster("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER_L4_200901_200912_V2.3.tif")

# Also import codes legend for land use
land.use.legend <- read.csv("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER2009_legend.csv")

# Human modification
anthromod.data <- raster("/Users/sassen/Desktop/Fish_GeoData/lulc-human-modification-terrestrial-systems-geographic-geotiff/lulc-human-modification-terrestrial-systems_geographic.tif")

# Anthromes
anthromes_2000.data <- raster("/Users/sassen/Desktop/Fish_GeoData/anthromes-v2-2000-global-geotif/a2000_global.tif")
anthromes_1900.data <- raster("/Users/sassen/Desktop/Fish_GeoData/anthromes-v2-1900-global-geotif/a1900_global.tif")

# Population
pop.data <- raster("/Users/sassen/Desktop/Fish_GeoData/gpw-v4-population-density-rev11_2015_1_deg_tif/gpw_v4_population_density_rev11_2015_1_deg.tif")

# Functions

anthro.impact.function <- function(geo.timeseries.sf, buffer_size){
  
  # Extract land use data from raster file (GlobCover)
  temp.output <- raster::extract(land.use.data, geo.timeseries.sf, buffer =buffer_size)
  
  # Calculate percentage of human land use in a 25 km area
  output <- lapply(temp.output, function(area){
    
    
    human.impact.pix <- length(which(area %in% c(11,14,20,190,200) ) )
    
    total.pix <- length(area)
    
    per.impact <- human.impact.pix/total.pix
    return(per.impact)
  })
  
  # Bind into data frame
  output<-as.data.frame.matrix(do.call(rbind, output))$V1
  
  
  return(output)
}

# Calculate percentage of land use as per Comte et al.
geo.timeseries.sf$Land_use_impact<-anthro.impact.function(geo.timeseries.sf, 
                                                          buffer_size = 10000)

# Combine novelty data with environmental data in a spatial data frame.
geo.timeseries.full <- cbind(geo.timeseries.sf, 
                             rbindlist(lapply(geo.timeseries.sf$Site, 
      function(site){
        # Add up all the novelty metrics for binomial model
        back <- length(which(full.novel.mat.season$site_ID == site & 
                               full.novel.mat.season$cat == "back"))
        
        instant <- length(which(full.novel.mat.season$site_ID == site & 
                                  full.novel.mat.season$cat == "instant"))
        
        cumul <- length(which(full.novel.mat.season$site_ID == site & 
                                full.novel.mat.season$cat == "cumul"))
        
        novel <- length(which(full.novel.mat.season$site_ID == site & 
                                full.novel.mat.season$cat == "novel"))
        
        
        
        # Include novelty classes based on persistence length for
        # model variation.
        if(novel > 0){
          index <- (which(full.novel.mat.season$site_ID == site & 
                            full.novel.mat.season$cat == "novel"))[1]
          print(index)
          
          novelty.class <- full.novel.mat.season[[index, "novel.class"]]
        }
        else{
          novelty.class <- "NONE"
        }
        
        # Add Mean invader increase in each time series
        indices <- which(full.novel.mat.season$site_ID== site)
        
        
        Country <- full.novel.mat.season$country[indices][1]
        
        # Add in total length of timeseries as a covariate
        length <- back + instant + cumul +novel
        
        # Return a clean df for modelling
        df <- data.frame("back" =back, "instant" = instant, 
                         "cumul" = cumul, "novel"= novel, "length" = length, "Class" = novelty.class,
                         "Country" = Country)
        return(df)
      })))

# Add aanthropogenic modification data
geo.timeseries.full$Anthromod <- raster::extract(anthromod.data, geo.timeseries.full)

# Add population
geo.timeseries.full$pop.2015 <- raster::extract(pop.data, geo.timeseries.full)

# Add antrhomes
geo.timeseries.full$Anthrome_2000<- raster::extract(anthromes_2000.data, geo.timeseries.full)
geo.timeseries.full$Anthrome_1900<- raster::extract(anthromes_1900.data, geo.timeseries.full)

# This could be a bit tidier
for(i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$Anthrome[i] %in% c(11,12)){
    geo.timeseries.full$Anthrome_factor[i] <- "Dense Settlement"
  }
  if(geo.timeseries.full$Anthrome[i] %in% c(21,22,23,24)){
    geo.timeseries.full$Anthrome_factor[i] <- "Villages"
  }
  if(geo.timeseries.full$Anthrome[i] %in% c(31,32,33,34)){
    geo.timeseries.full$Anthrome_factor[i] <- "Croplands"
  }
  if(geo.timeseries.full$Anthrome[i] %in% c(41,42,43)){
    geo.timeseries.full$Anthrome_factor[i] <- "Rangelands"
  }
  if(geo.timeseries.full$Anthrome[i] %in% c(51,52,53,54)){
    geo.timeseries.full$Anthrome_factor[i] <- "Seminatural"
  }
  if(geo.timeseries.full$Anthrome[i] %in% c(61,62)){
    geo.timeseries.full$Anthrome_factor[i] <- "Wildlands"
  }

}
for(i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(11,12)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Dense Settlement"
  }
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(21,22,23,24)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Villages"
  }
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(31,32,33,34)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Croplands"
  }
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(41,42,43)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Rangelands"
  }
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(51,52,53,54)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Seminatural"
  }
  if(geo.timeseries.full$Anthrome_1900[i] %in% c(61,62)){
    geo.timeseries.full$Anthrome_factor_1900[i] <- "Wildlands"
  }
  
}


geo.timeseries.full$Anthrome_factor <- as.factor(geo.timeseries.full$Anthrome_factor)
geo.timeseries.full$Anthrome_factor<-factor(geo.timeseries.full$Anthrome_factor, levels = c("Wildlands", "Croplands","Dense Settlement","Rangelands", "Seminatural", "Villages"))
geo.timeseries.full$Anthrome_factor_1900<-factor(geo.timeseries.full$Anthrome_factor_1900, levels = c("Wildlands", "Croplands","Dense Settlement","Rangelands", "Seminatural", "Villages"))

for(i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$Anthrome_factor_1900[i] != geo.timeseries.full$Anthrome_factor[i]){
    
    geo.timeseries.full$land_use_change[i] <- paste0(geo.timeseries.full$Anthrome_factor_1900[i], "-"
                                                     , geo.timeseries.full$Anthrome_factor[i])
  }
  else{
    geo.timeseries.full$land_use_change[i] <- "No Change"
  }
}
# Add a quarter column
geo.timeseries.full$Quarter <- NA

for (i in 1:nrow(geo.timeseries.full)){
 geo.timeseries.full$Quarter[i] <- strsplit(as.character(geo.timeseries.full$ID[i]), ".",
         fixed = TRUE)[[1]][2]
}



# Investigate effects

for (i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$novel[i] > 0){
    geo.timeseries.full$binary_novel[i] <- 1
  }
  else{
    geo.timeseries.full$binary_novel[i] <- 0
  }
  
}

no.blip.data <- subset(geo.timeseries.full, Class != "END" & Class != 'BLIP')

no.blip.data$pop.2000 <- as.numeric(scale(no.blip.data$pop.2000, scale = T, center = T))
no.blip.data$Anthromod <- as.numeric(scale(no.blip.data$Anthromod, scale = T, center = T))

mod <- (glm(binary_novel~length+pop.2000, data= no.blip.data, family = 'binomial'))

summary(mod)
visreg(mod, "cat", scale='response', partial=F, rug=F, 
       bty='l', xlab = 'Local Originations during Emergence', ylab = 'Persistence Probability of Novel Community',
       points = list(col='black'), line = list(col = 'gold'), axes=T)

mod<-glm(novel~INC +bin_lag+position, data = full.nov.no.lag, family = 'binomial')

summary(mod)

# Diversity
mod<-glm(cbind(orig, diversity)~ novel.class , data = full.nov.no.lag,family = 'binomial')

summary(mod)


points((novel)~(INC), data=full.novel.mat.season, col=alpha(colour = 'black', alpha = 0.4), pch = 19)

boxplot(Anthromod~binary_novel, data=no.blip.data)
points(mean(no.blip.data$Anthromod, na.rm=T))


plot(no.blip.data$pop.2000, log(no.blip.data$Anthromod+.001), "p")
