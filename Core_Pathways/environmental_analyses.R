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

path.69 <- "/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_tmax_1960-1969/wc2.1_2.5m_tmax_"

path.18 <- "/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_tmax_2010-2018/wc2.1_2.5m_tmax_"


year.list.69 <- list("1961-", "1962-", "1963-", "1964-", "1965-",
                  "1966-", "1967-", "1968-", "1969-")

year.list.18 <- list("2010-", "2011-", "2012-", "2013-", "2014-",
                     "2015-", "2016-", "2017-", "2018-")


month.list <- list("01.tif", "02.tif","03.tif", "04.tif",
                  "05.tif", "06.tif","07.tif", "08.tif",
                  "09.tif", "10.tif","11.tif", "12.tif")

temp.loader <- function(path, year.list, month.list){
  
  test <- lapply(year.list, function(year){
    
    lapply(month.list, function(month){
      file.name <- paste0(path, year, month)
      print(file.name)
      raster.data <- raster(file.name)
  
      raster::extract(raster.data, geo.timeseries.sf)
    })
  })
  
  k.2<-lapply(test, function(index){
      
    k<-as.data.frame(do.call(cbind, index))
    k$mean <- (rowSums(k)/ncol(k))
      
  })
      
  k.3<-as.data.frame(do.call(cbind, k.2))
    
  k.3$mean <- (rowSums(k.3)/ncol(k.3))
  
  return(k.3$mean)
}

#Test
test <- temp.loader(path.69, year.list.69, month.list)

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
geo.timeseries.sf$tmax.mean.1969 <- temp.loader(path.69, year.list.69, month.list)
geo.timeseries.sf$tmax.mean.2019 <- temp.loader(path.18, year.list.18, month.list)

geo.timeseries.sf$mean.delta <- geo.timeseries.sf$tmax.mean.2019 - geo.timeseries.sf$tmax.mean.1969


# Land use data in simple categories
geo.timeseries.sf$Land_use <- raster::extract(land.use.data, geo.timeseries.sf)

# Function for extracting adjacent values in land use data.
adj.cell.finder <- function(geo.timeseries.sf, land.use.data, directions = 8){
  
  # Filter those entries that are 'water body' as these are the ones
  # we want to convert to a terrestrial classification.
  
  filtered.df <- subset(geo.timeseries.sf, Land_use == 210)
  
  # Reduce it to just the necessary values for now
  
  filtered.df <- filtered.df[, c("TimeSeriesID", "geometry", "Land_use")]
  
  # Extract cell numbers using land use data
  filtered.df[c("cells", "GLOBCOVER_INDEX")] <- as.data.frame(raster::extract(land.use.data, 
                                                                              filtered.df, cellnumbers = T))
  
  # Re-extract values from adjacent cells with the chosen number of directions. def = 8
  cell.freq <-as.data.frame(raster::adjacent(land.use.data, 
                                             filtered.df$cells, directions = directions, sorted = T))
  
  # Assign values from cell indices 
  cell.freq$code <- land.use.data[cell.freq$to]
  
  # Format the cell frequencies as a tidy data frame
  freq.tb <- as.data.frame.matrix(table(cell.freq[,1], cell.freq[,3]))
  
  # Materialize rownames as a column for merging
  freq.tb$cells <- as.numeric(rownames(freq.tb))
  
  # Merge the two data frame according to cell index!
  merged.df <- merge(filtered.df, freq.tb ,by=c("cells"))
  
  
  # Nice! Now we know most of the adjacent cells will be water,
  # but that's ok. We are just looking for the most common terrestrial
  # classification. Need some code that computes this for us:
  
  # The water code is 210, let's assign that.
  water <- "210"
  
  col.index<-which(colnames(merged.df) == water)
  
  # Remove the water index
  merged.df <- merged.df[,-col.index]
  merged.df$new.code <- NA
  for (i in 1:nrow(merged.df)){
    a<-names(which.max(merged.df[i,5:ncol(merged.df)]))
    merged.df$new.code[i] <- a
    
  }
  
  return(merged.df)
}



# Calculate percentage of land use as per Comte et al



anthro.impact.function <- function(geo.timeseries.sf, buffer_size){

  # Extract land use data from raster file (GlobCover)
  temp.output <- raster::extract(land.use.data, geo.timeseries.sf, buffer = 25000)
  
  # Calculate percentage of human land use in a 25 km area
  output <- lapply(temp.output, function(area){
    
    print(area)
    human.impact.pix <- length(which(area %in% c(11,14,20,190,200) ) )
    
    total.pix <- length(area)
    
    per.impact <- human.impact.pix/total.pix
    return(per.impact)
  })
  
  # Bind into data frame
  output<-as.data.frame.matrix(do.call(rbind, output))
  
  
  return(output)
}
  
  
# Buffer size in meters 
  
geo.timeseries.sf$Land_use_impact<-anthro.impact.function(geo.timeseries.sf, buffer_size = 25000)





# Use output to assign new values

geo.timeseries.sf$Land_use_mod <- NA
adj.land.use<-adj.cell.finder(geo.timeseries.sf, land.use.data, directions = 8)

# Some categories are more interesting than others, i.e. we are not that concerned
# with the difference between an open forest or closed forest for our questions.
# Therefore, we will merge some groups which will allow for a more meaningful comparison.
# Here, I create a new df which holds info on original classes and new classes.

land.use.legend$New_Class <- c(rep("Cropland", 4), rep("Natural Vegetation", 14), "Urban and Artificial", "Bare area", "Water body", "Snow and Ice", "No data")

# Generate a modified land-use column
for(i in 1:nrow(geo.timeseries.sf)){
  
  if(geo.timeseries.sf$Land_use[i] == 210){
    
    for (j in 1:nrow(adj.land.use)){
      
      index <- which(geo.timeseries.sf$TimeSeriesID == adj.land.use$TimeSeriesID[j])
      
      geo.timeseries.sf$Land_use_mod[index] <- adj.land.use$new.code[j]
      
    }
    
  }else{
    
    geo.timeseries.sf$Land_use_mod[i] <- geo.timeseries.sf$Land_use[i]
    
  }
  
  
}

# Iterate to replace codes with category names
for(i in 1:nrow(geo.timeseries.sf)){
  
  index <- which(land.use.legend$Value == geo.timeseries.sf$Land_use_mod[i])
  
  geo.timeseries.sf$Land_use_mod[i] <- land.use.legend$New_Class[index]
  
}

# Set as factor
geo.timeseries.sf$Land_use_mod <- as.factor(geo.timeseries.sf$Land_use_mod)

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
  Country <- full.novel.mat$country[indices][1]
  
  # Add in total length of timeseries as a covariate
  length <- back + instant + cumul +novel
  
  # Return a clean df for modelling
  df <- data.frame("back" =back, "instant" = instant, 
                   "cumul" = cumul, "novel"= novel, "length" = length,
                   "Mean_Invader_Increase" =mean.INV.inc, "Class" = novelty.class,
                   "Country" = Country)
  return(df)
})))




#### Step 2. Modelling novelty emergence as an effect of environment ####

geo.timeseries.full.copy <- geo.timeseries.full

geo.timeseries.full.copy$Mean_Invader_Increase <- scale(geo.timeseries.full.copy$Mean_Invader_Increase, center = T, scale = T)

geo.timeseries.full.copy$mean.delta <- scale(geo.timeseries.full.copy$mean.delta, center = T, scale = T)
geo.timeseries.full.copy$tmin.mean.1969 <- scale(geo.timeseries.full.copy$tmin.mean.1969, center = T, scale = T)
geo.timeseries.full.copy$Anthromod <- scale(geo.timeseries.full.copy$Anthromod, center = T, scale = T)
geo.timeseries.full.copy$Elevation <- scale(geo.timeseries.full.copy$Elevation, center = T, scale = T)
geo.timeseries.full.copy$length <- geo.timeseries.full$length
geo.timeseries.full.copy$Land_use_impact <- scale(geo.timeseries.full.copy$Land_use_impact, center = T, scale = T)

# Add a binary response

geo.timeseries.full.copy$binary <- 0

for (i in 1:nrow(geo.timeseries.full.copy)){
  
  if(geo.timeseries.full.copy$novel[i] > 0){
  
    geo.timeseries.full.copy$binary[i] <- 1
  }
}

fit <- glm(cbind(novel, (length-novel))~
             Land_use_impact+mean.delta+Mean_Invader_Increase, 
           data = subset(geo.timeseries.full.copy, Class == "BLIP" | Class == "NONE"), family = "binomial")

fit <- glm(binary~length+
             Land_use_impact+
             Mean_Invader_Increase+mean.delta, 
           data = subset(geo.timeseries.full.copy, Class == "BLIP" | Class == "NONE"), family = "binomial")

no.outliers <- subset(geo.timeseries.full.copy, Country %in% c("SWE", "FRA", "GBR", "AUS", "USA"))
summary(fit)

# Streamflow
fit.no.blip <- glm(cbind(novel, (length-novel))~ Land_use_impact+mean.delta+Mean_Invader_Increase, 
                   data = subset(geo.timeseries.full.copy , Class !="BLIP" & Class != "END"), family = "binomial")

fit.no.blip <- glm(binary~length+mean.delta+ Land_use_impact+Mean_Invader_Increase , 
                   data = subset(geo.timeseries.full.copy, Class !="BLIP" & Class != "END"), family = "binomial")


summary(fit.no.blip)



# Turn land use into a continuous variable as per compte!


# Let's tidy this up and create an output csv.




#### Step 3. Plot Model Results ####
 
plot_model(fit.no.blip, type = "pred",  terms=c("Land_use_impact [all]"), 
           title = "Probability of Novelty Emergence explained x",
           axis.title = c("Net change in x","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green", show.data = F)



#### End of Analyses ####




# Function for extracting adjacent values in land use data.

adj.cell.finder <- function(geo.timeseries.sf, land.use.data, directions = 8){
  
  # Filter those entries that are 'water body' as these are the ones
  # we want to convert to a terrestrial classification.
  
  filtered.df <- subset(geo.timeseries.sf, Land_use == "Water body")
  
  # Reduce it to just the necessary values for now
  
  filtered.df <- filtered.df[, c("TimeSeriesID", "geometry", "Land_use")]
  
  # Extract cell numbers using land use data
  filtered.df[c("cells", "GLOBCOVER_INDEX")] <- as.data.frame(raster::extract(land.use.data, 
                                         filtered.df, cellnumbers = T))
  
  # Re-extract values from adjacent cells with the chosen number of directions. def = 8
  cell.freq <-as.data.frame(raster::adjacent(land.use.data, 
                                    filtered.df$cells, directions = directions, sorted = T))
  
  # Assign values from cell indices 
  cell.freq$code <- land.use.data[cell.freq$to]
  
  # Format the cell frequencies as a tidy data frame
  freq.tb <- as.data.frame.matrix(table(cell.freq[,1], cell.freq[,3]))
  
  # Materialize rownames as a column for merging
  freq.tb$cells <- as.numeric(rownames(freq.tb))
  
  # Merge the two data frame according to cell index!
  merged.df <- merge(filtered.df, freq.tb ,by=c("cells"))
  
  
  # Nice! Now we know most of the adjacent cells will be water,
  # but that's ok. We are just looking for the most common terrestrial
  # classification. Need some code that computes this for us:
  
  # The water code is 210, let's assign that.
  water <- "210"
 
  col.index<-which(colnames(merged.df) == water)
  
  # Remove the water index
  merged.df <- merged.df[,-col.index]
  merged.df$new.code <- NA
  for (i in 1:nrow(merged.df)){
    a<-names(which.max(merged.df[i,5:ncol(merged.df)]))
    merged.df$new.code[i] <- a
    
  }
  
  return(merged.df)
}

# Use output to assign new values

geo.timeseries.sf$Land_use_mod <- NA
adj.land.use<-adj.cell.finder(geo.timeseries.sf, land.use.data,directions = 16)

# Generate a modified land-use column
for(i in 1:nrow(geo.timeseries.sf)){
  
  if(geo.timeseries.sf$Land_use[i] == 210){
    
    for (j in 1:nrow(adj.land.use)){
      
      index <- which(geo.timeseries.sf$TimeSeriesID == adj.land.use$TimeSeriesID[j])
      
      geo.timeseries.sf$Land_use_mod[index] <- adj.land.use$new.code[j]
      
    }
    
  }else{
    
    geo.timeseries.sf$Land_use_mod[i] <- geo.timeseries.sf$Land_use[i]
    
  }
  
  
}

# Iterate to replace codes with category names
for(i in 1:nrow(geo.timeseries.sf)){
  
  index <- which(land.use.legend$Value == geo.timeseries.sf$Land_use_mod[i])
  
  geo.timeseries.sf$Land_use_mod[i] <- land.use.legend$New_Class[index]
  
}

# Set as factor
geo.timeseries.sf$Land_use_mod <- as.factor(geo.timeseries.sf$Land_use_mod)




