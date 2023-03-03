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
sites <- as.character(unique(full.novel.mat.season$site))

# Cross-reference with survey data to obtain coordinates
index <- which(time_series_data$TimeSeriesID %in% sites)

# Find relevant data and convert to SF
geo.timeseries <- time_series_data[index, c("TimeSeriesID", "Longitude", "Latitude")]

geo.timeseries.sf <- st_as_sf(geo.timeseries, coords = c("Longitude", "Latitude"), crs = WG84) 


#### Step 1. Extracting geographic data from multiple sources ####


# Elevation, from the WORLDCLIM database
elevation.data <- raster("/Users/sassen/Desktop/Fish_GeoData/wc2.1_2.5m_elev.tif")

# Land-use Globcover 2009
land.use.data <- raster("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER_L4_200901_200912_V2.3.tif")

# Also import codes legend for land use
land.use.legend <- read.csv("/Users/sassen/Desktop/Fish_GeoData/GLOBCOVER2009_legend.csv")

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

# Calculate percentage of land use as per Comte et al.

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
  
# Buffer size in meters, controls the area around the point which we collect pixels for.
  
geo.timeseries.sf$Land_use_impact<-anthro.impact.function(geo.timeseries.sf, 
                                                          buffer_size = 25000)


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
  
  mean.INV.inc <- sum(full.novel.mat[indices, "INC_increase"])/length(indices)
  
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

# This is just a placeholder function which contains everything.

enviro.models <- function(geo.timeseries.full){


  geo.timeseries.full.copy <- geo.timeseries.full
  
  # Scale all the variables
  geo.timeseries.full.copy$Mean_Invader_Increase <- scale(geo.timeseries.full.copy$Mean_Invader_Increase, center = T, scale = T)
  geo.timeseries.full.copy$mean.delta <- scale(geo.timeseries.full.copy$mean.delta, center = T, scale = T)
  geo.timeseries.full.copy$tmin.mean.1969 <- scale(geo.timeseries.full.copy$tmin.mean.1969, center = T, scale = T)
  geo.timeseries.full.copy$Anthromod <- scale(geo.timeseries.full.copy$Anthromod, center = T, scale = T)
  geo.timeseries.full.copy$Elevation <- scale(geo.timeseries.full.copy$Elevation, center = T, scale = T)
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


  summary(fit)

  # Exclude bips
  fit.no.blip <- glm(cbind(novel, (length-novel))~ Land_use_impact+Mean_Invader_Increase, 
                     data = subset(geo.timeseries.full.copy , Class !="BLIP" & Class != "END"), family = "binomial")

  fit.no.blip <- glm(binary~length+ Land_use_impact+Mean_Invader_Increase , 
                   data = subset(geo.timeseries.full.copy, Class !="BLIP" & Class != "END"), family = "binomial")
  

  summary(fit.no.blip)

  # Add coefficients to output to make inspection easier
  mod.1 <- list("model" = fit.no.blip, "summary" = as.data.frame(summary(fit.no.blip)$coef))
  
  
  # Write a csv file for later use.
  
  write.csv(as.data.frame(summary(fit.no.blip)$coef), paste0("./outputs/enviro_glm.csv"))

  return(mod.1)
}

enviro.models(geo.timeseries.full)

# Model per timeseries (turn into df)

timeseries.level.rates <- glm(binary~1, 
           data = (geo.timeseries.full.copy), family = "binomial")

test.data<-subset(geo.timeseries.full, Class == "Persister" | Class=="NONE" & Country == 'FRA')

mod<-(glm(cbind(novel, (length))~Pop.2015+Anthromod+Builtup.2015, 
          data = test.data, family = "binomial"))
summary(mod)
visreg(mod, c('mean.delta'), scale='response', partial=F, rug=F, 
       bty='l', xlab = 'Local Originations during Emergence', ylab = 'Persistence Probability of Novel Community',
       points = list(col='black'), line = list(col = 'green'), ylim = c(0,.15), axes=T)

points(novel~Anthromod, data = geo.timeseries.full, pch=19, col = alpha('black', alpha=0.2), cex=0.4)
vioplot(log(Pop.2015+1)~novel,data = geo.timeseries.full )


#### Step 3. Plot Model Results ####
 
plot_model(fitter, type = "pred",  terms=c("Anthromod [All]"), 
           title = "Probability of Novelty Emergence explained x",
           axis.title = c("Net change in x","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green", show.data = F)



#### End of Analyses ####










## Deprecated ##

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


# Use output to assign new values

geo.timeseries.sf$Land_use_mod <- NA
adj.land.use<-adj.cell.finder(geo.timeseries.sf, land.use.data, directions = 8)

# Some categories are more interesting than others, i.e. we are not that concerned
# with the difference between an open forest or closed forest for our questions.
# Therefore, we will merge some groups which will allow for a more meaningful comparison.
# Here, I create a new df which holds info on original classes and new classes.

land.use.legend$New_Class <- c(rep("Cropland", 4), rep("Natural Vegetation", 14), "Urban and Artificial", "Bare area", "Water body", "Snow and Ice", "No data")







########### Surface Stream flow and Temperature #######
######## At community level ##############

# Load ncdf4 library to work with netcdf files
library(ncdf4) 


# Map nc file as a RasterBrick, where bands denote months.
surface.temp.brick <- brick("/Users/sassen/Desktop/waterTemperature_monthly_1981-2014.nc")

surface.flow.brick <- brick("/Users/sassen/Desktop/discharge_Global_monthly_1979-2014.nc")

# Add spatial data to modelling frame
for(i in 1:nrow(full.novel.mat)){
  print(i)
  # Simply cross-reference sites and add in latitude and longitude
  index<-which(geo.timeseries$TimeSeriesID==full.novel.mat$site[i])
  full.novel.mat$lat[i] <- geo.timeseries$Latitude[index]
  full.novel.mat$long[i] <- geo.timeseries$Longitude[index]
}

# Add a CRS
full.novel.mat.sf<- st_as_sf(full.novel.mat, coords = c("long", "lat"), crs = WG84) 

# Step 1. Create new object with yearly average rasters.

create.avg.raster <- function(brick.data, n.years){
  
  # Create empty list
  output <- list()
  
  # Create month vector
  vec <- 1:12
  
  # Based on how many years the data spans, compute
  # averages using terra package.

  for (i in 1:n.years){

    print(i)
    output[[i]] <- terra::mean(brick.data[[vec]])
      
    vec <- vec + 12
  }
  # Name layers such that they match our community rows.
  names(output) <- 2021 - ((2015-n.years):2014)
  
  return(output)
}

output <- create.avg.raster(surface.temp.brick, 34)

output_flow <- create.avg.raster(surface.flow.brick, 36)

# Step 3. Extract surface temperatures from raster files

# Instantiate a new column
full.novel.mat.sf$surface_temp <- 0
full.novel.mat.sf$surface_flow <- 0

for (i in 1:nrow(full.novel.mat.sf)){
  print(i)
  # Temporal positioning, easier with full.novel.mat
  index <- as.character(full.novel.mat.sf$bins[i])
  
  # Assign the correct raster band to a temporary object
  current.ras.temp <- output[[index]]
  # Control for non-existent data (we have only 1981-2014)
  if(is.null(current.ras.temp)){
    full.novel.mat.sf$surface_temp[i] <- NA
    
  }
  else{
    # Extract temperatures for community within a given year
    full.novel.mat.sf$surface_temp[i] <- as.double(raster::extract(current.ras.temp, full.novel.mat.sf[i,]))
  }
}

for (i in 1:nrow(full.novel.mat.sf)){
  print(i)
  # Temporal positioning, easier with full.novel.mat
  index <- as.character(full.novel.mat.sf$bins[i])
  
  # Assign the correct raster band to a temporary object
  current.ras.flow <- output_flow[[index]]
  # Control for non-existent data (we have only 1981-2014)
  if(is.null(current.ras.flow)){
    full.novel.mat.sf$surface_flow[i] <- NA
    
  }
  else{
    # Extract temperatures for community within a given year
    full.novel.mat.sf$surface_flow[i] <- as.double(raster::extract(current.ras.flow, full.novel.mat.sf[i,]))
  }
}



# Step 4. Create a DELTA temperature variable for modelling

full.novel.mat.sf$surface_temp_delta <- NA
full.novel.mat.sf$surface_flow_delta <- NA

for (i in unique(full.novel.mat.sf$site)){
  print(i)
  
  # Isolate the time series within the full frame
  temp <- subset(full.novel.mat.sf, site == i)
  temp$delta.temp <- NA
  temp$delta.temp.sq <- NA
  temp$delta.temp.rel <- NA
  
  # Compute delta temperature where applicable
  for (j in 2:nrow(temp)){
    
    temp$delta.temp[j] <- abs(temp$surface_temp[j] - temp$surface_temp[j-1])
    temp$delta.temp.sq[j] <- (temp$surface_temp[j] - temp$surface_temp[j-1])^2
    temp$delta.temp.rel[j] <- (abs(temp$surface_temp[j] - temp$surface_temp[j-1]))/max(temp$surface_temp, na.rm = T)
  }
  
  # Feed this subset of delta's back into the full data frame.
  indices <- which(full.novel.mat.sf$site == temp$site)
  
  for (k in 1:nrow(temp)){
    full.novel.mat.sf$surface_temp_delta[indices[k]] <- temp$delta.temp[k]
    full.novel.mat.sf$surface_temp_delta.sq[indices[k]] <- temp$delta.temp.sq[k]
    full.novel.mat.sf$surface_temp_delta.rel[indices[k]] <- temp$delta.temp.rel[k]
  }
}

# Flow
for (i in unique(full.novel.mat.sf$site)){
  print(i)
  
  # Isolate the time series within the full frame
  temp <- subset(full.novel.mat.sf, site == i)
  temp$delta.flow <- NA
  temp$delta.flow.sq <- NA
  temp$delta.flow.rel <- NA
  
  # Compute delta temperature where applicable
  for (j in 2:nrow(temp)){
    
    temp$delta.flow[j] <- abs(temp$surface_flow[j] - temp$surface_flow[j-1])
    temp$delta.flow.sq[j] <- (temp$surface_flow[j] - temp$surface_flow[j-1])^2
    temp$delta.flow.rel[j] <- (abs(temp$surface_flow[j] - temp$surface_flow[j-1]))/max(temp$surface_flow, na.rm = T)
  }
  
  # Feed this subset of delta's back into the full data frame.
  indices <- which(full.novel.mat.sf$site == temp$site)
  
  for (k in 1:nrow(temp)){
    full.novel.mat.sf$surface_flow_delta[indices[k]] <- temp$delta.flow[k]
    full.novel.mat.sf$surface_flow_delta.sq[indices[k]] <- temp$delta.flow.sq[k]
    full.novel.mat.sf$surface_flow_delta.rel[indices[k]] <- temp$delta.flow.rel[k]
  }
}

# Step 5. Run a model
full.novel.mat.sf.copy <- full.novel.mat.sf

full.novel.mat.sf.copy$surface_flow_delta.sq<-scale(full.novel.mat.sf$surface_flow_delta.sq, center = T, scale = T)
full.novel.mat.sf.copy$surface_flow_delta<-scale(full.novel.mat.sf$surface_flow_delta, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp_delta.sq<-scale(full.novel.mat.sf$surface_temp_delta.sq, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp_delta<-scale(full.novel.mat.sf$surface_temp_delta, center = T, scale = T)

full.novel.mat.sf.copy$bin_lag<-scale(full.novel.mat.sf$bin_lag, center = T, scale = T)
full.novel.mat.sf.copy$position<-scale(full.novel.mat.sf$position, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp<-scale(full.novel.mat.sf$surface_temp, center = T, scale = T)
full.novel.mat.sf.copy$surface_flow<-scale(full.novel.mat.sf$surface_flow, center = T, scale = T)
full.novel.mat.sf.copy$INC_increase<-scale(full.novel.mat.sf$INC_increase, center = T, scale = T)

temp.mod<-glmer(novel~bin_lag+position+surface_temp_delta*surface_flow_delta+(1|site), family = "binomial", 
            data= full.novel.mat.sf.copy)

summary(temp.mod)


plot(novel~Land_use_impact, data= geo.timeseries.full.copy)



plot_model(temp.mod, type = "pred",  terms=c("INC_increase [all]"), 
           title = "Probability of Novelty Emergence explained x",
           axis.title = c("Net change in x","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green", show.data = T)


# Incorporate seasonality

create.sbc.mat <- function(Survey_Data){
  
  # Isolate a single time series
  test<-subset(Survey_Data, TimeSeriesID == n.ID )
  # Create a species by time matrix 
  test.2<-as.data.frame(tapply(test$Abundance,
         list(test$Year, test$Species),
         sum, na.rm = T))
  # Convert absences (NA) to 0's
  test.2[is.na(test.2)]<-0

  # Convert to relative abundance
  test.3<- test.2/rowSums(test.2)
  
  # Convert rownames to time in past
  rownames(test.3) <- 2021-as.numeric(rownames(test.2))
  
return(test.3)

}

mat.1<-create.sbc.mat(Survey_Data)



# Isolate a single time series

create.seasonal.ssmat <- function(n.ID){
  df<-subset(Survey_Data, TimeSeriesID == n.ID )
  output<-lapply(unique(df$Quarter), function(quarter){
  print(quarter)
  test <- subset(df, Quarter == quarter)

  test.2<-as.data.frame(tapply(test$Abundance,
                               list(test$Year, test$Species),
                               sum, na.rm = T))
  # Convert absences (NA) to 0's
  test.2[is.na(test.2)]<-0
  
  # Convert to relative abundance
  test.3<- test.2/rowSums(test.2)
  
  # Convert rownames to time in past
  rownames(test.3) <- 2021-as.numeric(rownames(test.2))
  if(nrow(test.3) < 10){
    test.3 <- NULL
  }
  return(test.3)
})
  names(output) <- paste0(n.ID,".", unique(df$Quarter))
  return(output)
}

view<-create.seasonal.ssmat('G1043')



view<-lapply(c('G1043', 'G1044'), function(x){
  create.seasonal.ssmat(x)
})

check<-unlist(view, recursive =F)








