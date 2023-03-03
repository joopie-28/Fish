# Environmental analyses version 2.

# Exclude time series outside timeframe ? 

exclusion_vec <- NULL
for (i in 1:nrow(full.novel.mat.sf)){
  print(i)
  if(full.novel.mat.sf$bins[i] > 40 | full.novel.mat.sf$bins[i] < 7 ){
    exclusion_vec <- c(exclusion_vec, as.character(full.novel.mat.sf$site[i]))
  }
}

remove.indices <- which(full.novel.mat.sf$site %in% unique(exclusion_vec))

full.novel.mat.filtered <- full.novel.mat.sf[-remove.indices,]



# Load ncdf4 library to work with netcdf files
library(ncdf4) 

# Map nc file as a RasterBrick, where bands denote months.
surface.temp.brick <- brick("/Users/sassen/Desktop/waterTemperature_monthly_1981-2014.nc")

surface.flow.brick <- brick("/Users/sassen/Desktop/discharge_Global_monthly_1979-2014.nc")

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

create.sd.raster <- function(brick.data, n.years){
  
  # Create empty list
  output <- list()
  
  # Create month vector
  vec <- 1:12
  
  # Based on how many years the data spans, compute
  # averages using terra package.
  
  for (i in 1:n.years){
    
    print(i)
    output[[i]] <- calc(brick.data[[vec]], fun=sd)
    
    vec <- vec + 12
  }
  # Name layers such that they match our community rows.
  names(output) <- 2021 - ((2015-n.years):2014)
  
  return(output)
}

# Calculate the sd for each year

output_temp_sd <- create.sd.raster(surface.temp.brick, 34)

output_flow_sd <- create.sd.raster(surface.flow.brick, 36)

# Instantiate a new column
full.novel.mat.sf$surface_temp_sd <- 0

full.novel.mat.sf$surface_flow_sd <- 0

# Add to modelling frame

for (i in 1:nrow(full.novel.mat.sf)){
  print(i)
  # Temporal positioning, easier with full.novel.mat
  index <- as.character(full.novel.mat.sf$bins[i])
  
  # Assign the correct raster band to a temporary object
  current.ras.temp.sd <- output_temp_sd[[index]]
  # Control for non-existent data (we have only 1981-2014)
  if(is.null(current.ras.temp.sd)){
    full.novel.mat.sf$surface_temp_sd[i] <- NA
    
  }
  else{
    # Extract temperatures for community within a given year
    full.novel.mat.sf$surface_temp_sd[i] <- as.double(raster::extract(current.ras.temp.sd, full.novel.mat.sf[i,]))
  }
}



for (i in 1:nrow(full.novel.mat.sf)){
  print(i)
  # Temporal positioning, easier with full.novel.mat
  index <- as.character(full.novel.mat.sf$bins[i])
  
  # Assign the correct raster band to a temporary object
  current.ras.flow <- output_flow_sd[[index]]
  # Control for non-existent data (we have only 1981-2014)
  if(is.null(current.ras.flow)){
    full.novel.mat.sf$surface_flow_sd[i] <- NA
    
  }
  else{
    # Extract temperatures for community within a given year
    full.novel.mat.sf$surface_flow_sd[i] <- as.double(raster::extract(current.ras.flow, full.novel.mat.sf[i,]))
  }
}

full.novel.mat.sf.copy <- full.novel.mat.filtered
full.novel.mat.sf.copy$bin_lag<-scale(full.novel.mat.filtered$bin_lag, center = T, scale = T)
full.novel.mat.sf.copy$position<-scale(full.novel.mat.filtered$position, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp_sd<-scale(full.novel.mat.filtered$surface_temp_sd, center = T, scale = T)
full.novel.mat.sf.copy$surface_flow_sd<-scale(full.novel.mat.filtered$surface_flow_sd, center = T, scale = T)
full.novel.mat.sf.copy$INC_increase<-scale(full.novel.mat.filtered$INC_increase, center = T, scale = T)

full.novel.mat.sf.copy <- full.novel.mat.sf
full.novel.mat.sf.copy$bin_lag<-scale(full.novel.mat.sf$bin_lag, center = T, scale = T)
full.novel.mat.sf.copy$position<-scale(full.novel.mat.sf$position, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp_sd<-scale(full.novel.mat.sf$surface_temp_sd, center = T, scale = T)
full.novel.mat.sf.copy$surface_flow_sd<-scale(full.novel.mat.sf$surface_flow_sd, center = T, scale = T)
full.novel.mat.sf.copy$INC_increase<-scale(full.novel.mat.sf$INC_increase, center = T, scale = T)
full.novel.mat.sf.copy$surface_flow_delta<-scale(full.novel.mat.sf$surface_flow_delta, center = T, scale = T)
full.novel.mat.sf.copy$surface_temp_delta<-scale(full.novel.mat.sf$surface_temp_delta, center = T, scale = T)

full.mod<-lme4::glmer(novel~bin_lag+position+surface_temp_sd+surface_flow_sd+INC_increase+(1|country), family = "binomial", 
                data= full.novel.mat.sf.copy)


res.mod<-lme4::glmer(novel~bin_lag+position+surface_temp_sd+surface_flow_sd+INC_increase+(1|country), family = "binomial", 
                      data= full.novel.mat.sf.copy)

# Complete timeseries removed
summary(res.mod)

# Only missing values removed
summary(full.mod)



