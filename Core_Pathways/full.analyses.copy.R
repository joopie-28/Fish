#################################################################
# Novel Communities of Freshwater Fishes: a global investigation #
# into their emergence and persistence over the last 50 years	####
##################################################################

################ Reproducible Main Analyses ######################

# J.M. Sassen 
# 14-12-2022 

#### Set-up Environment ####

# source functions from 'functions' sub-folder
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

# Packages
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
                 "gridExtra", "clustsig", "dendextend"))



# Load all the exotic species databases (constructed manually from literature)
exotic.list <- sapply(list.files("./Exotics_databases_countries", pattern="\\.rds", full.names=TRUE), 
                      readRDS)

# Assign appropriate names for all the files (41...)
for (i in 1:length(exotic.list)){
  
  name <- names(exotic.list)[i]
  
  split.1 <- str_split(name, pattern = "/")[[1]][3]
  split.2 <- str_split(split.1, pattern = ".rds")[[1]][1]
  
  object <- exotic.list[[i]]
  
  assign(split.2, object)
  
  rm(object, name, split.1,split.2)
}
rm(exotic.list)

# Import data from RivFishTime
time_series_data <- read.csv("inputs/1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("inputs/1873_2_RivFishTIME_SurveyTable.csv")

# Import environmental data

# Population
pop.data <- raster("/Users/sassen/Desktop/Fish_GeoData/gpw-v4-population-density-rev11_2015_1_deg_tif/gpw_v4_population_density_rev11_2015_1_deg.tif")

# Human modification of the environment
anthromod.data <- raster("/Users/sassen/Desktop/Fish_GeoData/lulc-human-modification-terrestrial-systems-geographic-geotiff/lulc-human-modification-terrestrial-systems_geographic.tif")

# Surface temperature data 1981-2014
surface.temp.brick <- brick("/Users/sassen/Desktop/waterTemperature_monthly_1981-2014.nc")

# Water flow data 1979-2014
surface.flow.brick <- brick("/Users/sassen/Desktop/discharge_Global_monthly_1979-2014.nc")



# Load CRS
WG84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

###### PHASE 1 - PRE-PROCESSING AND MANIPULATION ######

#### Pre-processing 1. Construct community matrices and detect novelty ####

# These are the countries with sufficient data for bin range = 1 year, skipping countries with
# no usable time series saves us some time.

countries.suf.data <- list("BEL", "BRA", "BWA", "CAN", "CIV", "COL", 
                           "ESP", "FIN", "FRA", "GBR", "HUN", "JPN", 
                           "SWE", "USA", "AUS")

# Create a grand list of Time Series which can be used for analysis

full.ID.list <- do.call(c, lapply(countries.suf.data, function(country){
  test <- country_list_assigner(country)
  return(test)
}))

ID.list <- as.list(time_series_data$TimeSeriesID)

matrix_list_seasonality <- list_matrix_seasonality_function(ID.list)

nov.output<-novelty.detection.gam(matrix_list_seasonality)

#### Pre-processing 2. Tag species as 'Invader', 'Non-native', or 'Native'####

full.stat.matrices.season <- assign.stat.country.V2(names(matrix_list_seasonality))


#### Pre-processing 3. Compute turnover metrics for each group of species ####

full.nnc.matrices.season <- mat.nnc.ass.V2(full.stat.matrices.season)


#### Pre-processing 4. Create a data frame that can be used for modelling ####

full.novel.mat.season <- inv.frame.builder.V2(full.nnc.matrices.season)

# Partition site and season names for efficient cross-referencing later on
temp_df <- data.frame(do.call("rbind", strsplit(as.character(full.novel.mat.season$site), ".",
                                                fixed = TRUE)))
names(temp_df) <- c("site_ID", "Quarter")
full.novel.mat.season<-cbind(temp_df, full.novel.mat.season)
full.novel.mat.season$site <- as.character(full.novel.mat.season$site)


#### Pre-processing 5. Classifying novel communities using Similarity Profiles ####

# Filter matrices where novelty occurred
nov.matrices <- matrix_list_seasonality[full.novel.mat.season$site[which(full.novel.mat.season$cat == 'novel')]]

# Run persistence framework. Here we apply hierarchical clustering 
# and use type I SIMPROF to find which clusters are significantly
# different at alpha = 0.05.

novelty.pers <- do.call(c, 
                        lapply(1:length(nov.matrices), function(x){
                          print(x)
                          nov.cluster.id.V6(nov.matrices[x])
                          
                        })) 

# Write the results back into our main data frame.
full.novel.mat.season$novel.length <- NA
full.novel.mat.season$novel.class <- NA

for(i in 1:length(novelty.pers)){
  indices <- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel")
  if(length(indices) > 1){
    indices<- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel" & full.novel.mat.season$bins == novelty.pers[[i]][[1]]$begin)
  }
  
  full.novel.mat.season$novel.length[indices] <- novelty.pers[[i]][[1]]$length 
  full.novel.mat.season$novel.class[indices] <-  novelty.pers[[i]][[1]]$Class
  
}

# Set non-novel lengths to 0 instead of NA
full.novel.mat.season$novel.length[is.na(full.novel.mat.season$novel.length)] <- 0
full.novel.mat.season$novel.class[is.na(full.novel.mat.season$novel.class)] <- "NONE"


#### Pre-processing 6. Computing demographic processes and how they relate to novelty ####

# Add transition groups to the main data frame
full.novel.mat.season$novel.groups <- "Non-Novel"
full.novel.mat.season$prev.group <- "Non-Novel"

for (i in 1:length(novelty.pers)){
  site <- names(novelty.pers[[i]])
  print(i)
  print(site)
  matrix<-full.novel.mat.season[which(full.novel.mat.season$site == site),]
  nov.data <- novelty.pers[[i]][[site]]
  matrix.expanded<-classify.compositions(matrix, pers.data=nov.data)
  
  full.novel.mat.season[which(full.novel.mat.season$site == site),] <- matrix.expanded
}

full.novel.mat.season$transition <- paste0(full.novel.mat.season$prev.group, ":", full.novel.mat.season$novel.groups)

# Create a new data frame with all turnover parameters
demography.frame <- rbindlist(lapply(1:length(matrix_list_seasonality), function(x){
  print(x)
  
  demog.df<-local.extorig.RivFishTime(ssmat=matrix_list_seasonality[[x]],novel=nov.output[[x]],
                                      site.name = names(nov.output[x]))
  demo.inv.df <- invaders.extorig.RivFishTime(ssmat=full.stat.matrices.season[[x]],novel=nov.output[[x]],
                                              site.name = names(nov.output[x]))
  demo.inv.df<-demo.inv.df[,-c(1:7)]
  colnames(demo.inv.df) <- paste0(colnames(demo.inv.df), ".inv" )
  demo.df<-cbind(demog.df,demo.inv.df)
  
  return(demo.df)
}))

# Add data to our main modelling frame
full.novel.mat.season<-cbind(full.novel.mat.season, demography.frame[,-c(1,2,3,5,6)])
full.novel.mat.season$basin <- as.factor(full.novel.mat.season$basin)


#### Pre-processing 7. Importing and extracting environmental drivers at time series level #####

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

# Combine novelty data with environmental data in a spatial data frame.
geo.timeseries.full <- cbind(geo.timeseries.sf, 
                             rbindlist(lapply(geo.timeseries.sf$ID, 
                                              function(site){
                                                print(site)
                                                # Add up all the novelty metrics for binomial model
                                                back <- length(which(full.novel.mat.season$site == site & 
                                                                       full.novel.mat.season$cat == "back"))
                                                
                                                instant <- length(which(full.novel.mat.season$site == site & 
                                                                          full.novel.mat.season$cat == "instant"))
                                                
                                                cumul <- length(which(full.novel.mat.season$site == site & 
                                                                        full.novel.mat.season$cat == "cumul"))
                                                
                                                novel <- length(which(full.novel.mat.season$site == site & 
                                                                        full.novel.mat.season$cat == "novel"))
                                                
                                                # Include novelty classes based on persistence length for
                                                # model variation.
                                                if(novel > 0){
                                                  index <- (which(full.novel.mat.season$site == site & 
                                                                    full.novel.mat.season$cat == "novel"))[1]
                                                  print(index)
                                                  
                                                  novelty.class <- full.novel.mat.season[[index, "novel.class"]]
                                                }
                                                else{
                                                  novelty.class <- "NONE"
                                                }
                                                
                                                # Add Mean invader increase in each time series
                                                indices <- which(full.novel.mat.season$site== site)
                                                
                                                
                                                Country <- full.novel.mat.season$country[indices][1]
                                                BioRealm <- full.novel.mat.season$BioRealm[indices][1]
                                                Basin <- full.novel.mat.season$basin[indices][1]
                                                
                                                # Add in total length of timeseries as a covariate
                                                length <- back + instant + cumul +novel
                                                
                                                # Return a clean df for modelling
                                                df <- data.frame("back" =back, "instant" = instant, 
                                                                 "cumul" = cumul, "novel"= novel, "length" = length, "Class" = novelty.class,
                                                                 "Country" = Country, 'BioRealm' = BioRealm, 
                                                                 'Basin' = Basin)
                                                return(df)
                                              })))

# Extract aanthropogenic modification data
geo.timeseries.full$Anthromod <- raster::extract(anthromod.data, geo.timeseries.full)

# Extract population
geo.timeseries.full$pop.2015 <- raster::extract(pop.data, geo.timeseries.full)


# We're interested in persistent novelty rather than blips
geo.timeseries.full$Class_Cleaned <- geo.timeseries.full$Class

for(i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$Class_Cleaned[i] == "BLIP"){
    geo.timeseries.full$Class_Cleaned[i] <- "NONE"
  }
}

# Add a binary variable denoting whether novelty occurs within a time series at all
for (i in 1:nrow(geo.timeseries.full)){
  if(geo.timeseries.full$novel[i] > 0 ){
    geo.timeseries.full$binary_novel[i] <- 1
  }
  else{
    geo.timeseries.full$binary_novel[i] <- 0
  }
  
}

# Focus will be population and anthropogenic impact. May include builtup?


#### Pre-processing 8. Importing and extracting environmental drivers at community level ####

create.avg.raster <- function(brick.data, n.years){
  
  n.quarter<-n.years*4
  # Create empty list
  output <- list()
  
  # Create quarter vector
  vec <- 1:3
  
  # Based on how many years the data spans, compute
  # averages using terra package.
  
  for (i in 1:n.quarter){
    
    print(i)
    output[[i]] <- terra::mean(brick.data[[vec]])
    
    vec <- vec + 3
    
  }
  # Name layers such that they match our community rows.
  names_vec <- NULL
  for(i in 2014:1981){
    temp<-rep(i,4)
    temp[1]<-paste0(temp[1], ".4")
    temp[2]<-paste0(temp[2], '.3')
    temp[3]<-paste0(temp[3], '.2')
    temp[4]<-paste0(temp[4], '.1')
    names_vec <- c(names_vec, temp)
  }
  names(output) <- names_vec
  
  return(output)
}

temp_data_comm<-create.avg.raster(surface.temp.brick, 34)
flow_data_comm<-create.avg.raster(surface.flow.brick, 36)
temp_test<-flow_data_comm
saveRDS(temp_data_comm, "./outputs/SURFACE_TEMP.rds")
saveRDS(flow_data_comm, "./outputs/SURFACE_TEMP.rds")
# Which time series have data beyond 1981-2014?
# add lat and long to this

env.comm.level <- rbindlist(lapply(names(nov.output), function(ID){
  
  mat <- identify.novel.gam(site.sp.mat = matrix_list_seasonality[[ID]], 
                            alpha = 0.05,
                            metric = "bray",
                            plot =F, 
                            site = ID,
                            plot.data = FALSE,
                            gam.max.k = -1)
  print(ID)
  mat$surface_temp <- NA
  mat$suface_temp_prev_season <- NA
  #mat$surface_temp_delta <- NA
  # mat$surface_temp_lag1 <- NA
  mat$Longitude <- subset(time_series_data, TimeSeriesID == strsplit(ID, split="[.]")[[1]][1])$Longitude
  mat$Latitude <- subset(time_series_data, TimeSeriesID == strsplit(ID, split="[.]")[[1]][1])$Latitude
  mat$timeseries <- strsplit(ID, split="[.]")[[1]][1]
  
  mat <- st_as_sf(mat, coords = c('Longitude', "Latitude"), crs=WG84)
  for (i in 1:nrow(mat)){
    
    # Temporal positioning, easier with full.novel.mat
    index.bin <- as.character(2021-as.numeric(mat$bins[i]))
    index.quarter <- strsplit(ID, split="[.]")[[1]][2]
    
    if(index.quarter == "4/1"){
      index.quarter <- "4"
    }
    lag.quarter <- as.character(as.numeric(index.quarter) - 1)
    if(lag.quarter == "0"){
      lag.quarter <- "4"
      index.bin <- as.character(as.numeric(index.bin)-1)
    }
    
    index<-paste0(index.bin, ".", index.quarter)
    lag_index <- paste0(index.bin, ".", lag.quarter)
    # Assign the correct raster band to a temporary object
    current.ras.temp <- temp_test[[index]]
    current.lag.temp <- temp_test[[lag_index]]
    # Control for non-existent data (we have only 1981-2014)
    if(is.null(current.ras.temp)){
      mat$surface_temp[i] <- NA
      
    }
    else{
      # Extract temperatures for community within a given year
      mat$surface_temp[i] <- as.double(raster::extract(current.ras.temp, mat[i,]))
      
    }
    # Control for non-existent data (we have only 1981-2014)
    if(is.null(current.lag.temp)){
      mat$surface_temp_prev_season[i] <- NA
      
    }
    else{
      # Extract temperatures for community within a given year
      mat$surface_temp_prev_season[i] <- as.double(raster::extract(current.lag.temp, mat[i,]))
      
    }
    
  }
  
  #for(i in 2:nrow(mat)){
  # mat$surface_temp_delta[i] <- mat$surface_temp[i] - mat$surface_temp[i-1]
  #  mat$surface_temp_delta_lag1[i] <- mat$surface_temp_delta[i-1]
  # mat$surface_temp_lag1[i] <- mat$surface_temp[i-1]
  # }
  mat <- mat[-c(1:5),c("site","bins", "surface_temp","surface_temp_prev_season", "cat", "novel", 'cumul', 'instant', 'timeseries', 'bin.lag')]
  return(mat)
}))

# Remove timeseries with less than x points
env.test<- env.comm.level
for(ID in unique(env.comm.level$site)){
  print(ID)
  if(nrow(subset(env.comm.level, site == ID & !is.na(surface_temp))) < 5){
    env.test[which(env.test$site == ID),] <- NA
  }
}







# Need to add lat and long

full.env.sf <- st_as_sf( filter.1, coords = c('Longitude', "Latitude"), crs=WG84)

# Instantiate a new column
full.env.sf$surface_temp <- 0
full.env.sf$surface_temp_delta <- 0
# Add to modelling frame

for (i in 1:nrow(full.env.sf)){
  print(i)
  # Temporal positioning, easier with full.novel.mat
  index.bin <- as.character(2021-full.env.sf$bins[i])
  index.quarter <- as.character(full.env.sf$Quarter[i])
  index<-paste0(index.bin, ".", index.quarter)
  # Assign the correct raster band to a temporary object
  current.ras.temp <- temp_test[[index]]
  # Control for non-existent data (we have only 1981-2014)
  if(is.null(current.ras.temp)){
    full.env.sf$surface_temp[i] <- NA
    
  }
  else{
    # Extract temperatures for community within a given year
    full.env.sf$surface_temp[i] <- as.double(raster::extract(current.ras.temp, full.env.sf[i,]))
  }
}

for (i in unique(full.env.sf$site)){
  print(i)
  
  # Isolate the time series within the full frame
  temp <- subset(full.env.sf, site == i)
  temp$delta.temp <- NA
  
  
  # Compute delta temperature where applicable
  for (j in 2:nrow(temp)){
    
    temp$delta.temp[j] <- abs(temp$surface_temp[j] - temp$surface_temp[j-1])
  }
  
  # Feed this subset of delta's back into the full data frame.
  indices <- which(full.env.sf$site == temp$site)
  
  for (k in 1:nrow(temp)){
    full.env.sf$surface_temp_delta[indices[k]] <- temp$delta.temp[k]
    
  }
}

# filter the timeseries that have NA's
# what about increase lag?
mod<-(glmer(novel+0~surface_temp +(1|site) , data=env.test, family='binomial'))
summary(mod)
#### Pre-processing 9. Computing lag between novel community and succeeding state ####

# Add a new feature representing the lag between the novel community T and the 
# succeeding community T + 1. This could  affect the size of the novel cluster.
# We thus filter out time series where the lag between a novel community and 
# the succeeding community seriously affects our ability to infer anything about
# the persistence of that community.

# Compute the lags for each community (slightly different to earlier defined bin lag variable)

full.nov.no.lag <- filter.by.lag(full.novel.mat.season, 2)

full.nov.no.lag <- full.nov.no.lag[-which(full.nov.no.lag$lag_to_next == 1000),]

full.nov.no.emergence.lag <- filter.by.lag.novel(full.novel.mat.season, 2)

full.nov.no.emergence.lag <- full.nov.no.lag[-which(full.nov.no.emergence.lag$lag_to_next == 1000),]

# Remove time series where the novel community emerges with signfiicant lag.



###### PHASE 2 - MODELLING AND ANALYSES ######

#### Modelling 1. Rates of Novelty Emergence around the globe ####

# Three separate models for each novelty type, at the community level
fixed.emergence.nov.mod <- glmer(novel ~ bin_lag +position+ (1|Quarter/site_ID), data = full.novel.mat.season, family= 'binomial')

fixed.emergence.cumul.mod <- glmer(cumul ~  bin_lag +position + (1|Quarter/site_ID), data = full.novel.mat.season, family= 'binomial')

fixed.emergence.instant.mod <- glmer(instant ~bin_lag +position + (1|Quarter/site_ID), data = full.novel.mat.season, family= 'binomial')

emergence.mod.list <- list(fixed.emergence.instant.mod, fixed.emergence.cumul.mod, fixed.emergence.nov.mod)

# One model at the time series level
broad.emergence.mod <- glmer(binary_novel ~ (1|Country), data = geo.timeseries.full, family = 'binomial')


#### Modelling 2. Effects of invader presence on probability of emergence ####

# Demographic turnover of invaders during transitions

mod.NAC <-glm(novel~NAC_increase+bin_lag+position, data = full.novel.mat.season, family ="binomial")
mod.INV <-glmer(novel~INC_increase+bin_lag+position + (1|basin), data = full.novel.mat.season, family ="binomial")

# Species no more likely, but relative abundance more likely yo trigger novelty framework than native

#### Modelling 3. Demographic Processes: Turnover during transitions #####

# Extinction and Origination Models
ext.mod <- glm(cbind(ext, (ext.rich-ext))~transition, data=full.novel.mat.season, family='binomial')
orig.mod <- glm(cbind(orig, (orig.rich-orig))~transition, data=full.novel.mat.season, family='binomial')

# Emigration and Immigration models
emig.mod <- glm(cbind(emig, (ext.rich-emig))~transition , data=full.novel.mat.season, family='binomial')
immig.mod <- glm(cbind(immig, (orig.rich-immig))~transition , data=full.novel.mat.season, family='binomial')

# Full taxa loss and gain models
loss.mod <- glm(cbind(loss, (ext.rich-loss))~transition , data=full.novel.mat.season, family='binomial')
gain.mod <- glm(cbind(gain, (orig.rich-gain))~transition , data=full.novel.mat.season, family='binomial')



#### Modelling 4. Correlation of human activity with Novelty ####

environ.driver.mod <- glm(binary_novel~Anthromod+length, data= geo.timeseries.full, family = 'binomial')



###### PHASE 3 - PRODUCING TABLES AND FIGURES ######

#### Figure 1. Map of Novelty for Austrlalasia, Palearctic and Nearctic, with rates. ####
pdf(file = "/Users/sassen/Desktop/Figure_1a-d.pdf",
    width = 15,
    height = 12)

map.realm.plotter(full.novel.mat.season)

dev.off()

#### Figure 2. Visualising model predictions for invader presence ####
pdf(file = "/Users/sassen/Desktop/Figure_4.pdf",
    width = 12,
    height = 10)
visreg(mod.INV, 'INC_increase' ,scale = 'response', rug=F, ylim = c(0,1),ylab="Probability of Novel Community Emergence",
       xlab="Change in relative abundance of exotics")
points(novel~INC_increase, data= full.novel.mat.season,pch=19, col = alpha('grey', alpha=0.4), cex=0.8)


dev.off()
#### Figure 3. Scatter plot showing transitions from novelty to the next community on an axis of demographic processes ####

# Plot 3a. Immigration versus Origination

pdf(file = "/Users/sassen/Desktop/Figure_3a-b.pdf",
    width = 9,
    height = 15)

figure.3.turnover(emig.mod,ext.mod, immig.mod, orig.mod)

dev.off()


#### Figure 4. Visualising model predictions for population/anthropogenic modification ####
pdf(file = "/Users/sassen/Desktop/Figure_5.pdf",
    width = 12,
    height = 10)

visreg(environ.driver.mod, 'Anthromod' ,scale = 'response', rug=F, ylim = c(0,.2), ylab = "Probability of persistent novelty emergence",
       xlab = "Anthropogenic Modification Index" ,line=(list(col = "blue")))
points(binary_novel~Anthromod, data= geo.timeseries.full,pch=19, col = alpha('grey', alpha=0.4), cex=0.8)
points(y=rep(0.2, nrow(subset(geo.timeseries.full,binary_novel == 1))), x= subset(geo.timeseries.full,binary_novel == 1)$Anthromod, pch=19, col = alpha('grey', alpha=0.4), cex=0.8)
dev.off()

###### END OF MAIN ANALYSES ######


#### Table 1. 

#### Table 2.

#### Table 3.

#### Table 4.
# To do #
# Complete this document
# Formalize the environmental drivers - make a decision
# Time to flesh out this origination stuff - make a decision wether or not it it worth having.
# Write up results and the rest of the paper. That's it!



table(full.novel.mat.season$novel)

test<- aggregate(full.novel.mat.season$novel, 
                 by = list(full.novel.mat.season$site_ID), 
                 FUN = sum)
test$binary <- NA
nrow(subset(test,x >0))/nrow(test)
for(i in 1:nrow(test)){
  if(test$x[i] > 0){
    test$binary[i] <- 1
  }else{
    test$binary[i] <- 0
  }
}

mod<-glm(binary~1, data=test, family = 'binomial')
plogis(confint(mod))

