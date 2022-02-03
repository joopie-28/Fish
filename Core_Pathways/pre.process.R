##### This is the script that includes the pre-processing steps #####

#### Characterizing the drivers of novel freshwater fish communities ####

###########################################################################
### This script includes main analyses as well as data import & tidying ###
### it requires functions to be loaded in from the function file, as well #
### as access to certain databases. #######################################
###########################################################################

# J.M. Sassen 31-01-2022 

##################################################
#### Step 0 . Load in packages and databases #####
##################################################

# Clear environment (if wanted) and set your working directory
rm(list = ls())

# source functions from 'functions' sub-folder and load them
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

package.loader(c("rgdal", "raster", "rgeos", "sf", 
                 "tidyverse", "rfishbase","mgcv",
                 "vegan", "lme4", "nlme", 
                 "DHARMa", "merTools", "shape",
                 "multcomp", "maptools", "sp", 
                 "divDyn", "plotrix", "raster",
                 "rgeos", "fun", "analogue",
                 "brms", "data.table", "data.table", "stringi"))

# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv") # RIVFishTime 
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv") # RIVFishTime

invasives_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Occurrence_Table.csv") # Tedesco et Al. Invasives database
colnames(invasives_data) <- c("Basin", "Species", "Status", "TSN_ITIS_Code", "FishBase_code", "Valid_name", "Sighting") # Tedesco et Al. Invasives database

basin_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Drainage_Basins_Table.csv") # Tedesco et Al. Basins geodatabase
colnames(basin_data) <- c("Basin", "Country", "BioRealm", "Endorheic", "Lon", "Lat", "Med_Lon", "Med_Lat", "Surface_Area") # Tedesco et Al. Basins geodatabase

occurence_shapefile <- read_sf(dsn = "/Users/sassen/Desktop/datatoFigshare/Basin042017_3119.shp") # Tedesco et AL. Shapefile

france_country_level <- france_species_status # Country level species data (Fishbase)


# Tidy data where necessary

invasives_data$Species <- gsub(invasives_data$Species, 
                               pattern ="\\.", 
                               replacement = " ")

invasives_data$Valid_name <- gsub(invasives_data$Valid_name, 
                                  pattern ="\\.", 
                                  replacement = " ")




##################################################
#### Step 1. Tagging species status  #############
##################################################

# Here, we will use a combination of the Tedesco et Al.
# database and our own custom species tags. 

# First we create a list of timeseries for all the countries we're 
# interested in.

fra.ID.list <- as.list(subset(time_series_data, Country == "FRA")$TimeSeriesID) # France
gbr.ID.list <- as.list(subset(time_series_data, Country == "GBR")$TimeSeriesID) # Great Britain
swe.ID.list <- as.list(subset(time_series_data, Country == "SWE")$TimeSeriesID) # Sweden
usa.ID.list <- as.list(subset(time_series_data, Country == "USA")$TimeSeriesID) # United States

# This function generates a dataframe matching a HydroBasin code to
# the name used in Tedesco et Al. This allows for tagging species
# on the basin level.

basin_name_code <- basin_name_match_function(occurence_shapefile)

# Use these data to add a basin_name column to the original 
# timeseries data

time_series_data$Basin_name <- NA
for (i in 1:nrow(time_series_data)) {
  print(i)
  for (j in 1:nrow(basin_name_code)) {
    if (time_series_data$HydroBasin[i] == basin_name_code$HydroBasin[j]){
      time_series_data$Basin_name[i] <- basin_name_code$Basin_name[j]
      
    }
  }
}

# Use this function to tag all species in a list of TimeSeries

fra.stat.matrices <- assign.stat.country(fra.ID.list, country = "FRA") # uses Tedesco basin data as well as country data.
gbr.stat.matrices <- assign.stat.country(gbr.ID.list, country = "GBR")
swe.stat.matrices <- assign.stat.country(swe.ID.list, country = "SWE")
usa.stat.matrices <- assign.stat.country(usa.ID.list, country = "USA")




###############################################################
#### Step 2. Computing invasive turnover metrics  #############
###############################################################

# Compute contribution of natives

fra.nnc.matrices <- mat.nnc.ass(fra.stat.matrices) # calculates relative contributions of natives and exotics.
gbr.nnc.matrices <- mat.nnc.ass(gbr.stat.matrices)
swe.nnc.matrices <- mat.nnc.ass(swe.stat.matrices)
usa.nnc.matrices <- mat.nnc.ass(usa.stat.matrices)



###############################################################
#### Step 3. Building a frame containing community metrics ####
###############################################################

fra.novel.mat <- inv.frame.builder(fra.nnc.matrices, country = "FRA") # this function returns the data frame for modelling.
gbr.novel.mat <- inv.frame.builder(gbr.nnc.matrices, country = "GBR")
swe.novel.mat <- inv.frame.builder(swe.nnc.matrices, country = "SWE")
usa.novel.mat <- inv.frame.builder(usa.nnc.matrices, country = "USA")



# Scale the data frame containing all countries (needs to be done after rbind..)
#test <- rbind(fra.novel.mat, gbr.novel.mat)

test$NNC <- scale(test$NNC, center = T, scale = T)
test$NNC_increase <- scale(test$NNC_increase, center = T, scale = T)
test$NAC <- scale(test$NAC, center = T, scale = T)
test$NAC_increase <- scale(test$NAC_increase, center = T, scale = T)
test$INC <- scale(test$INC, center = T, scale = T)
test$INC_increase <- scale(test$INC_increase, center = T, scale = T)
test$bin_lag <- scale(as.numeric(test$bin_lag, center = T, scale = T ))
test$position <- scale(test$position, center = T, scale = T)
test$basin <- as.factor(test$basin)
test$site <- as.factor(test$site)
test[is.na(test)] <- 0

#################################################################
#### Step 4. Modelling emergence of novelty by turnover ####
#################################################################

# I'm now thinking it makes more sense to use absolute abundance data

novel_model <- glmer(novel ~ bin_lag + position + NAC_increase* NNC_increase*INC_increase + (1|basin),
                     data = fra.novel.mat, 
                     family = binomial)


instant_model <- glmer(instant ~ bin_lag + position + NAC_increase*NNC_increase*INC_increase + (1|basin),
                       data = test, 
                       family = binomial)

cumul_model <- glmer(cumul ~ bin_lag + position + NAC_increase+NNC_increase+INC_increase + (1|basin),
                     data = test, 
                     family = binomial)

# Model Diagnostics
sim.res <- simulateResiduals(novel_model)
disp.test <- testDispersion(sim.res)
print(ifelse(disp.test$p.value <=0.05, 
             paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),
             "Dispersal okay :)"))

# Summarize model output
pred.df_a <- as.data.frame(summary(instant_model)$coefficients)
pred.df_a$taxa.rand <- summary(instant_model)$varcor$taxa[1,1]
pred.df_a$length.rand <- summary(instant_model)$varcor$length[1,1]

# Plot models"INC_increase [all]"
library(sjPlot)

plot_model(novel_model, type = "pred",  terms=c( "INC_increase [all]"), 
           title = "Probability of Novelty Emergence explained by Invader dynamics",
           axis.title = c("Net change in Invader Relative Abundance","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green")

# Nicer to do it in ggPlot
df_novel <- ggpredict(novel_model, type = "fe",  terms="INC_increase [all]")
df_instant <- ggpredict(instant_model, type = "fe" ,terms= "INC_increase [all]")
df_cumul <- ggpredict(cumul_model, type = "fe" ,terms= "INC_increase [all]")

# Plot INC holding all else equal
p1 <- ggplot(df_novel, aes(x, predicted)) +
  geom_line(color = "gold", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "yellow") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
          labs(y = "True Novel State (%)", x = "Invader Relative Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(novel_model)$coefficients)[6,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")
  
p2 <- ggplot(df_instant, aes(x, predicted)) +
  geom_line(color = "red1", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "red3") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Instantaneous Novel State (%)", x = "Invader Relative Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(instant_model)$coefficients)[6,4], 
                                                                           digits = 5), " *"), family = "Times New Roman")

p3 <- ggplot(df_cumul, aes(x, predicted)) +
  geom_line(color = "steelblue4", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "mediumturquoise") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Cumulative Novel State (%)", x = "Invader Relative Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(cumul_model)$coefficients)[6,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")

multiplot(p1, p2, p3)








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








