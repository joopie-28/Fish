##################################################################
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
                 "gridExtra", "clustsig", "dendextend",
                 'betareg', 'car', 'visreg'
                 ))



# Load all the exotic species databases (constructed from literature and fishbase)
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

# Load CRS
WG84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

###### PHASE 1 - PRE-PROCESSING AND MANIPULATION ######

#### Pre-processing 1. Construct community matrices and detect novelty ####

# These are the countries with sufficient data for bin range = 1 year, skipping countries with
# no usable time series saves us some redundant processing.

countries.suf.data <- list("BEL", "BRA", "CAN", "CIV", "COL", 
                           "ESP", "FIN", "FRA", "GBR", "HUN", "JPN", 
                           "SWE", "USA", "AUS")

# Create a grand list of Time Series which can be used for analysis

full.ID.list <- do.call(c, lapply(countries.suf.data, function(country){
  test <- country_list_assigner(country)
  return(test)
}))

ID.list <- as.list(time_series_data$TimeSeriesID)

# Store all matrices in a named list
matrix_list_seasonality <- list_matrix_seasonality_function(ID.list)

# Apply the NDF to all matrices and store results in a new list
nov.output <- novelty.detection.gam(matrix_list_seasonality)

#### Pre-processing 2. Tag species as 'Invader', 'Non-native', or 'Native'####

full.stat.matrices.season <- assign.stat.country.V2(names(matrix_list_seasonality)) # please check anguilla australis australis


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
                        lapply(1:length(nov.matrices), function(m){
                          print(m)
                          nov.cluster.id.V7(nov.matrices[m],
                                            method_clus = 'single',
                                            alpha_clust = 0.05)
                        })) 

# Write the results back into our main data frame.
full.novel.mat.season$novel.length <- NA
full.novel.mat.season$novel.class <- NA
full.novel.mat.season$consistency <- NA
full.novel.mat.season$median_seq_dis <- NA
full.novel.mat.season$median_cum_dis <- NA
full.novel.mat.season$nspp <- NA
full.novel.mat.season$FirstPostNov <- NA

for(i in 1:length(novelty.pers)){
  indices <- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel")
  if(length(indices) > 1){
    indices<- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel" & full.novel.mat.season$bins == novelty.pers[[i]][[1]]$begin)
  }
  
  full.novel.mat.season$consistency[indices] <- novelty.pers[[i]][[1]]$consistency
  
  full.novel.mat.season$novel.length[indices] <- ifelse(is.null(novelty.pers[[i]][[1]]$length.bins), NA, novelty.pers[[i]][[1]]$length.bins)
  full.novel.mat.season$novel.class[indices] <-  ifelse(is.null(novelty.pers[[i]][[1]]$Class), NA, novelty.pers[[i]][[1]]$Class)
  full.novel.mat.season$median_seq_dis[indices] <- ifelse(is.null(novelty.pers[[i]][[1]]$median_seq_dis), NA, novelty.pers[[i]][[1]]$median_seq_dis)
  full.novel.mat.season$median_cum_dis[indices] <- ifelse(is.null(novelty.pers[[i]][[1]]$median_cum_dis), NA, novelty.pers[[i]][[1]]$median_cum_dis)
  full.novel.mat.season$nspp[indices] <- ifelse(is.null(novelty.pers[[i]][[1]]$nspp), NA, novelty.pers[[i]][[1]]$nspp)
  full.novel.mat.season$FirstPostNov[indices] <- ifelse(is.null(novelty.pers[[i]][[1]]$FirstPostNov), NA, novelty.pers[[i]][[1]]$FirstPostNov)
}



# Set non-novel lengths to 0 instead of NA
full.novel.mat.season$novel.length[is.na(full.novel.mat.season$novel.length)] <- 0
full.novel.mat.season$novel.class[is.na(full.novel.mat.season$novel.class)] <- "NONE"
full.novel.mat.season$novel.class <- as.factor(full.novel.mat.season$novel.class)

#### Pre-processing 6. Computing demographic processes and how they relate to novelty ####

# Create a new data frame with all turnover parameters
demography.frame <- rbindlist(lapply(1:length(matrix_list_seasonality), function(x){
  print(x)
  
  demo.df<-local.extorig.RivFishTime(ssmat=matrix_list_seasonality[[x]],novel=nov.output[[x]],
                                                   site.name = names(nov.output[x]))
  return(demo.df)
}))

# Add data to our main modelling frame
full.novel.mat.season<-cbind(full.novel.mat.season, demography.frame[,-c(1,2,3,5,6)])
full.novel.mat.season$basin <- as.factor(full.novel.mat.season$basin)


# Lastly, we filter out all time series with inconsistencies between SIMPROF and NDF,
# Following criteria outlined in the main document. We want to be sure that what is 
# detected by NDF is of a magnitude fit for reality.

LagConsistencyFilter <- function(full.novel.mat.season,nlag){
  
  for(site in full.novel.mat.season$site){
    print(site)
    # Subset the data by site
    sub <- full.novel.mat.season[which(full.novel.mat.season$site == site),]
    nov.lag <- sub$bin_lag[which(sub$cat == 'novel')]
    # Check if we have a consistency to inspect
    if(all(is.na(sub$consistency))){
      next
    }
    else{
      # If true, keep the data. We also remove those time series with lag of +3!
      if(na.omit(unique(sub$consistency)) & any(nov.lag < nlag)){
        next
      }
      else{
        # If false, the novel community is unreliable and we remove it and the time series
        print('Removing these rows')
        full.novel.mat.season <- full.novel.mat.season[-which(full.novel.mat.season$site == site),]
      }
    }
  }
  return(full.novel.mat.season)
}

full.novel.mat.season <- LagConsistencyFilter(full.novel.mat.season, nlag=4)

#### Pre-processing 7. Quantifying invaders at basin level ####

# Create a new data frame holding the number of encountered invasive species in a whole basin, per year.

InvByBasin <- rbindlist(lapply(unique(time_series_data$HydroBasin), 
                               function(basin){
                                 print(basin)
                                 countInv_By_Basin(time_series_data, basin)
                                 
                               }))

# Now need to add this on to the modelling data frame, probably with a join

fullNovFrame_complete <- merge(full.novel.mat.season, 
                               InvByBasin, by=c("basin","bins"))


#### Pre-processing 8. Importing and extracting environmental drivers at time series level #####

# Extract all pfafstetter codes at levels 1-12 from the HydroBasins dataset.
# This is used for extracting environmental data at varying spatial scales
# in later analyses.
HYBAS_scheme <- create_basin_TS(time_series_data)

environmental_variables = data.frame('Variable' = c("Natural_Discharge_Annual",
                                                    "Land_Surface_Runoff_Annual",
                                                    "Degree_Regulation",
                                                    "River_Area",
                                                    "Elevation",
                                                    "Terrain_slope",
                                                    "Stream_Gradient",
                                                    "Aridity_Index",
                                                    "Cropland_Extent",
                                                    "Pasture_Extent",
                                                    "Protected_Area_Extent",
                                                    "Population_Density",
                                                    "Urban_Extent",
                                                    "Human_Footprint"),
                                     'Code' = c('dis_m3_pyr',
                                                "run_mm_syr",
                                                "dor_pc_pva",
                                                "ria_ha_ssu",
                                                "ele_mt_sav",
                                                "slp_dg_sav",
                                                "sgr_dk_sav",
                                                "ari_ix_sav",
                                                "crp_pc_sse",
                                                "pst_pc_sse",
                                                "pac_pc_sse",
                                                "ppd_pk_sav",
                                                "urb_pc_sse",
                                                "hft_ix_s09"))

saveRDS(HYBAS_scheme, file = "./outputs/HYBAS_Scheme_1-12.rds")

#### Pre-processing 9. Creating a site level base dataframe ####

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

# This function extracts environmental data from the HydroAtlas, and adds it to our modelling frame.
# We opted to go one step further and extract vairables for individual river segments, however
# there is not a massive difference when considering envrionemntal variables at level 12.
EnviroByTS_L12 <- create_ENV_frame(geo.timeseries.full, HYBAS_Level = 12, 
                                   HYBAS_scheme, environmental_variables$Code) |>
  mutate("binary_novel" =  ifelse(novel > 0, 1, 0))


#### Pre-processing 10. Creating finished data frames for modelling at site and community level ####

# Import spatial data for all relevant areas (HydroRivers database)
RiverID_eu <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")
RiverID_aus <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_au_shp/HydroRIVERS_v10_au_shp")
RiverID_af <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af_shp")
RiverID_sa <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_sa_shp/HydroRIVERS_v10_sa_shp")
RiverID_na <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp")
RiverID_as <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as_shp")

# Extract relevant rivers from the HydroRivers database
RiverID_Global <- rbindlist(list(RiverID_eu, RiverID_af, RiverID_as, 
                                 RiverID_aus, RiverID_na, RiverID_sa))

# This may take a while - we are working large volumes of spatial data.
# a completed file can be found at './outputs/globalRivers_Extracted.rds'
# and can be used for downstream analyses.

globalRivers_Extracted <- rbindlist(lapply(countries.suf.data, function(country){
  print(country)
  ExtractRiverFromSurvey(EnviroByTS_L12, country, RiverID_Global)
}))

# This was a very intensive computation, due to the size of the HydroRiver Database,
# so we will be saving the obejct in an .RDS file.

saveRDS(globalRivers_Extracted, file = './outputs/globalRivers_Extracted.rds')


# Now we need to import all the environmental variables.

AUSRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_au")
EURiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                        lay = "RiverATLAS_v10_eu")
NARiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                        lay = "RiverATLAS_v10_na")
SANRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_sa_north")
SASRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_sa_south")
AFRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                        lay = "RiverATLAS_v10_af")
ASIRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_as")

# Declare our list of variables that we are interested in: doing this here saves some computing time.
environmental_variables_river = data.frame('Variable' = c("Natural_Discharge_Annual",
                                                          "Land_Surface_Runoff_Annual",
                                                          "Degree_Regulation",
                                                          "River_Area",
                                                          "Elevation",
                                                          "Terrain_slope",
                                                          "Stream_Gradient",
                                                          "Aridity_Index",
                                                          "Cropland_Extent",
                                                          "Pasture_Extent",
                                                          "Protected_Area_Extent",
                                                          "Population_Density",
                                                          "Urban_Extent",
                                                          "Human_Footprint"),
                                           'Code' = c('dis_m3_pyr',
                                                      "run_mm_cyr",
                                                      "dor_pc_pva",
                                                      "ria_ha_csu",
                                                      "ele_mt_cav",
                                                      "slp_dg_cav",
                                                      "sgr_dk_rav",
                                                      "ari_ix_cav",
                                                      "crp_pc_cse",
                                                      "pst_pc_cse",
                                                      "pac_pc_cse",
                                                      "ppd_pk_cav",
                                                      "urb_pc_cse",
                                                      "hft_ix_c09")) 

# Create this massive list with all environmental variables for all possible rivers.
data_list <- list(AUSRiver_Data, EURiver_Data, NARiver_Data , 
                  SANRiver_Data, SASRiver_Data, AFRiver_Data,
                  ASIRiver_Data)

# Filter out rows and variables to create a workable amount of data
for (i in 1:length(data_list)){
  print(i)
  data_list[[i]] <- data_list[[i]] |> 
    select(c(1:14,environmental_variables_river$Code)) |>
    filter(HYRIV_ID %in% globalRivers_Extracted$HYRIV_ID)
}

# Collate into a large dataframe and clean up.
RiverData_Global <- rbindlist(data_list)
rm(data_list)

# Save it to the environemnt so we do not have to recompute these data.
saveRDS(RiverData_Global, file = './outputs/RiverData_Global.rds')

# Join all data, such that we have environmnetal variables
# linked with a timeseries ID.
MergedENV_by_RivID <- globalRivers_Extracted |>
  left_join(RiverData_Global, by = 'HYRIV_ID') |>
  rename(ID = site)

# Join all environmental data to our time series 
# data frame; at both site and community levels.

# This is the site-level data frame CONSISTENCY?????
FullGeoFrame <- geo.timeseries.full |>
  left_join(MergedENV_by_RivID, by = ("ID")) |>
  mutate(binary_novel = ifelse(novel > 0, 1, 0))|>
  separate(ID, c('Site', 'Quarter'), remove = F) 

# This is the community-level data frame
FullEnvFrame <- fullNovFrame_complete |>
  rename(ID = site) |>
  left_join(MergedENV_by_RivID, by = c('ID'))

# This frame includes only novel communities and is used
# for persistence length modelling
PersistenceFrame <- FullEnvFrame|>
  filter(novel.class == 'Persister' | novel.class == 'BLIP') |>
  filter(consistency == T) |>
  mutate(binary_pers = ifelse(novel.class == "Persister", 1, 0))



###### PHASE 2 - MODELLING AND ANALYSES ######

#### Modelling 1. Rates of Novelty Emergence around the globe ####

# Three separate models for each novelty type, at the community level
fixed.emergence.nov.mod <- glmer(novel ~ bin_lag + position+ (1|site_ID/Quarter), 
                                 data = FullEnvFrame, family= 'binomial')

fixed.emergence.cumul.mod <- glmer(cumul ~  bin_lag + position + (1|Quarter/site_ID), 
                                   data = FullEnvFrame, family= 'binomial')

fixed.emergence.instant.mod <- glmer(instant ~ bin_lag + position + (1|Quarter/site_ID), 
                                     data = FullEnvFrame, family= 'binomial')

# One model at the time series level
broad.emergence.mod <- glmer(binary_novel ~ (1|Basin), 
                             data = FullGeoFrame, family = 'binomial')

#### Modelling 2. Modelling persistence length and chance of blip versus persistant state ####

# Simple poisson glm to understand the variables associated with persistence time.
persLengthMod <- glmer(novel.length ~ n.from.end +loss + gain + evenness + INC_increase*total_inv + (1|site_ID/Quarter) , 
                data = PersistenceFrame, family='poisson')

# Intercept-only random effects model to get estimate for the proportion of persisters v non-persisters.
persBinaryNullMod <- glmer(binary_pers ~ 1 + (1|site_ID/Quarter) , data = PersistenceFrame,
              family = 'binomial')

# Binomial regression to understand factors that contribute to whether or not a community persists at all.
persBinaryFullMod <- glmer(binary_pers ~ n.from.end+total_inv +(1|site_ID/Quarter) , data = PersistenceFrame,
              family = 'binomial')

#### Modelling 3. Drivers of emergence #####
test<-FullEnvFrame %>% 
  mutate_at(c('position','loss', 'gain','bin_lag' ,"run_mm_cyr",'delta_eveness',"ari_ix_cav", "INC_increase","total_inv", "DIST_DN_KM","DIST_UP_KM","run_mm_cyr", 'ppd_pk_cav','hft_ix_c09',"DIST_DN_KM","DIST_UP_KM",
              "dis_m3_pyr",'ria_ha_csu',"ele_mt_cav","urb_pc_cse", "crp_pc_cse", "pac_pc_cse", "pst_pc_cse"), 
            ~(scale(., center =T, scale =T) %>% as.vector)) 
 
# Binomial glmm looking at community level emergence using ecological and environemntal variables
EmergCommMod <- glmer(cumul~position + gain+loss+delta_eveness +INC_increase*total_inv+ dis_m3_pyr +DIST_DN_KM+run_mm_cyr*ppd_pk_cav+
                        ele_mt_cav*DIST_DN_KM+(1|site_ID/Quarter), 
              data = test, family = 'binomial')

forestPlotter <- function(model, labels){
  
  # Set theme
  set_theme(
    geom.outline.color = "antiquewhite4", 
    geom.outline.size = 1, 
    geom.label.size = 2,
    title.size = 1.5, 
    axis.angle.x = 0, 
    axis.textcolor = "black",
    axis.linecolor.x = 'black',
    axis.linecolor.y = 'black',
    base = theme_classic(),
    title.color = 'white'
      
  )
  # plot model
  plot_model(model, show.p = T, axis.labels = NULL,
             col = c('red', 'steelblue'), vline.color = 'red')
  
}

forestPlotter(EmergCommMod)

# Binomial glm looking at site-level emergence proportions using only environmental variables,
# essentially inspecting whether or characteristics of sites are associated with novelty.
test.geo<-FullGeoFrame %>% 
  mutate_at(c("run_mm_cyr", 'ppd_pk_cav','hft_ix_c09',"DIST_DN_KM","DIST_UP_KM","dis_m3_pyr",'ria_ha_csu',"ele_mt_cav","urb_pc_cse", "crp_pc_cse", "pac_pc_cse", "pst_pc_cse"), 
            ~(scale(., center =T, scale =T) %>% as.vector)) |>
  filter()
EmergSiteMod <- glm(binary_novel ~ dis_m3_pyr + run_mm_cyr + ele_mt_cav +DIST_DN_KM+
                      ppd_pk_cav+pac_pc_cse*hft_ix_c09 + ele_mt_cav*DIST_DN_KM, 
              data = test.geo, family = 'binomial') 



# Modelling 4. First state post novelty.

# Simply looking at the averages and if they differ at all.
test <- FullEnvFrame |>
  mutate(binary_postnov = ifelse(FirstPostNov == 'New_exploratory_State', 1, 0)) |>
  filter(novel.class == 'Persister' | novel.class == 'BLIP')

summary(glm(binary_postnov~position+novel.class , data = test, family = "binomial")) # Is this a result? consult Simon.

#### Combine all models in a list for ease of analysis ####

NovelFishModels <- list('Emergence Rates' = list("fixed.emergence.nov.mod" = fixed.emergence.nov.mod,
                              "fixed.emergence.cumul.mod" = fixed.emergence.cumul.mod,
                              "fixed.emergence.instant.mod" = fixed.emergence.instant.mod,
                              "broad.emergence.mod" = broad.emergence.mod),
     'Persistence Models' = list("persLengthMod" = persLengthMod,
                                 "persBinaryNullMod" = persBinaryNullMod,
                                 "persBinaryFullMod" = persBinaryFullMod),
     'Emergence Drivers' = list("EmergCommMod" = EmergCommMod,
                                "EmergSiteMod" = EmergSiteMod))

###### PHASE 3 - MODEL DIAGNOSTICS AND FILE EXPORT ######

# Test for overdispersion in each model using DHARma package.

CheckDispersion(models)

# vif

# Create output csv's for all models using a custom function.
lapply(1:length(flatten(NovelFishModels)), function(x){
  extract.coefs(flatten(NovelFishModels)[x])
})

###### PHASE 4 - PRODUCING TABLES AND FIGURES ######

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
extract.coefs(broad.emergence.mod)

#### Table 2.

#### Table 3.

#### Table 4.
# To do #
# Complete this document

# Write up results and the rest of the paper. That's it!




