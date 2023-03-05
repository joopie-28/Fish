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
                 'betareg'
                 ))



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

# Load CRS
WG84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

###### PHASE 1 - PRE-PROCESSING AND MANIPULATION ######

#### Pre-processing 1. Construct community matrices and detect novelty ####

# These are the countries with sufficient data for bin range = 1 year, skipping countries with
# no usable time series saves us some redundant processing.

countries.suf.data <- list("BEL", "BRA", "BWA", "CAN", "CIV", "COL", 
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
                          nov.cluster.id.V7(nov.matrices[69],
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

# Create a persistence frame for length modelling

persistence_frame <- fullNovFrame_complete |>
  filter(novel.class == 'Persister' | novel.class == 'BLIP') |>
  filter(consistency == T) |>
  mutate(pers.proportion = ifelse(novel.class == 'Persister', novel.length/(total.n - position), 0)) |>
  # Transform to satisfy beta regression
  mutate(pers.prop.trans = (pers.proportion*(total.n - 1) + 0.5)/total.n,
         pers_binary = ifelse(novel.class == 'Persister', 1, 0),
         delta_spp = gain-loss,
         prop_delta_spp = delta_spp/orig.rich)


#### Pre-processing 8. Importing and extracting environmental drivers at time series level #####

# Extract all pfafstetter codes at levels 1-12 from the HydroBasins dataset.
# This is used for extracting environmental data at varying spatial scales
# in later analyses.
HYBAS_scheme <- create_basin_TS(time_series_data)

# This function extracts environmental data from the HydroAtlas, and adds it to our modelling frame.
EnviroByTS_L12 <- create_ENV_frame(fullNovFrame_complete, HYBAS_Level = 12, HYBAS_scheme, c('HYBAS_ID','run_mm_syr',
                                                                           'dis_m3_pyr', 'riv_tc_ssu',
                                                                           'dor_pc_pva',
                                                                           "crp_pc_sse", 'pst_pc_sse',
                                                                           'pac_pc_sse', 'hft_ix_s09','ppd_pk_sav'))


###### PHASE 2 - MODELLING AND ANALYSES ######

#### Modelling 1. Rates of Novelty Emergence around the globe ####

# Three separate models for each novelty type, at the community level
fixed.emergence.nov.mod <- glmer(novel ~ bin_lag + position+ (1|site_ID/Quarter), data = full.novel.mat.season, family= 'binomial')

fixed.emergence.cumul.mod <- glmer(cumul ~  bin_lag + position + (1|Quarter/site_ID), data = full.novel.mat.season, family= 'binomial')

fixed.emergence.instant.mod <- glmer(instant ~bin_lag + position + (1|Quarter/site_ID), data = full.novel.mat.season, family= 'binomial')

emergence.mod.list <- list(fixed.emergence.instant.mod, fixed.emergence.cumul.mod, fixed.emergence.nov.mod)

# One model at the time series level
broad.emergence.mod <- glmer(binary_novel ~ (1|Country), data = geo.timeseries.full, family = 'binomial')


#### Modelling 2. Modelling persistence length ####

# We use a beta regression to model persistence proportion. We use the transformed data.

summary(betareg(pers.prop.trans ~ n.from.end +loss + gain +shannon.d + evenness + INC_increase*total_inv , 
        data = persistence_frame))
# or
summary(glmmTMB(pers.prop.trans ~ n.from.end +loss + gain + evenness + INC_increase*total_inv + (1|site_ID/Quarter) , 
        data = persistence_frame, family='beta_family'))
# or
summary(glm(novel.length ~ position+loss + gain + evenness + INC_increase*total_inv , 
                data = persistence_frame, family='poisson'))

#### Modelling 3. Drivers of emergence #####

summary(glmer(novel~position+bin_lag + gain+loss +INC_increase*total_inv + (1|site_ID/Quarter), 
              data = fullNovFrame_complete, family = 'binomial'))

#### Modelling 4. Drivers of persistence #####
# add environmntal 
summary(glmer(pers_binary ~n.from.end +delta_spp+total_inv +(1|basin/site_ID/Quarter) , data = persistence_frame,
    family = 'binomial'))


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

# Write up results and the rest of the paper. That's it!




