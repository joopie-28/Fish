#### Characterizing the drivers of novel freshwater fish communities #######

############################################################################
### This script includes main analyses as well as data import & tidying ####
### it requires functions to be loaded in from the function folder, as well 
### as access to certain databases. 


# J.M. Sassen 31-01-2022 

##################################################
#### Step 0 . Load in packages and databases #####
##################################################

# Clear environment (if wanted) and set your working directory
rm(list = ls())

# source functions from 'functions' sub-folder and load them
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

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
                 "gridExtra", "clustsig", "dendextend"))

# Load data and make sure names are in orderly fashion.
time_series_data <- read.csv("./inputs/1873_2_RivFishTIME_TimeseriesTable.csv") # RIVFishTime 
Survey_Data <- read.csv("./inputs/1873_2_RivFishTIME_SurveyTable.csv") # RIVFishTime

invasives_data <- read.csv2("./inputs/datatoFigshare/Occurrence_Table.csv") # Tedesco et Al. Invasives database
colnames(invasives_data) <- c("Basin", "Species", "Status", "TSN_ITIS_Code", "FishBase_code", "Valid_name", "Sighting") # Tedesco et Al. Invasives database

basin_data <- read.csv2("./inputs/datatoFigshare/Drainage_Basins_Table.csv") # Tedesco et Al. Basins geodatabase
colnames(basin_data) <- c("Basin", "Country", "BioRealm", "Endorheic", "Lon", "Lat", "Med_Lon", "Med_Lat", "Surface_Area") # Tedesco et Al. Basins geodatabase

occurence_shapefile <- read_sf(dsn = "./inputs/datatoFigshare/Basin042017_3119.shp") # Tedesco et AL. Shapefile


# Tidy Species names to prevent errors

invasives_data$Species <- gsub(invasives_data$Species, 
                               pattern ="\\.", 
                               replacement = " ")

invasives_data$Valid_name <- gsub(invasives_data$Valid_name, 
                                  pattern ="\\.", 
                                  replacement = " ")

# Set up the HydroBasin names and codes

## NOTE: We are likely to not assign on basin level as basin level data is 
## not uniformly available.

# This function generates a data frame matching a HydroBasin code to
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


##################################################
#### Step 1. Tagging species status  #############
##################################################

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

# This function takes in all species abundance matrices, and tags each
# species featured in those matrices as a native or an invader.

full.stat.matrices <- assign.stat.country_nn(full.ID.list) 

full.stat.matrices.season <- assign.stat.country.V2(names(matrix_list_seasonality))

###############################################################
#### Step 2. Computing invasive turnover metrics  #############
###############################################################

# Computes certain ecological metrics of invaders and natives

full.nnc.matrices <- mat.nnc.ass(full.stat.matrices.season)
full.nnc.matrices.season <- mat.nnc.ass.V2(full.stat.matrices.season)
###############################################################
#### Step 3. Building a frame containing community metrics ####
###############################################################

# This function creates a all-encompassing data frame for
# suitable for use in modelling 

full.novel.mat <- inv.frame.builder(full.nnc.matrices)
full.novel.mat.season <- inv.frame.builder.V2(full.nnc.matrices.season)

temp_df <- data.frame(do.call("rbind", strsplit(as.character(full.novel.mat.season$site), ".",
                                                fixed = TRUE)))

names(temp_df) <- c("site_ID", "Quarter")

full.novel.mat.season<-cbind(temp_df, full.novel.mat.season)
###########################################################################
#### Step 4. Modelling emergence of novelty by invasives turnover #########
###########################################################################

# Create a df that excludes post-novel communities as this will bias the result.
# We are only interested in the first part
temp_list<-NULL
for (site.ID in unique(full.novel.mat.season$site)){
  print(site.ID)
  temp<-subset(full.novel.mat.season, site == site.ID)
  if(any(temp$cat == 'novel')){
   index <- which(temp$cat =='novel')[1]
   temp<- temp[1:index,]
  }
  else{
    temp<-temp
  }
  temp_list <- c(temp_list, list(temp))

}
inv.df<-rbindlist(temp_list)




# This function runs GLMM's for the invasives-novelty investigations. If plot = T,
# Additional plots of the invader effect on novelty probability are saved. 
# CSV's with model output are simultaneously created

inv.models <- inv.nov.glmm(full.novel.mat, plot = T)


inv.2 <- subset(full.novel.mat.season, country == "FRA" )
full.novel.mat.season$test <- abs(full.novel.mat.season$INC_spec_increase)
inv.df$position <- scale(inv.df$position, center = T, scale = T)
inv.df$bin_lag <- scale(inv.df$bin_lag, center = T, scale = T)
inv.df$INC_increase <- scale(inv.df$INC_increase, center = T, scale = T)
inv.df$INC <- scale(inv.df$INC, center = T, scale = T)
inv.df$NAC <- scale(inv.df$NAC, center = T, scale = T)

mod<-glm(novel~position+bin_lag+test,
            data =full.novel.mat.season , family ="binomial")

mod<-glm(novel~position+bin_lag+INC_spec+NNC_spec+NAC_spec,
         data =full.nov.no.lag, family ="binomial")




summary(mod)

visreg(mod.list[[1]], 'BioRealm' ,scale = 'response', rug=F, ylim = c(0,1))
visreg(mod.INV, 'INC_increase' ,scale = 'response', rug=F, ylim = c(0,1))
visreg(environ.driver.mod, 'Anthromod' ,scale = 'response', rug=F, ylim = c(0,.2))

points(novel ~Anthromod , data =geo.timeseries.full, pch=19, col = alpha('black', alpha=0.4), cex=1)


# Random allocation of species as a null model




#################################################################
#### Step 5. Investigating Novelty Persistence  #################
#################################################################

# Here, we use ANOSIM to further analyse what happens to a novel state
# after it emerges.

# We now have a method of 1) visualizing differences between communities in two-dimensional
# space and 2) determining whether the difference between pre- and post-novel states is greater
# than the difference within these groups.

# Now, we want to use these techniques to discover how long a novel state is maintained (in years).
# Because our ANOSIM test simply computes differences between groups, there might be some communities
# that are in fact quite different from the novel community (i.e. on their way back to a pre-novel state)
# in the post-novel group.


# Filter matrices where novelty occurred
test<-full.novel.mat.season$site[which(full.novel.mat.season$cat == 'novel')]
nov.matrices<-matrix_list_seasonality[full.novel.mat.season$site[which(full.novel.mat.season$cat == 'novel')]]

# Run persistence framework

novelty.pers <- do.call(c, 
                        sapply(1:length(nov.matrices), function(x){
                          print(x)
                          nov.cluster.id.V6(nov.matrices[x])
                          
                        })) 

# Add the calculated lengths to our main data frame
full.novel.mat.season$novel.length <- NA
full.novel.mat.season$novel.class <- NA

for(i in 1:length(novelty.pers)){
  indices <- which(full.novel.mat.season$site_ID == novelty.pers[[i]][[1]]$ID & full.novel.mat.season$cat == "novel")
  if(length(indices)> 1){
   indices<- which(full.novel.mat.season$site_ID == novelty.pers[[i]][[1]]$ID & full.novel.mat.season$cat == "novel" & full.novel.mat.season$bins == novelty.pers[[i]][[1]]$begin)
  }
  
  full.novel.mat.season$novel.length[indices] <- novelty.pers[[i]][[1]]$length 
  full.novel.mat.season$novel.class[indices] <-  novelty.pers[[i]][[1]]$Class
  
}

# Set non-novel lengths to 0 instead of NA
full.novel.mat.season$novel.length[is.na(full.novel.mat.season$novel.length)] <- 0
full.novel.mat.season$novel.class[is.na(full.novel.mat.season$novel.class)] <- "NONE"

### Investigate percentage of timeseries spent in novel state post-emergence ###
# This looks at the proportion of a time series spent in a novel state post-emergence
proportion.persistence.calculator <- function(full.novel.mat){
  
  nov.duration.df <- c('site'=c(NULL), 'proportion' = c(NULL))
  for (site2 in unique(subset(full.novel.mat, cat =='novel')$site)){
    print(site2)
    test <- subset(full.novel.mat, site == site2)
    
    nov.ind<-which(test$cat == 'novel')[1]
    nov.pos<-test$bins[nov.ind]
    nov.end.pos<-nov.pos-test$novel.length[nov.ind]
    end.pos <- test$bins[nrow(test)]
    
    nov.proportion<-(nov.pos-nov.end.pos)/(nov.pos-end.pos)

    nov.duration.df$site <- c(nov.duration.df$site, site2)
    nov.duration.df$proportion <- c(nov.duration.df$proportion, nov.proportion)
      
    }
    
  return(as.data.frame(nov.duration.df))
}

persistence.length.df <- proportion.persistence.calculator(full.novel.mat)

full.novel.mat$persistence.proportion <- NA
for(i in 1:nrow(persistence.length.df)){
  print(i)
 indices<-which(full.novel.mat$site == persistence.length.df$site[i])
 full.novel.mat$persistence.proportion[indices] <- persistence.length.df$proportion[i]
}

for(i in 1:nrow(full.novel.mat)){
  if(full.novel.mat$novel.class[i] == "BLIP"){
    full.novel.mat$persistence.proportion[i] <- 0
  }
}

##############################################################################
### Step 6. Understanding lag and missing data with respect to persistence ###
##############################################################################

# We wish to understand the role of missing data on the persistence classification
# of the novel community (i.e. longer lag, more likely to be blip?)

# Add a new feature representing the lag between the novel community T and the 
# succeeding community T + 1.

full.novel.mat$lag_to_next <- NA

# Compute the lags for each community (slightly different to earlier defined bin lag variable)
for (site.ID in unique(full.novel.mat$site)){
  print(site.ID)
  temp <- subset(full.novel.mat, site == site.ID)
  indices <- which(full.novel.mat$site == temp$site)
  for (i in 1:nrow(temp)){
    temp$lag_to_next[i] <- abs(temp$bins[i+1] - temp$bins[i])
  }
  full.novel.mat$lag_to_next[indices] <- temp$lag_to_next
}

# Define the model
lag.mod.data <- subset(full.novel.mat, novel.class != "END" & novel.class != "NONE" )
lag.mod <- glm(novel.class.binary~lag_to_next, data = lag.mod.data, family = 'binomial')

data.2<- subset(full.novel.mat.test, novel.class != "Persister")
visreg(immig.mod, 'transition', scale = "response", partial = F, rug = F, ylim = c(0,.5),
       ylab = 'persistence type (0 = blip, 1 = persistent)', xlab = 'lag between novel and succeeding community (years)')
points(novel.class.binary~lag_to_next, data = lag.mod.data, pch=19, col = alpha('black', alpha=0.4), cex=1)

hist(lag.mod.data$lag_to_next, xlim = c(0,15), xlab  = 'lag between novel and succeeding bin', main = '')

nrow(subset(full.novel.mat,lag_to_next < 3 & cat =='novel'))
##########################################################
### Step 7. Varying Tau for persistence classification ###
##########################################################

# Convert the novel length in years to a time bin value 
# (doing this here as it was not part of original analysis)
full.novel.mat$novel.length.bins <- 0

for (site.ID in unique(subset(full.novel.mat, cat == "novel" & novel.class == 'Persister')$site)){
  print(site.ID)
  temp <- subset(full.novel.mat, site == site.ID)
  indices <- which(full.novel.mat$site == temp$site)
  nov.ind <- which(temp$cat == 'novel')[1]
  
  end.bin <- temp$bins[nov.ind] - temp$novel.length[nov.ind]
  end <- temp$position[which(temp$bins == end.bin)]
  start <- temp$position[nov.ind]
  
  novel.length.bin <- end-start
  
  full.novel.mat$novel.length.bins[which(full.novel.mat$cat == 'novel' & full.novel.mat$site == site.ID)] <- novel.length.bin
  
}

# Now create a new classification for each level of Tau

# Probably turn this into an apply function of sorts,
# just tricky with the separate dataframe dependency.

tau = 3
full.novel.mat$class.T3 <- NA

for (i in 1:nrow(full.novel.mat)){
  print(i)
  if(full.novel.mat$novel.class[i] == "BLIP"){
    full.novel.mat$class.T3[i] <- "BLIP"
  }
  if(full.novel.mat$novel.class[i] == "NONE"){
    full.novel.mat$class.T3[i] <- "NONE"
  }
  if(full.novel.mat$novel.class[i] == "END"){
    full.novel.mat$class.T3[i] <- "END"
  }
  if(full.novel.mat$novel.class[i] == "Persister"){
    if(full.novel.mat$novel.length.bins[i] > tau){
      full.novel.mat$class.T3[i] <- "Persister"
    }else{
      full.novel.mat$class.T3[i] <- "BLIP"
    }
  }
}


######################################################################
### Step 8. Model persister probability as an effect of demography ###
######################################################################

# Binary classification of persistent (1) versus blips (0)
for (i in 1:nrow(full.novel.mat)){
  print(i)
  if(full.novel.mat$novel.class[i] == 'Persister'){
    full.novel.mat$novel.class.binary[i] <- 1
  }
  else{
    full.novel.mat$novel.class.binary[i] <- 0
  }
}

# Model the effect of origination, extinction and invaders on probability of persistent or not.
data.nov=subset(full.novel.mat,novel.class != 'END'& novel.class != 'NONE')
data.nov$shannon.d <- scale(data.nov$shannon.d , center = T, scale = T)
data.nov$orig <- scale(data.nov$orig , center = T, scale = T)
data.nov$ext <- scale(data.nov$ext , center = T, scale = T)

persistent.trajectory.mod<- glm(novel.class.binary~orig+ext+shannon.d , data=data.nov, family = 'binomial')

summary(persistent.trajectory.mod)

# Plot invader effect holding all else equal
demo.persistent.plots <- function(mod){
  layout(matrix(1:2, nrow=1))

  par(mar=c(4,4,2,0.5))
  visreg(mod, 'orig', scale='response', partial=F, rug=F, 
         bty='l', xlab = 'Local Originations during Emergence', ylab = 'Persistence Probability of Novel Community',
         points = list(col='black'), line = list(col = 'green'), ylim = c(0,1), axes=F)
  axis(side = 1)
  axis(side = 2)
  points(novel.class.binary~orig, data = data.nov, pch=19, col = alpha('black', alpha=0.6), cex=0.4)
  box()
  mtext("c)", side = 3, line=0.5, adj = c(-4,0))
 

  par(mar=c(4,0.5,2,4))
  visreg(mod, 'ext', scale='response', partial=F, rug=F, 
         bty='l', xlab = 'Local Extinctions during Emergence', ylab = 'Persistence Probability of Novel Community',
         points = list(col='black'), line = list(col = 'red'), ylim = c(0,1), axes=F)

  axis(side = 4)
  axis(side = 1)
  points(novel.class.binary~ext, data = data.nov, pch=19, col = alpha('black', alpha=0.6), cex=0.4)
  box()
  mtext("d)", side = 3, line=0.5, adj = c(-4,0))

}


demo.persistent.plots(persistent.trajectory.mod)


### Zero inflated alternative model ####
summary(zeroinfl(novel.class.binary~orig+ext+shannon.d +orig, data=data.nov))


# Alternatively fit a glmm with cbind binary variable

full.novel.mat$all_comm <- paste0(full.novel.mat$cat, "-", full.novel.mat$novel.class)
full.novel.mat$bin_to_end <- full.novel.mat$total.n - full.novel.mat$position

data.nov=subset(full.novel.mat, novel.class != 'END')
mod<-(glmer(cbind(ext, (diversity.previous)) ~ all_comm+ bin_lag+bin_to_end+ (1|basin), data=data.nov, family = 'binomial'))
mod<- glm(cbind(orig, (diversity.next)) ~ all_comm +bin_lag+ (position), data=data.nov, family = 'binomial')

visreg(mod, c('all_comm'), scale='response', partial=F, rug=F, 
       bty='l', xlab = 'Local Originations during Emergence', ylab = 'Persistence Probability of Novel Community',
       points = list(col='black'), line = list(col = 'green'), ylim = c(0,.15), axes=T)


#############################################
### Inspecting seasonality  #################
#############################################





#### End of Main Analyses ####


#############################################
### Modelling distance after emergence ######
#############################################

matrix_test <- matrices[[1]][["G7715"]]

cluster.assigner <- function(matrix_test, label_frame){
  if(!any(label_frame$cat == 'novel')){
    matrix_test$cluster_type <- 'non.novel'
  }else{
    test<-simprof(data = matrix_test, num.expected = 1000, undef.zero = TRUE,
                  num.simulated = 999, method.distance =vegdist, 
                  method.cluster = "average", alpha=0.05)
    
    names.vec <- unlist(list(test$significantclusters))
    new.vec <- NULL
    for(i in 1:length(names.vec)){
      for(j in 1:length(test$significantclusters)){
        if (names.vec[i] %in% test$significantclusters[[j]]){
          new.vec[i] <- j
          names(new.vec)[i] <- names.vec[i]
        }
      }
    }
    sorted.index <- as.character(sort(as.numeric(names(new.vec)), decreasing  =T))
    new.vec <- new.vec[sorted.index]
    matrix_test$cluster <- new.vec
    nov.bin<-label_frame$bins[which(label_frame$cat == 'novel')][1]
    matrix_test$cluster[which(rownames(matrix_test) == nov.bin)]
    nov.cluster <- matrix_test$cluster[which(rownames(matrix_test) == nov.bin)]
    clus.list<- NULL
    matrix_test$cluster_type <- "pre.novel"

    for(i in 1:nrow(matrix_test)){
      if(matrix_test$cluster[i] != nov.cluster & nov.cluster %!in% clus.list){
        matrix_test$cluster_type[i] <- 'pre.novel'
        clus.list <- c(clus.list,matrix_test$cluster[i])
      }
      if(matrix_test$cluster[i] %!in% clus.list & matrix_test$cluster[i] != nov.cluster & nov.cluster %in% clus.list){
        matrix_test$cluster_type[i] <- 'post.novel'
      }
      if(matrix_test$cluster[i] == nov.cluster){
        matrix_test$cluster_type[i] <- 'novel'
        clus.list <- c(clus.list,matrix_test$cluster[i])
      }
  
    }
  }
  return(matrix_test$cluster_type)
}


cluster.assigner(matrix, label_frame) 












################################################
#### Step X. Analyzing state-level metrics ####
###############################################

### deprecated ###

# Determination of lengths allows us to look at characteristics of the novel state as it evolves.
# We have an unbiased method that tells us when the state ends, so we know what it "novel" and
# what is not.

# This functions extract information from the novelty results and returns a data frame with summary
# statistics. It also returns a data frame with 'state-level metrics'; a summary of important 
# biological metrics aggregated at distinct periods within each time frame.

# NEED TO UPDATE THIS UNFINISHED PER 16-06-2022
# 11/08/22 still need to update this function
# to reflect changes in length calculatror.

nov.comm.summary <- novel.comm.analyzer(full.novel.mat)

# Save summary stats to a CSV file

write.csv(nov.comm.summary$summary, 
          "./outputs/nov_summary_stats.csv")


#### End of Main Analyses ####



##########################################
### Supp. Analysis 1: Simulated ANOSIM ###
##########################################

# 1. Ecological Blip

# Simulate a one-off novelty emergence followed by
# an immediate return to 'normality'

matrices <- generate.matrices("blip", 1)

# 2. Short persistence

# Simulate novelty emergence followed by a short
# persistence period, prior to reverting

matrices <- generate.matrices("short.pers", 1)

# 3 Full persistence

matrices <- generate.matrices("full.pers", 1)

# 4. Slow decay, the most random of simulations

matrices <- generate.matrices("slow.decay", 1)


#### Hierarchical clustering to identify novel community persistence ####

# Function to check for sequentiality
is.sequential <- function(x){
  all(diff(x) == diff(x)[1])
}

# want to change this to maybe all the communities
# anyhow, persistence should be just that - no subcategories

nov.cluster.id.V1 <- function(matrix){

  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  number.names <- as.numeric(rownames(matrix))[-c(1:5)]
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                         alpha = 0.05,
                         metric = "bray",
                         plot =F, 
                         site = ID,
                         plot.data = FALSE,
                         gam.max.k = -1)
  
 
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  novel_frame <- novel_frame[1,]
  
  # Find year and row data
  novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix[-c((nrow(matrix)-4):nrow(matrix)),]
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance ="czekanoswki", 
            method.cluster = "average", alpha=.05), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.05) })
  
  # Initialize parameters for the length calculations
  run <- TRUE
  no.clus <- 2 # always start on 2
  
  # Plot dendogram result
  par(mfrow = c(1,1))

  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <- FALSE
    # This will activate the BLIP module and assign the correct 
    # category
    length <- 1
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    test<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(test, hang=-1,
         main = "Blip Dendrogram", xlab = "")
    })

  # Start of the length calculation module
  
  while(run == TRUE){
    
    # Length counter module, reverse such that oldest is last
    trajectory <- rev(cutree(test$hclust, k=no.clus))
    
    # Find group allocation
    novel.group <- trajectory[[which(names(trajectory)==novel_bin)]]
    
    # Count consecutive length of novel group
    novel.slice <- (which(trajectory == novel.group))
    
    # Find start and end points
    extrema <- as.numeric(names(novel.slice)[c(1,length(novel.slice))])
    
    # Check that start corresponds to novel bin. If False,
    # increase number of supersets to better reflect strutcture.
    if(novel_bin != extrema[1]){
      no.clus <- no.clus +1
    }else{
      
      start <- extrema[1]
      
      # Check that numbers are sequential using custom 
      # function 
      
      if(is.sequential(novel.slice)){
        end <- extrema[2]
      }
      
      # If not sequential, adjust the end year accordingly
      else{
        bool.vec <- diff(novel.slice) == 1
        index <- which(bool.vec == FALSE)[1] - 1
        # correct for 0 index during blips
        end <- ifelse(index > 0,as.numeric(names(bool.vec[index])),novel_bin)
      }

      # Calculate Length
      length <- start-end
      
      # Solve for blips, as we want to set the time in years, 
      # not timebins
      if(length == 0){
        length <- 1
      }
      run <- FALSE
    }
    
    # Use SIMPROF to assess significance of clusters. Need 
    # to make sure that coarser groups are supersets of signif
    # divisions.
    
  }
  
  # Initialize a class variable
  Class <- "NONE"

  # Allocate blips
   if(length == 1){
  
      Class <- "BLIP"
      end <- as.numeric(rownames(matrix.temp[novel_index-1,]))
      length <- novel_bin - end
   }

  # Test for full persistence
  if(Class != "BLIP"){
  
    index.vector <- which(names(novel.slice) %in% names(trajectory[novel.slice[1]:length(trajectory)]))
    
    if (length(trajectory[novel.slice[1]:length(trajectory)]) == length(index.vector)){
    
      Class <- "Full Persistence"
    }
  }
  # Rest are short persisters
  if(Class != "BLIP" & Class != "Full Persistence"){
    
    Class <- "Short Persistence"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  


 # Return data
 return.data <- list(c(list("ID" = ID, "begin" =start, 
                     "end" =end, "length"= length, 
                     "Class" = Class, "clusters" =no.clus)))
 names(return.data) <- ID
 return(return.data)
 
}

# Updating this such that all communities are considered, makes more sense.
# Also include hclust check


nov.cluster.id.V2 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  number.names <- as.numeric(rownames(matrix))
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  novel_frame <- novel_frame[1,]
  
  # Find year and row data
  novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance ="czekanoswki", 
            method.cluster = "average", alpha=.05), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.05) })
  
  # Initialize parameters for the length calculations
  run <- TRUE
  no.clus <- 2 # always start on 2
  
  # Plot dendogram result
  par(mfrow = c(1,1))
  
  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <- FALSE
    # This will activate the BLIP module and assign the correct 
    # category
    length <- 1
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(clust.data, hang=-1,
         main = "Blip Dendrogram", xlab = "")
  })
  
  # Start of the length calculation module
  
  while(run == TRUE){
    
    # Length counter module, reverse such that oldest is last
    trajectory <- rev(cutree(test$hclust, k=no.clus))
    
    # Find group allocation
    novel.group <- trajectory[[which(names(trajectory)==novel_bin)]]
    
    # Count consecutive length of novel group
    novel.slice <- (which(trajectory == novel.group))
    
    # Find start and end points
    extrema <- as.numeric(names(novel.slice)[c(1,length(novel.slice))])
    
    # Check that start corresponds to novel bin. If False,
    # increase number of supersets to better reflect strutcture.
    if(novel_bin != extrema[1]){
      no.clus <- no.clus +1
    }else{
      
      start <- extrema[1]
      
      # Check that numbers are sequential using custom 
      # function 
      
      if(is.sequential(novel.slice)){
        end <- extrema[2]
      }
      
      # If not sequential, adjust the end year accordingly
      else{
        bool.vec <- diff(novel.slice) == 1
        index <- which(bool.vec == FALSE)[1] - 1
        # correct for 0 index during blips
        end <- ifelse(index > 0,as.numeric(names(bool.vec[index])),novel_bin)
      }
      
      # Calculate Length
      length <- start-end
      
      # Solve for blips, as we want to set the time in years, 
      # not timebins
      if(length == 0){
        length <- 1
      }
      run <- FALSE
    }

  }
  
  # Initialize a class variable
  Class <- "NONE"
  
  # Allocate blips
  if(length == 1){
    
    Class <- "BLIP"
    end <- as.numeric(rownames(matrix.temp[novel_index-1,]))
    length <- novel_bin - end
  }
  
  # Test for full persistence
  if(Class != "BLIP"){
    
    index.vector <- which(names(novel.slice) %in% names(trajectory[novel.slice[1]:length(trajectory)]))
    
    if (length(trajectory[novel.slice[1]:length(trajectory)]) == length(index.vector)){
      
      Class <- "Full Persistence"
    }
  }
  # Rest are short persisters
  if(Class != "BLIP" & Class != "Full Persistence"){
    
    Class <- "Short Persistence"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  
  
  
  # Return data
  return.data <- list(c(list("ID" = ID, "begin" =start, 
                             "end" =end, "length"= length, 
                             "Class" = Class, "clusters" =no.clus)))
  names(return.data) <- ID
  return(return.data)
  
}


# WORK BACK FROM CLUSTER NO. INSTEAD OF COUNTING UP!!!
# That way, we can guarantee signifivance of groups in 
# a nicer fashion. 

# Think about the stopping criterion...

nov.cluster.id.V3 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  number.names <- as.numeric(rownames(matrix))
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  if(nrow(novel_frame) > 1){
    novel_frame <- novel_frame[2,]
  }else{
    novel_frame <- novel_frame[1,]
  }
  # Find year and row data
  novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <<- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance ="czekanowski", 
            method.cluster = "average", alpha=.01), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.01) })
  
  # Initialize parameters for the length calculations
  # Initialize a class variable
  length <- 0
  Class <- "NONE"
  run <- FALSE
  # Plot dendogram result
  par(mfrow = c(1,1))
  
  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <<- TRUE
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(clust.data, hang=-1,
         main = "Blip Dendrogram", xlab = "")
  })
  
  # This will activate the BLIP module and assign the correct 
  # category
  if(run){
    length <- 1
    start <- novel_bin
    end <- novel_bin
    Class <- "BLIP"
  }
  
  
  # Start of the length calculation module
  if(length == 0){
    
    length.output<-length.calculator(test,novel_bin)
  
    length <- length.output$length 
    start <- length.output$start
    end <- length.output$end
  }
  # Initialize a class variable
  
  
  # Allocate blips needs to be done!!!
  if(start == end | length == 1){
    
   Class <- "BLIP"
    
 
  }
  
  # Test for full persistence
  if(Class != "BLIP"){
    
    Class <- "Persister"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  
  
  
  # Return data
  return.data <- list(c(list("ID" = ID, "begin" =start, 
                             "end" =end, "length"= length, 
                             "Class" = Class, "clusters" =no.clus)))
  names(return.data) <- ID
  return(return.data)
  
}



nov.cluster.id.V3(unlisted.mat['palearctic_mat_A.G10145'])




# get structure
length.calculator <- function(test, novel_bin){

  # Extract all dendrogram labels
  names.vec <- unlist(list(test$significantclusters))

  # Construct a new named vector with communities assigned to significant clusters
  new.vec <- NULL
  for(i in 1:length(names.vec)){
    for(j in 1:length(test$significantclusters)){
      if (names.vec[i] %in% test$significantclusters[[j]]){
        new.vec[i] <- j
        names(new.vec)[i] <- names.vec[i]
      }
    }
  }

  # Ensure community labels are sorted in descending order.
  # Use to reorder the cluster vector.
  sorted.index <- as.character(sort(as.numeric(names(new.vec)), decreasing  =T))
  new.vec <- new.vec[sorted.index]

  # Extract the vector index for the novel community
  novel_bin_index <- which(names(new.vec) == novel_bin)


  # Important part. If there are more than 3 significant clusters,
  # cut the Dendrogram into 3 groups as we are just interested in 
  # broad "before", "during" and "after" novelty groups. If there's
  # only 2 significant clusters than we just take the 2.
  cut.length <- ifelse(length(test$significantclusters) >= 3, 3, 
                     ifelse(length(test$significantclusters) == 2, 2, NA))
  #need to work
  if(is.na(cut.length)){
    return(list("length" = 1, "start" = novel_bin,
              "end" = names(new.vec[novel_bin_index+1])))
  }
  # Cut the dendrogram 
  cut.vec <- rev(cutree(test$hclust, k = cut.length))
  
  # Establish "before", "during" and "after" groupings
  group.relations <- vector(mode="list",length=cut.length)

  for(i in 1:length(new.vec)){
  
    group.relations[[as.numeric(cut.vec[i])]] <- c(group.relations[[as.numeric(cut.vec[i])]], new.vec[i] )
  
  }

  # Ensure lsit entries are just unique indices, not double ones
  unique.group.rel <- lapply(group.relations, FUN = unique)

  # Extract the novel group
  group.novel<-new.vec[[novel_bin_index]]

  # Find the group of communities the novel compositin belongs to.
  for(i in 1:length(unique.group.rel)){
    if(group.novel %in% unique.group.rel[[i]]){
      group.no <- i
    }
  }

  # Loop through the vector of clusters to find the start and end of the
  # novel composition
  start <- novel_bin
  end <- NULL
  
  for(i in 1:length(new.vec)){

    # This reads as: If we bump into a non-novel cluster after having
    # already passed the novel community, then the novel state has ended.
    if(new.vec[i] %!in% unique.group.rel[[group.no]] & as.numeric(names(new.vec)[i]) < start){
      end <- as.numeric(names(new.vec)[i-1])
      break
    }
    
    # If we dont meet these requirements, the novel community doesn't end
    if(is.null(end)){
      end <- as.numeric(names(new.vec)[length(new.vec)])
    }
  
  }
  
  # Correct for novel communities at very end
  
  if(novel_bin == as.numeric(names(new.vec[length(new.vec)]))){
    end <- novel_bin
  }
  
  
  # Persistence length is simply the start of the novel state until we 
  # enter something significantly different.
  length <- start-end

  return(list("length" = length,
              "start" = start,
              "end" = end))
}

# These new groups can simply be used to calculate length!
length.calculator(test, novel_bin)



identify.novel.gam(unlisted.mat[4])

blip <- 0
pers <- 0
end <- 0
for(i in (novelty.pers.V3)){
  
  if(i$Class == "BLIP"){
    blip = blip +1
    
  }
  if(i$Class == "Persister"){
    pers = pers +1
    
  }
  if(i$Class == "END"){
    end = end +1
    
  }
}
blip/503
pers/503
end



####################################
### Step X. Sensitivity Analyses ###
####################################

# We need to conduct a number of sensitivity analyses to improve validity of work done.

# 1. Varying alpha for detection of novelty

# 2. Random selection of species group to test if there is a significant association (invader models)

# This means we randomly assign the invader tag to the species in each country. We must take care
# to have the number of "fake" invaders in each country to match the real invaders.

# 3. Testing the ANOSIM method on non-novel time series; I believe this shows that the novelty framework
# is actually quite conservative. Big changes are found throughout the time series, but a big change is not
# always the same as a novel community. It seems the method is quite robust indeed. 

ANOSIM.supp.test <- function(number){
  
  id.vector <- NULL

  while (length(id.vector) < number){
  
  # subset all analysed timeseries
  unique.timeseries <- unique(full.novel.mat[, c("site", "cat")])
  
  # create a random sample of size 10
  rand.index <- sample(seq(1, nrow(unique.timeseries)), size = 1)
  
  # isolate the data for each randomly selected timeseries
  
  ID <- unique.timeseries$site[rand.index]
  
  matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
  
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[ID]]
    
  }
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = FALSE, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)

  # check if it does not contain novelty types (necessary for this analysis)
  if (any(label_frame$cat %in% c("novel", "instant", "cumul"))){
    
    id.vector <- c(id.vector, NULL)
    
  } else id.vector <- c(id.vector, ID) # only add if there is no novelty
  
}

  # Next part, run ANOSIM algorithm on these time series. But, we need to randomly assign a 
  # 'novel community' so the ANOSIM has groups to work with.

  sim.trajectories <- ANOSIM.sensitivity.lists(id.vector)

  lengths <- ANOSIM.checker(sim.trajectories, cut.off = 0.9)
  

    

  
  return((lengths.fil))

}

shouldbezero <- ANOSIM.supp.test(10) 


# 4. Does decreasing the number of communities within each ANOSIM group automatically lead to a 
# higher R statistic? Possibly, but it will also reduced the significance of the comparison.








#######################################################################
##### Step X. Publication Figures #####################################
#######################################################################

##########
# Fig 1.##
##########

# It will be important to elegantly visualize the rate of novelty community emergence
# within Freshwater Fish communities globally.

# Plot novelty emergence over time using GAM's, input is any data frame with a bins and novelty column

emergence.plotter <- function(full.novel.mat, category){
  par(mfrow=c(1,1))
  n <- as.data.frame.matrix(table(full.novel.mat$cat, full.novel.mat$bins))
  colnames(n) <- 2020 - as.numeric(colnames(n))

  

  for (i in 1:ncol(n)){
  
    n["total", i] <- sum(n[1:4, i])
    n["per.novelty", i] <- round(n[category, i]/n["total", i], digits = 3)
    n["bins",i] <- colnames(n[i])
  }

  n <- as.data.frame(t(n))
  rownames(n) <- NULL

  n$cex <- as.numeric(n$total)/max(as.numeric(n$total))

  v1 <- as.numeric(n$per.novelty)
  v2 <- as.numeric(n$bins)

  a <- as.data.frame(cbind(v1,v2))
  print("Running GAM for plot")
  test_gam <- gam(v1 ~ s(v2, bs = "tp", k = -1), weights = as.numeric(n$total), data = n,
      family = betar(link = "logit"), 
      method = "REML")

  plot(v1~v2, type="n",
       ylim=c(0, 0.1),
       xlim=c((min(v2)+5),(max(v2)-1)),
       axes=T, ylab = "Proportion of Novelty", xlab = "Time (Years)")

  # Confidence intervals
  
  se <- predict(test_gam, se = T)$se.fit
  fit <- plogis(predict(test_gam , se = TRUE )$fit)
  
  lcl <- fit - 1.96 * se
  ucl <- fit + 1.96 * se
  
  polygon(x=c(v2, rev(v2)),
          y=c(lcl, rev(ucl)), 
          col=ifelse(category=="novel","gold",
                         ifelse(category=="instant","coral1", "darkslategray1")), border="black", lty = "dashed")
  
  # Add border
  lines(plogis(predict(test_gam)) ~ v2, col="black", lwd = 5)
  lines(plogis(predict(test_gam)) ~ v2, col=ifelse(category=="novel","orange",
                                                   ifelse(category=="instant","red", "darkturquoise")), lwd=3)
 
  points(v1~v2, bg = ifelse(category=="novel","orange",
                            ifelse(category=="instant","red", "lightblue")), 
         pch = 21,
         cex = (n$cex*1.5),col = "black")
  
  
return("Novelty emergence through time")
}

# Save to desktop
pdf(file = "/Users/sassen/Desktop/plot_1.pdf",
    width = 10,
    height = 5)

emergence.plotter(full.novel.mat, "novel")
dev.off()

##########
# Fig 2. #
##########

# A fancy histogram or something like that which summarizes the persistence lengths.
# It would be cool to have a graph that has time on the x-axis and some form of density
# to show when novelty was occurring.

persistence.plotter <- function(nov.summary){

  nov.time <- subset(nov.summary, state == "novel")

  # find the durations of the novel communities
  start.vector <- 2019 - nov.time$begin
  end.vector <- 2019 - nov.time$end
  x = c(1960:2019)
  nov.duration <- NULL

  timeframe <- min(start.vector):max(end.vector)

  counts <- as.list(rep(0,length(x)))
  names(counts) <- x

  for (i in 1:length(start.vector)){
  print(i)
  nov.duration <- start.vector[i]:end.vector[i]
  
  for (j in 1:length(counts)){
    
    if(names(counts)[j] %in% nov.duration){
      
      counts[j] = counts[[j]] + 1
      
    }
  }
  
  
}

  # Find the timeseries we actually used
  
  used.sites <- unique(full.novel.mat$site)
  unique.combos <- unique(Survey_Data[, c("Year", "TimeSeriesID")])
  non.nov.time <- unique.combos[unique.combos$TimeSeriesID %in% used.sites,]
  
  
  #non.nov.time <- unique(Survey_Data[, c("Year", "TimeSeriesID")])

  start.non.nov <- min(non.nov.time$Year)
  end.non.nov <- max(non.nov.time$Year)

  timeframe.non <- min(start.non.nov ):max(end.non.nov)
  counts.non <- as.list(rep(0,length(timeframe.non)))
  names(counts.non) <- timeframe.non

  for (i in 1:length(counts.non)){
  print(i)
  
  year_sub <- subset(non.nov.time, Year == as.numeric(names(counts.non)[i]))
  
  counts.non[[i]] <- nrow(year_sub)
  
}

  whole.val <- counts.non[names(counts.non) %in% as.character(c(1960:2019))]

  y <- (as.numeric(counts)/as.numeric(whole.val))
  y2 <- (as.numeric(whole.val))

  par(mar = c(5.1, 4.1, 4.1, 5.1))
  barplot(height = y2, ylab = "", xlab = "", yaxs="i", axes = F, col = "steelblue")
  axis(side = 4, at = pretty(range(0, 2000)))
  
  mtext(side = 4, "Number of Time Series", line = 3)
  
  par(new = TRUE)
  
  
  
  plot(x = x, y = y, type = "l", xlim = c(1960,2019), ylim = c(0,0.06),
       ylab = "Proportion of sampled communities ", xlab = "Time (years)", bty = "U"
       , yaxs="i", axes = F)
  polygon(c(2019, 1960, x), c(0,0, y), col = "orange", border = "gold", yaxs="i", lwd = 1.5)
  axis(side = 1, at = pretty(range(1960,2019)))
  axis(side = 2, at = pretty(range(0.00,0.06)))
  
  
  gc()
  
  
}

# Save to desktop
pdf(file = "/Users/sassen/Desktop/plot_2.pdf",
    width = 10,
    height = 5)

persistence.plotter(nov.summary)
dev.off()


##########
# Fig 3. #
##########

# We will need to include a map of the world with our time series locations and probably also 
# HydroBasins. 


create.map <- function(){
# Make a map

  # get map
  worldmap <- getMap(resolution = "coarse")
  # plot world map
  plot(worldmap, col = "black", bg = "lightblue",
       xlim =c(-160, 160), ylim = c(-50, 100),
       asp = 1)

  novels <- subset(time_series_data, TimeSeriesID %in% full.novel.mat$site)[,c("Latitude","Longitude")]
  points(novels$Longitude, novels$Latitude,col = "black", cex = 0.6, pch = 21, bg = "orange")
}

pdf(file = "/Users/sassen/Desktop/plot_3.pdf",
    width = 15,
    height = 10)

create.map()

dev.off()

##########
# Fig 4. #
##########

# We will include an example nMDS (k = 2) plot of our novel community transitions paired with
# the classic NDF plots from (Pandolfi, Staples and Kiessling). Golden polygon surrounding the 
# communities within a novel state.



##########
# Fig 5. #
##########

# A figure that somehow presents the results of the invader - persistence investigation is necessary.
#



require(devtools)
install_version("clustsig", version = "1.1", repos = "http://cran.us.r-project.org")

mod<-glm(novel~BioRealm, data = full.novel.mat.season, family = 'binomial')
summary(mod)



types <- NA
for (site in full.novel.mat.season$site_ID){
  index<- which(Survey_Data$TimeSeriesID == site)
  types <- c(unique(Survey_Data$UnitAbundance[index]), types)
}
  
unique(types)



