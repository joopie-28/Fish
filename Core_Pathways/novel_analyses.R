##### This script includes all pre-processing steps and major analyses #####

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
                 "gridExtra"))

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

###############################################################
#### Step 2. Computing invasive turnover metrics  #############
###############################################################

# Computes certain ecological metrics of invaders and natives

full.nnc.matrices <- mat.nnc.ass(full.stat.matrices)

###############################################################
#### Step 3. Building a frame containing community metrics ####
###############################################################

# This function creates a all-encompassing data frame f
# suitable for use in modelling 

full.novel.mat <- inv.frame.builder(full.nnc.matrices)

###########################################################################
#### Step 4. Modelling emergence of novelty by invasives turnover #########
###########################################################################

# This function runs GLMM's for the invasives-novelty investigations. If plot = T,
# Additional plots of the invader effect on novelty probability are saved. 
# CSV's with model output are simultaneously created

inv.models <- inv.nov.glmm(full.novel.mat, plot = T)

#################################################################
#### Step 5. Investigating Novelty Persistence  #################
#################################################################

# Here, we use ANOSIM to further analyse what happens to a novel state
# after it emerges.

# This function applies the ANOSIM algorithm to all novel comm.
# and returns a list of summary statistics which help us to
# objectively classify 'blips'. It takes the novelty framework
# output as an input. 

anosim.plots.list <- anosim.analyzer(Fish_Communities_A)

# We now have a method of 1) visualizing differences between communities in two-dimensional
# space and 2) determining whether the difference between pre- and post-novel states is greater
# than the difference within these groups.

# Now, we want to use these techniques to discover how long a novel state is maintained (in years).
# Because our ANOSIM test simply computes differences between groups, there might be some communities
# that are in fact quite different from the novel community (i.e. on their way back to a pre-novel state)
# in the post-novel group.

# I propose that we remove one community (starting with the youngest) and examine whether or not this 
# leads to an R statistic closer to 1. If it does, that means that the removed community was actually closer
# to the pre-novel states than the original R statistic made it appear.

# We will compute the length of the novel state, assuming that a dissimilarity of 0.85 is sufficient.
# Essentially, when a dissimilarity of >0.85 is encountered, the length of that period is validated.
# A group of novel states may span quite some time, but if the R statistic is low, these 'novel states'
# are simply extensions of the pre-novel states and thus the novel community is assigned a duration of 1 
# (which is just the actual novel event).

# This function executes everything discussed above
full.novelty.lengths <- novel.length.calculator(anosim.plots.list, cut.off = 0.85)

# Add the calculated lengths to our main data frame
full.novel.mat$novel.length <- NA

for(i in 1:nrow(full.novel.mat)){
  print(i)
  
  for (j in 1:length(full.novelty.lengths)){
    
    if(full.novel.mat$cat[i] == "novel" & full.novel.mat$site[i] == full.novelty.lengths[[j]]$ID){
      
      full.novel.mat$novel.length[i] <- full.novelty.lengths[[j]]$length
      
    }
  }
}

full.novel.mat$novel.length[is.na(full.novel.mat$novel.length)] <- 0


################################################
#### Step 6. Analyzing state-level metrics ####
###############################################

# Determination of lengths allows us to look at characteristics of the novel state as it evolves.
# We have an unbiased method that tells us when the state ends, so we know what it "novel" and
# what is not.

# This functions extract information from the novelty results and returns a data frame with summary
# statistics. It also returns a data frame with 'state-level metrics'; a summary of important 
# biological metrics aggregated at distinct periods within each time frame.

nov.comm.summary <- novel.comm.analyzer(full.novel.mat)

# Save summary stats to a CSV file

write.csv(nov.comm.summary$summary, 
          "./outputs/nov_summary_stats.csv")


#### End of Main Analyses ####




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



