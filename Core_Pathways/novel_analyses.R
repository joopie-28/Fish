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

# We now have a method of 1) visualizing differences between communities in two-dimensional
# space and 2) determining whether the difference between pre- and post-novel states is greater
# than the difference within these groups.

# Now, we want to use these techniques to discover how long a novel state is maintained (in years).
# Because our ANOSIM test simply computes differences between groups, there might be some communities
# that are in fact quite different from the novel community (i.e. on their way back to a pre-novel state)
# in the post-novel group.

# Isolate novelty matrices
matrices <- lapply(Fish_Communities_A$BioRealm_Matrices_A[], 
                   function(x){
 
   indices <- names(x) %in% subset(full.novel.mat, cat == "novel")$site
  
   return(x[indices])
})

# Execute ANOSIM method 
novelty.lengths <- do.call(c, 
                           lapply(1:5, function(x){
                             
  novelty.lengths <- novel.length.algo(matrices[[x]])
  
}))

# Add the calculated lengths to our main data frame
full.novel.mat$novel.length <- NA
full.novel.mat$novel.class <- NA

for(i in 1:length(novelty.lengths)){
  indices <- which(full.novel.mat$site == novelty.lengths[[i]]$Simulation & full.novel.mat$cat == "novel")[1]
  full.novel.mat$novel.length[indices] <- novelty.lengths[[i]]$Length_years
  full.novel.mat$novel.class[indices] <-  novelty.lengths[[i]]$Class
  
}

# Set non-novel lengths to 0 instead of NA
full.novel.mat$novel.length[is.na(full.novel.mat$novel.length)] <- 0
full.novel.mat$novel.class[is.na(full.novel.mat$novel.class)] <- "NONE"



################################################
#### Step 6. Analyzing state-level metrics ####
###############################################

# Determination of lengths allows us to look at characteristics of the novel state as it evolves.
# We have an unbiased method that tells us when the state ends, so we know what it "novel" and
# what is not.

# This functions extract information from the novelty results and returns a data frame with summary
# statistics. It also returns a data frame with 'state-level metrics'; a summary of important 
# biological metrics aggregated at distinct periods within each time frame.

# NEED TO UPDATE THIS UNFINISHED PER 16-06-2022

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

nov.cluster.id <- function(matrix){

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

# Check if overarching groups are significant (PART OF nov.cluster.id)

# Need to make some iterations which test that groups are composed of significant subgroups
# rather than mixed.

test$significantclusters

nov.cluster.id(matrices[[1]][13])
# Execute persistence clustering over all novel time series
# Need to correct for inability to cluster some..

unlisted.mat <- unlist(matrices, recursive = F)

novelty.pers <- do.call(c, 
                           lapply(1:length(unlisted.mat), function(x){
                             print(x)
                             nov.cluster.id(unlisted.mat[x])
                             
                           }))

# Add the calculated lengths to our main data frame
full.novel.mat$novel.length <- NA
full.novel.mat$novel.class <- NA

for(i in 1:length(novelty.pers)){
  indices <- which(full.novel.mat$site == novelty.pers[[i]]$ID & full.novel.mat$cat == "novel")[1]
  full.novel.mat$novel.length[indices] <- novelty.pers[[i]]$length # need to fix end point comm.
  full.novel.mat$novel.class[indices] <-  novelty.pers[[i]]$Class
  
}








# Plot NMDS next to Dendrogram

mds.cluster.plotter <- function(matrix){
  
  # Set plotting params
  par(mfrow = c(1,2))
  
  
  
  
  # Extract labels
  label <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                  alpha = 0.05,
                                  metric = "bray",
                                  plot = F, 
                                  site = "NA",
                                  plot.data = FALSE,
                                  gam.max.k = -1)
  
  matrix$category <- label$cat
  matrix$colour <- NA
  
  # Assign colours
  for (i in 1:nrow(matrix)) {
    if(matrix$category[i] == "cumul"){
      matrix$colour[i] <- "skyblue"
    }
    if(matrix$category[i] == "instant"){
      matrix$colour[i] <- "red1"
    }
    if(matrix$category[i] == "novel"){
      matrix$colour[i] <- "orange"
    }
    if(matrix$category[i] == "back"){
      (matrix$colour[i] <- "grey")}
  }
  
  # Run MDS 
  NMDS=metaMDS(matrix[,-c(ncol(matrix), (ncol(matrix)-1))], 
               k=2, trymax = 10000)
  plot(NMDS, type = "n", ylab = "MDS2", xlab = "MDS1")
  points(x = NMDS$points[,"MDS1"], 
         y = NMDS$points[, "MDS2"], 
         bg = matrix$colour,
         pch = 21,
         col = "black",
         cex = 1.4)
  
  for (i in 1:(nrow(NMDS$points)-1)){
    
    arrows(x0 = NMDS$points[i,"MDS1"], 
           y0 = NMDS$points[i, "MDS2"], 
           x1 = NMDS$points[i+1,"MDS1"], 
           y1 = NMDS$points[i+1, "MDS2"], length = 0.05, lwd = 1)
  }
  
  print("Clustering")
  
  test <- simprof(data = matrix[,-c(ncol(matrix), (ncol(matrix)-1))], num.expected = 1000,
                  num.simulated = 999, method.distance ="braycurtis", 
                  method.cluster = "centroid"
                  ,alpha=0.05, undef.zero = T)
  
  temp <- simprof.plot(test, plot = F)
  
  dendro.col.df <- data.frame(labels = as.numeric(labels(temp)))
  dendro.col.df$colour <- NA
  dendro.col.df$cat <- NA

  for(i in 1:nrow(dendro.col.df)){
    
    index <- which(as.numeric(label$bins) == dendro.col.df$labels[i])
    dendro.col.df$cat[i] <- label$cat[index]
    
    if(dendro.col.df$cat[i] == "cumul"){
      dendro.col.df$colour[i] <- "skyblue"
    }
    if(dendro.col.df$cat[i] == "instant"){
      dendro.col.df$colour[i] <- "red1"
    }
    if(dendro.col.df$cat[i] == "novel"){
      dendro.col.df$colour[i] <- "orange"
    }
    if(dendro.col.df$cat[i] == "back"){
      dendro.col.df$colour[i] <- "grey"}
  }
  
  # Plot the Dendrogram
 
  labels_colors(temp) <- dendro.col.df$colour
  
  labels(temp) <- paste0(dendro.col.df$cat, "-", dendro.col.df$labels)
  
  plot(temp, ylab = "Height")
  
}

pdf(file = "/Plots/clustering_persistence.pdf",
    width = 12,
    height = 6)

mds.cluster.plotter(matrices[[1]][[295]])
dev.off()


matrices[[1]][295]
blip <- 0
full <- 0
short <- 0
for(i in 1:length(novelty.pers)){
  
  if(novelty.pers[[i]]$Class == "BLIP"){
    blip = blip + 1
    
  }
  if(novelty.pers[[i]]$Class == "Short Persistence"){
    short = short + 1
    
  }
  if(novelty.pers[[i]]$Class == "Full Persistence"){
    full = full + 1
    
  }
  
}

full/503
short/503
blip/503
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








