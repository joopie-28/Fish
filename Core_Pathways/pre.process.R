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


##################################################
#### Step 1. Tagging species status  #############
##################################################

# Here, we will use a combination of the Tedesco et Al.
# database and our own custom species tags. 

# First we create a list of TimeSeries for all the countries we're 
# interested in.

country_list_assigner <- function(country){

  country.lower <- tolower(country)

  a <- as.list(subset(time_series_data, Country == country)$TimeSeriesID)
  names(a) <- rep(subset(time_series_data, Country == country)$Country[1], times = length(a))

  return(a)
}
  
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

# Using alternate method (no-non-natives, only invaders and natives)
# Species will be tagged and a list of matrices is returned.

full.stat.matrices <- assign.stat.country_nn(full.ID.list) 

###############################################################
#### Step 2. Computing invasive turnover metrics  #############
###############################################################

# Compute contribution of natives and invasives 

full.nnc.matrices <- mat.nnc.ass(full.stat.matrices)

###############################################################
#### Step 3. Building a frame containing community metrics ####
###############################################################

# This function creates a data frame for modelling 

full.novel.mat <- inv.frame.builder(full.nnc.matrices)

###########################################################################
#### Step 4. Modelling emergence of novelty by invasives turnover #########
###########################################################################

model.plotter.6 <- function(invader.models, points){
  
  if (points){   
    plot_model(invader.models$model$novel, type = "pred",  terms=c("INC_increase [all]"), 
               title = "Probability of Novelty Emergence explained by Invader dynamics",
               axis.title = c("Net change in Invader Relative Abundance","Probability of Novel State (%)"),
               pred.type = "fe", colors = "green", show.data = T)
  }else{
    
    
    # Plot effects with ggplot
    
    df_novel <- ggpredict(invader.models$novel$model, type = "fe",  terms="INC_increase [all]")
    df_instant <- ggpredict(invader.models$instant$model, type = "fe" ,terms= "INC_increase [all]")
    df_cumul <- ggpredict(invader.models$cumul$model, type = "fe" ,terms= "INC_increase [all]")
    
    # Plot invader effect holding all else equal
    
    p1 <- ggplot(df_novel, aes(x, predicted)) +
      geom_line(color = "gold", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "yellow") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "True Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    p2 <- ggplot(df_instant, aes(x, predicted)) +
      geom_line(color = "red1", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "red3") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "Instantaneous Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    p3 <- ggplot(df_cumul, aes(x, predicted)) +
      geom_line(color = "steelblue4", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "mediumturquoise") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "Cumulative Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    # 3-way plot
    grid.arrange(p1,p2,p3, nrow = 3)
  }
}

inv.nov.glmm <- function(full.novel.mat, plot){

  cat.var <- c("instant", "cumul", "novel")
  
  # This function runs a GLM for each novelty category
  invader.models <- lapply(cat.var, function(cat.1){
    
    # Had to do it this way, bit clunky but works fine
    full.novel.mat$response <-  full.novel.mat[,..cat.1]  

    print(paste0("Running ", cat.1, " invader model"))  
    
    # GLMM with interaction term and a nested random effect
    mod <- glmer(response ~ bin_lag + position + INC_increase + (1|country/basin),
                 data = full.novel.mat, family = binomial)
    
    full.novel.mat$response <- NULL
    
    # Model Diagnostics with DHARMA
    disp.test <- testDispersion(simulateResiduals(mod))
    print(ifelse(disp.test$p.value <= 0.05, 
                 paste0("Dispersal not OK ", round(disp.test$statistic, 2)),
                 "Dispersal OK"))
  
    # Add coefficients to output to make inspection easier
    mod.1 <- list("model" = mod, "summary" = as.data.frame(summary(mod)$coef))
    
    return(mod.1)
  
  })
  
  names(invader.models) <- cat.var
  
  if(plot){
    
    # if true, model plots for the invader effect are returned.
    model.plotter.6(invader.models, points = F)
  }
  
  return(invader.models)
}

# This function runs GLMM's for the invasives-novelty investigations. If plot = T,
# Additional plots of the invader effect on novelty probability are plotted. 

inv.models <- inv.nov.glmm(full.novel.mat, plot = T)


#################################################################
#### Step 5. Investigating Novelty Persistence  #################
#################################################################

# Generate a list of Time Series where novelty was discovered

rel.abu.1.com <- rbindlist(lapply(list(Fish_Communities_A$BioRealm_Novelty_A$palearctic_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$nearctic_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$afrotropics_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$neotropics_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$australasia_novelty_A), FUN = "rbindlist"))

# List of sites in our community list
novel.list <- as.list(rel.abu.1.com[rel.abu.1.com$cat == "novel"]$site)

# We will use ANOSIM to compare pre and post novel groups. Once we have 
# established which
# communities are not ecological blips, we can investigate how long novelty
# persists.

# Applies novelty.trajectory.plotter to a list of Time Series
anosim.plots.list <- novelty.trajectory.lists(novel.list)

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


# Compute the length of the novel state period, assuming that a dissimilarity of 0.85 is sufficient.
# Essentially, when a dissimilarity of 0.85 is encountered, the length of that period is validated.
# A group of novel states may span quite some time, but if the R statistic is low, these 'novel states'
# are simply extensions of the pre-novel states and thus the novel states are assigned a duration of 1 
# (which is just the actual novel event).


# Filter out novel communities that show a lack of difference. It also 
# returns the data in a useful format. 

true.novelty.output <- cut.off.generator(anosim.plots.list)

# Calculate lengths of novel communities, if they exceed the threshold value. Returns length and R stat. 

full.novelty.lengths <- novel.length.checker(true.novelty.output, cut.off = 0.85)

# Add on to big df straight away!
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

## Create new data frame for analysis of the novel communities.##
## look at average ecological metrics for the novel community instead of just the transition. ##


# Create a frame with only novel communities for further analysis
nov.comm.all <- subset(full.novel.mat, cat == "novel")

# Isolate non-blips
nov.comm.no.blip <- subset(full.novel.mat, cat == "novel" & novel.length > 0)

# New state-level df for analyses of the non-blip novel communities through time
nov.summary <- rbindlist(lapply(nov.comm.no.blip$site, function(ID){
  
  # Identify the novel points using the length we calculated
  print(ID)
  id <- ID # start with a single to see how it works
  
  # Extract length of the novel state
  state.length <- as.double(unique(nov.comm.no.blip[which(nov.comm.no.blip$site == id), "novel.length"]))
  
  # Use that length to extract all the points that make up the novel state
  
  year.start <- as.double(min(nov.comm.no.blip[which(nov.comm.no.blip$site == id), "bins"]))
  
  year.end <- year.start - state.length # Easy way to find the endpoints
  
  # Add a little variable that tracks whether or not the state is "still counting"
  # We compare the nov. comm. end bin to the timeseries end bin.
  
  pres.stat <- ifelse(year.end == min(full.novel.mat[full.novel.mat$site == id, "bins"]), "Still counting", "Over")
  
  
  # Now we extract the points from our df.
  
  years.nov <- c(year.start:year.end) # this will include values with no communities, but that's fine
  
  nov.stat.comm <- subset(full.novel.mat, site == id & bins %in% years.nov) # There are the communities within a novel state!
  
  # There is always a transition period. Sometimes this is a bit longer than one year,
  # I think it would be fair to add this to the novel community, just to count ext and orig.
  temp <- subset(full.novel.mat, site == id)
  trans.gap <- temp$bins[which(temp$cat == "novel")-1] - temp$bins[which(temp$cat == "novel")] - 1
  transition.period <- ifelse(is_empty(trans.gap), 0, trans.gap)
  
  # compute average diversity metrics, we want PER YEAR, NOT PER TIMEJUMP for extinctions and originations
  
  mean.evenness.n <- mean(nov.stat.comm$evenness) # normal mean is ok here
  mean.shannon.n <- mean(nov.stat.comm$shannon.d)
  mean.ext.n <- sum(nov.stat.comm$ext, na.rm = T)/(state.length +transition.period) # more accurate this way
  mean.orig.n <- sum(nov.stat.comm$orig, na.rm = T)/(state.length +transition.period)
  
  # I think we also need to extract the non-novel communities to do comparisons
  
  # Save total years non-novel in case we want it as a covariate
  
  years.non.nov <- max(temp$bins) - min(temp$bins) - state.length - transition.period
  
  # Non-novel comms
  back.stat.comm <- subset(full.novel.mat, site == id & bins %!in% years.nov)
  
  mean.evenness.b <- mean(back.stat.comm$evenness)
  mean.shannon.b <- mean(back.stat.comm$shannon.d)
  mean.ext.b <- sum(back.stat.comm$ext, na.rm = T)/years.non.nov
  mean.orig.b <- sum(back.stat.comm$orig, na.rm = T)/years.non.nov
  
  # Invaders
  mean.invaders.n <- mean(nov.stat.comm$INC)
  mean.invaders.b <- mean(back.stat.comm$INC)
  
  # Return a new df with state-level metrics
  
  state.level.metrics <- data.frame(matrix(data = NA, nrow =2))
  
  state.level.metrics$state <- NA
  state.level.metrics$length <- NA
  state.level.metrics$mean.even <- NA
  state.level.metrics$mean.shannon <- NA
  state.level.metrics$mean.orig <- NA
  state.level.metrics$mean.ext <- NA
  state.level.metrics$site <- id
  state.level.metrics[,1] <- NULL # remove the auto column 
  state.level.metrics$country <- unique(temp$country)
  state.level.metrics$invaders <- NA
  state.level.metrics$begin <- year.start
  state.level.metrics$end <- year.end
  state.level.metrics$status <- pres.stat
  
  # There's one more fun thing we can add, the ratio between non-novel and novel for orig and ext.
  # It will just be entered for the novel community in the df
  
  # We can now fill in the new df based on some distinctions between states,
  # this looks a bit tedious but is actually the quickest way.
  
  state.level.metrics[1, c("state", "length", "mean.even", 
                           "mean.shannon", "mean.orig", "mean.ext", "invaders")] <- c("novel", state.length, mean.evenness.n,
                                                                                      mean.shannon.n, mean.orig.n, mean.ext.n, mean.invaders.n)
  state.level.metrics[2, c("state", "length", "mean.even", 
                           "mean.shannon", "mean.orig", "mean.ext", "invaders")] <- c("back", years.non.nov, mean.evenness.b,
                                                                                      mean.shannon.b, mean.orig.b, mean.ext.b, mean.invaders.b)
  
  
  
  return(state.level.metrics)
  
}))

# Add a status variable for persistence type AND for if the state ended before the time series
nov.comm.all$type <- NA
nov.comm.all$status <- NA

# Loops to fill in status and type data
for (i in 1:nrow(nov.comm.all)){
  if (nov.comm.all$site[i] %in% nov.comm.no.blip$site){
    nov.comm.all$type[i] <- "persistent"
  } else {nov.comm.all$type[i] <- "blip"}
}
for (i in 1:nrow(nov.comm.all)){
  print(i)
  for (j in 1:nrow(nov.summary)){
    if(nov.comm.all$site[i] == nov.summary$site[j]){
      nov.comm.all$status[i] <- nov.summary$status[j]
    }
  }
}

# We now have a complete data frame with all the information needed to do some evaluations

# How many are blips?
perc.blip <- nrow(subset(nov.comm.all, type == "blip"))/(nrow(subset(nov.comm.all, type == "persistent")) + nrow(subset(nov.comm.all, type == "blip")))
perc.no.blip <- 1-perc.blip

# How many persist till the end of the time series?
no.end.comm <- nrow(subset(nov.comm.all, status == "Still counting"))/(nrow(subset(nov.comm.all, status == "Still counting")) + nrow(subset(nov.comm.all, status == "Over")))





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



