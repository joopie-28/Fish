# Identifiying Novel Communities (as per the Method of Pandolfi et al., 2020)

# Clear environment (if wanted) and set your working directory
rm(list = ls())
setwd("YOUR OWN CHOICE HERE")

# Packages
install.packages(c("mgcv", "vegan", "lme4", "nlme", 
                 "DHARMa", "merTools", "shape",
                 "multcomp", "maptools", "sp", 
                 "divDyn", "plotrix", "raster",
                 "rgeos", "fun", "analogue",
                 "brms"))
# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")

# Here, we use the two-metric approach to identify a novel community, using a freshwater fish
# time-series as a test. Much of the code is directly taken from the 2020 publication by
# Pandolfi et al. (Increased extinction in the emergence of novel ecological communities). We
# have repurposed it to (hopefully) identify novel fish communities on a scale of decades rather than 
# millenia. 

# Step 1: Transforming our data into a species-specific presence abundance matrix, including columns denoting
# bin timeframe (we will work with 1 year steps in this demo) and site location (kind of pointless in this "one site"
# demo, but necessary if we expand to more timeseries). The function used for this transformation was made by me (Joop).

# I have made a quick alteration here which migth need to be addressed if we keep using this method; some years have multiple surveys
# at different points within a year. As per Compte et al (2021), I have decided to only include surveys that were taken at the same point
# in the year (e.g. quarter 3). If there was more than one survey at this time, I randomly dropped one to avoid any bias.

binary_converter_function <- function(survey_identifier){
  
  # Here I create some reference lists 
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  
  # Here I create the empty dataframe which will be populated with presence/absence values
  binary_df <- as.data.frame(matrix(data = NA, nrow = nrow(unique(surveys_unique)), ncol = nrow(ref_list)))
  colnames(binary_df) <- ref_list[, 1]
  rownames(binary_df) <- surveys_unique[, 1]
  
  # Let's populate this data frame. I will use the column and row names to loop through a subset of the data
  # a couple times. This will essentially crosscheck if a species is present in the set or not. 
  for (i in 1:nrow(binary_df)) {
    sub_temp <- subset(sub_series, SurveyID == rownames(binary_df)[i])
    for (j in 1:ncol(binary_df)) {
      if (colnames(binary_df[j]) %in% sub_temp$Species){
        binary_df[i,j] <- 1
      }
      else{
        binary_df[i,j] <- 0
      }
    }
  }
  
  # Making sure we get time values in the dataframe
  unique_df <- subset(unique(sub_series[,c("SurveyID","Year", "Quarter")]))
  unique_df$real_time <- 0
  for (i in 1:nrow(unique_df)) {
    # This ensures no NA produced due to silly things like Q2/3 
    if (typeof(unique_df[i,3]) != "character"){
      unique_df[i, 4] <- (as.numeric(unique_df[i,2]) + as.numeric(unique_df[i,3])*0.225)
    }
    else {
      unique_df[i, 4] <- as.numeric(unique_df[i,2])
    }
  }
  
  binary_df$time <- 0
  binary_df$quarter <- 0
  for (i in 1:nrow(binary_df)) {
    for (j in 1:nrow(unique_df)) {
      if (rownames(binary_df[i,]) == unique_df[j,1]) {
        binary_df$time[i] <- unique_df[j,4]
        binary_df$quarter[i] <- unique_df[j, 3]
      }
    }
  }
  
  # Selecting only those that occur in the MOST CONSISTENT QUARTER, then deleting the quarter column. But only if we have enough
  # Timepoints to justify this (no point deleting them if we only have 2 or 3)
  if (nrow(binary_df) >= 5) {
  binary_df <- subset(binary_df, quarter == as.numeric(tail(names(sort(table(binary_df$quarter))), 1)))
  }
  binary_df <- binary_df[, !(names(binary_df) == "quarter")]
  
  # Arbitrarily removing duplicate surveys based on time frame. 11-12-21 I HAVE A BETTER IDEA FOR THIS. SELECT THE ONE WITH 
  # THE HIGHEST NUMBER OF SPECIES RATHER THAN JUST DOING IT RANDOMLY. THIS MAKES IMPROVES VALIDITY OF NOVELTY DETECTION.
  binary_df <- binary_df[!duplicated(binary_df[,c("time")]),]
  
  # Change row names to bin (time)
  row.names(binary_df) <- (2021-as.numeric(binary_df$time))
  
  # Now remove the time column
  binary_df <- binary_df[, !(names(binary_df) == "time")]
  
  # Very small samples will be returned as a string of numbers (type = double) rather than data frames, so will instruct the function to skip this survey if
  # That happens. The timeseries is just too short basically.
  if (typeof(binary_df) == "double"){
    return("Not enough data")
  }
  
  # Make sure the order is OK, it will influence the model results (our priamry reference community HAS to be the oldest one for the
  # model to make sense).
  binary_df <- binary_df[order(as.numeric(row.names(binary_df)), decreasing = FALSE),]
  
  return(binary_df)
}



