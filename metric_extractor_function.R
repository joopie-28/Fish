# The ultimate metric-maker function

# Clear environment (if wanted) and set your working directory
rm(list = ls())
setwd("YOUR OWN CHOICE HERE")

# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")

# I've made this function special, it's a bit weird, but the input allows you to decide what you want from the function.
# the "survey_identifier" input is the name of the survey (e.g. "G1034"), the other inputs are binary variables. Pick the one
# you want to extract and input 1 for that variable and 0 for all others (note, if you put 1 for more than one decider variable,
# an error will occur). The output varies from single numbers, to strings, to complete dataframes.
metric_extractor <- function(survey_identifier, subset_dataframe, number_of_surveys, timeseries_timepoints, time_range, species_in_survey, useful_df) {
  
  # Create the subset we're interested in, aka the timeseries we're interested in.
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  
  # How many surveys does the timeseries consist of?
  unique_surveys <- length(unique(sub_series$SurveyID))  
  
  # How many distinct timepoints does the timeseries consist of (note, some surveys taken at same time)
  sub_series$real_time <- sub_series$Year + (as.numeric(sub_series$Quarter)*0.225)
  order(sub_series$real_time)
  sub_series[order(sub_series$real_time),]
  unique_timepoints <- length(unique(sub_series$real_time))  
  
  # What is the actual time range?
  a <- min(as.numeric(unique(sub_series$real_time)))
  b <- max(as.numeric(unique(sub_series$real_time)))
  output <- paste("For survey",survey_identifier,"From",a,"-",b)
  
  # How many species does the timeseries include, not correcting for surveys but looking at the timeseries as a whole.
  unique_species <- length(unique(sub_series$Species))  
  
  # Zoom in and look at the species per survey, or better yet, unique time point.
  df_new <- unique(sub_series[,c("SurveyID","real_time")])
  order(df_new$real_time)
  df_new <- df_new[order(df_new$real_time),]
  df_new$species <- NA
  
  # These conditionals are for the selection of output.
  
  if (subset_dataframe == 1) {
    return(sub_series)
  }
  
  if (number_of_surveys == 1) {
    return(unique_surveys)
  }
  
  if (timeseries_timepoints == 1) {
    return(unique_timepoints)
  }
  
  if (time_range == 1) {
    return(output)
  }
  
  if (species_in_survey == 1) {
    return(unique_species)
  }  
  
  if (useful_df == 1) {
    # Use a loop to create a new dataframe that will include species per survey per timeframe (to avoid counting species 2x)
    for (j in 1:nrow(df_new)) {
      count <- 0
      for (i in 1:nrow(sub_series)) {
        if (sub_series[i, 8] == df_new[j, 2] & sub_series[i, 2] == df_new[j, 1]) {
          count = count + 1
        }
      }
      df_new[j, 3] <- count 
    }  
    
    # Lastly, we can create a method to "count the unique species added at each timeframe".
    # This allows us to make a quick assessment on whether or not finding a novel community here
    # will actually be likely or not.
    
    `%!in%` <- Negate(`%in%`)
    species_checklist <- subset(sub_series, sub_series$SurveyID == min(sub_series$SurveyID))
    species_checklist <- matrix(unique(species_checklist[, c("Species")]))
    df_new$new_species <- 0
    df_new[1,4] <- nrow(species_checklist)
    
    for (j in 2:nrow(df_new)) {
      count = 0
      for (i in 1:nrow(sub_series)) {
        if (sub_series[i, 2] == df_new[j, 1]) {
          if (sub_series[i, 5] %!in% species_checklist) {
            count <- count + 1
            species_checklist <- rbind(sub_series[i,5], species_checklist)
          }
        }
      }
      df_new[j, 4] <- count
    } 
    
    return(df_new)
  }
}

# Example usage
test <- metric_extractor("G1037", 0, 0, 0, 0, 0, 1)
