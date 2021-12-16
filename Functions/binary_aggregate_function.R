### Binary_aggregate_function ###

# This function takes a survey identifier (e.g. "G1034") and a bin width (1 and 2 are the best, any higher and most TimeSeries
# become too short) to return a matrix of binary data.

binary_aggregate_function <- function(survey_identifier, bin_width){
  
  # Here I create some reference lists 
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  year_unique <- as.data.frame(unique(sub_series$Year))
  min_year <- min(sub_series$Year)
  max_year <- max(sub_series$Year)
  
  # We are now condensing the data quite significantly, so I'm going to write off smaller surveys
  if (bin_width == 4){
    if (nrow(year_unique) < 24){
      return("Not enough data when bin size is 4 years")
    }
  }
  
  if (bin_width == 2){
    if (nrow(year_unique) < 12){
      return ("Not enough data when bin size is 2 years")
    }
  }
  
  # Here I create the empty dataframe which will be populated with presence/absence values, based on the bin width
  binary_df <- as.data.frame(matrix(data = NA, nrow = (ceiling(((max_year-min_year+1)/bin_width))), ncol = nrow(ref_list)))
  colnames(binary_df) <- ref_list[, 1]
  
  # Select the correct years as a function of bin width
  if (bin_width != 1){
    name_list <- list()
    name_list <- seq(min_year, max_year, by = bin_width)
    name_list <- name_list[order(as.numeric((name_list)), decreasing = FALSE)]
    rownames(binary_df) <- name_list
  }
  
  # Let's populate this data frame. I will use the column and row names to loop through a subset of the data
  # a couple times. This will essentially crosscheck if a species is present in the set or not. 
  for (i in 1:nrow(binary_df)) {
    sub_temp <- subset(sub_series, Year >= as.numeric(rownames(binary_df)[i]) & Year < (as.numeric(row.names(binary_df)[i]) + (bin_width)))
    for (j in 1:ncol(binary_df)) {
      if (colnames(binary_df[j]) %in% sub_temp$Species){
        binary_df[i,j] <- 1
      }
      else{
        binary_df[i,j] <- 0
      }
    }
  }
  
  # Order dataframe so that novelty framework does not get confused
  binary_df <- binary_df[order(as.numeric((row.names(binary_df))), decreasing = FALSE),]
  
  # Very small samples will be returned as a string of numbers (type = double) rather than data frames, so will instruct the function to skip this survey if
  # That happens. The timeseries is just too short basically.
  if (typeof(binary_df) == "double"){
    return("Not enough data")
  }
  
  # Change row names to bin (time)
  row.names(binary_df) <- (2021-as.numeric(row.names(binary_df)))
  
  # Now remove rows that have sums of 0, these are rows of years that were actually not present in the original data.
  binary_df <- binary_df[rowSums(binary_df[])> 0, ]
  
  # Make sure the order is OK, it will influence the model results (our priamry reference community HAS to be the oldest one for the
  # model to make sense).
  
  return(binary_df)
}