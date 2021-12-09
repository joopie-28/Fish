# Identifier function abundance # Exactly the same as the binary function only now we use the abundance data (which
# requires a slightly different manipulation pathway.

abundance_identifier_function <- function(survey_identifier){
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  
  # Loading the previously used function so that we can use them separately
  abundance_df_function <- function(survey_identifier){
    
    # Here I create some reference lists which help us create the matrices
    sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
    ref_list <- as.data.frame(unique(sub_series$Species))
    surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
    
    # Here I create the empty dataframe which will be populated with abundance data
    abundance_df <- as.data.frame(matrix(data = 0, nrow = nrow(unique(surveys_unique)), ncol = nrow(ref_list)))
    colnames(abundance_df) <- ref_list[, 1]
    rownames(abundance_df) <- surveys_unique[, 1]
    
    # These loops will populate the empty data frame; it will simply add in the value from the main dataset into the
    # correct cell based on the row/column ID's.
    for (i in 1:nrow(abundance_df)) {
      sub_temp <- subset(sub_series, SurveyID == rownames(abundance_df)[i])
      for (j in 1:ncol(abundance_df)) {
        for (k in 1:nrow(sub_temp)) {
          if ((colnames(abundance_df[j])) == toString(sub_temp[k, 5])){
            abundance_df[i,j] <- sub_temp[k, 6]
          }
        }
      }
    }
    return(abundance_df)
  }
  
  # Changing the name here because I'm recycling code I wrote earlier, a bit messy but saves time. Binary_df is just a 
  # placeholder variable.
  binary_df <- abundance_df_function(survey_identifier)
  
  
  # This part creates a subset of ID, year and quarter which will be used to create distinct timepoints.
  unique_df <- subset(unique(sub_series[,c("SurveyID","Year", "Quarter")]))
  unique_df$real_time <- 0
  
  # Here I ensure no NA's are produced to entries like "quarter 2/3" where we would expect a single number
  for (i in 1:nrow(unique_df)) {
    if (typeof(unique_df[i, 3]) != "character"){
      unique_df[i, 4] <- (as.numeric(unique_df[i,2]) + as.numeric(unique_df[i,3])*0.225)
    }
    else {
      unique_df[i, 4] <- as.numeric(unique_df[i, 2])
    }
  }
  
  # Renaming the rownames from survey ID to respective timepoints (necessary for the Pandolfi model)
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
  if (nrow(binary_df) >= 5){
  binary_df <- subset(binary_df, quarter == as.numeric(tail(names(sort(table(binary_df$quarter))), 1)))
  }
  binary_df <- binary_df[, !(names(binary_df) == "quarter")]

  # Arbitrarily removing duplicate surveys based on time frame. In accordance with Compte et al 2021
  binary_df <- binary_df[!duplicated(binary_df[,c("time")]),]

  # Change row names to bin (time)
  row.names(binary_df) <- (2021-as.numeric(binary_df$time))

  # Now remove the time column
  binary_df <- binary_df[, !(names(binary_df) == "time")]
  
  # Very small samples will be returned as vectors rather than dataframes, so will instruct the function to skip this survey if
  # That happens
  if (typeof(binary_df) == "double"){
    return("Not enough data")
  }
  
  # Make sure the order is OK, it will influence the model results
  binary_df <- binary_df[order(as.numeric(row.names(binary_df)), decreasing = FALSE),]

  return(binary_df)
  
}  





