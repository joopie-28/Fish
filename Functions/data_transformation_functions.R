# This function will create a dataframe with binary entries denoting species presence or absence. The #rows is equal to the #surveys
# for a timeseries and the #columns is equal to the number of unique species in that timeseries. The function takes a timeseries ID as 
# input (e.g. G1034) and also requires the timeseries and survey datafiles to be loaded.

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
    sub_temp <- subset(sub_series, SurveyID == rownames(binary_df[i,]))
    for (j in 1:ncol(binary_df)) {
      if (colnames(binary_df[j]) %in% sub_temp$Species){
        binary_df[i,j] <- 1
      }
      else{
        binary_df[i,j] <- 0
      }
    }
  }
  
  return(binary_df)
}


# In addition, if we do not want to go with a binary approach after all, here is a similar function which simply produces
# a data frame with abundance. Works in a similar way but is a bit slower because it is assigning values through looping 
# through the survey data. Again, ensure the data files are loaded in properly.

abundance_df_function <- function(survey_identifier){
  
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  
  abundance_df <- as.data.frame(matrix(data = 0, nrow = nrow(unique(surveys_unique)), ncol = nrow(ref_list)))
  colnames(abundance_df) <- ref_list[, 1]
  rownames(abundance_df) <- surveys_unique[, 1]
  
  
  for (i in 1:nrow(abundance_df)) {
    sub_temp <- subset(sub_series, SurveyID == rownames(abundance_df[i,]))
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


# Finally, one last idea which might be useful down the track. Creating a data frame that shows the relative proportion of each species within
# a community. We do not consider biomass etc here, so our metrics are solely based on number of individuals (or whatever crazy metric the fishermen used.)
# It should make it easier to see if a community is dominated by one or a few species. Note; relative abundances given in percentages.

proportion_df_function <- function(survey_identifier){
  
  # Loading the previously used function so that we can use them separately
  abundance_df_function <- function(survey_identifier){
    
    sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
    ref_list <- as.data.frame(unique(sub_series$Species))
    surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
    
    abundance_df <- as.data.frame(matrix(data = 0, nrow = nrow(unique(surveys_unique)), ncol = nrow(ref_list)))
    colnames(abundance_df) <- ref_list[, 1]
    rownames(abundance_df) <- surveys_unique[, 1]
    
    
    for (i in 1:nrow(abundance_df)) {
      sub_temp <- subset(sub_series, SurveyID == rownames(abundance_df[i,]))
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
  proportion_df <- abundance_df_function(survey_identifier)
  sum_vector <- apply(proportion_df, 1, sum)
  
  # Loop to replace absolute counts with relative abundance
  for (i in 1:nrow(proportion_df)){
    for (j in 1:ncol(proportion_df)) {
      proportion_df[i,j] <- round((proportion_df[i,j]/sum_vector[i]) * 100, digits = 2)
    }
  }
  return(proportion_df)
}
 




