### abundance_aggregate_function ###

# This function takes a TimeSeriesID and a bin width as input and returns a matrix of abundance data. Where bin width is larger
# than 1, values within the matrix are average abundance over that time period.

abundance_aggregate_function <- function(survey_identifier, bin_width){
  
  # Here I create some reference lists which help us create the matrices
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
  
  if (bin_width == 1){
    if (nrow(year_unique) < 10){
      return ("Not enough data when bin size is 1 years (less than 10 years)")
    }
  }
  
  # Here I create the empty dataframe which will be populated with abundance data
  
  abundance_df <- as.data.frame(matrix(data = NA, 
                                       nrow = (ceiling(((max_year-min_year+1)/bin_width))), 
                                       ncol = nrow(ref_list)))
  
  colnames(abundance_df) <- ref_list[, 1]
  
  # Select the correct years as a function of bin width
  
  name_list <- list()
  name_list <- seq(min_year, max_year, by = bin_width)
  
  name_list <- name_list[order(as.numeric((name_list)), 
                               decreasing = FALSE)]
  
  rownames(abundance_df) <- name_list
  
  
  # These loops will populate the empty data frame; it will simply add in the value from the main dataset into the
  # correct cell based on the row/column ID's.
  
  for (i in 1:nrow(abundance_df)) {
    sub_temp <- subset(sub_series, Year >= as.numeric(rownames(abundance_df)[i]) & 
                         Year < (as.numeric(row.names(abundance_df)[i]) + (bin_width)))
    
    for (j in 1:ncol(abundance_df)) {
      abundance_list <- vector()
      
      for (k in 1:nrow(sub_temp)) {
        
        if ((colnames(abundance_df[j])) == toString(sub_temp[k, 5])){
          abundance_list <- append(abundance_list, sub_temp[k, 6])
          
        }
      }
      abundance_list <- as.numeric(abundance_list)
      average <- mean(abundance_list)
      abundance_df[i,j] <- average
    }
  }
  
  # Order dataframe so that novelty framework does not get confused
  
  abundance_df <- abundance_df[order(as.numeric((row.names(abundance_df))), 
                                     decreasing = FALSE),]
  
  # Very small samples will be returned as vectors rather than dataframes, so will instruct the function to skip this survey if
  # That happens Shouldn't be necessary based on initial filter conditionals.
  
  if (typeof(abundance_df) == "double"){
    return("Not enough data")
  }
  
  # Change row names to bin (time)
  
  row.names(abundance_df) <- (2021-as.numeric(row.names(abundance_df)))
  
  # Replace NA with 0 (this is what they signify)
  
  abundance_df[is.na(abundance_df)] = 0
  
  # Now remove rows that have sums of 0, these are rows of years that were actually not present in the original data.
  
  abundance_df <- abundance_df[rowSums(abundance_df[])> 0, ]
  
  
  return(abundance_df)
  
}
