#### Julian days bins ####

library(timeDate)
library(insol)


# Ensure quarters are made numeric. We jsut replace 2/3 for 2.3 as this approximates it quite well.
Survey_Data$Quarter <- gsub(Survey_Data$Quarter, 
                            pattern ="\\/", 
                            replacement = ".")

Survey_Data$Quarter <- as.numeric(Survey_Data$Quarter)

# Replace NA with 0 (this is what they signify)

Survey_Data$Quarter[is.na(Survey_Data$Quarter)] = 0

# This function will use julian days as bins and thus preserve more data.

abundance_julian_function <- function(survey_identifier){
  
  # Convert bin width to julian days.
  
  
  # Here I create some reference lists which help us create the matrices
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  
  # Create a column which signifies Julian days.
  sub_series$True_Date <- NA
  sub_series$Julian_Day <- NA

  for (i in 1:nrow(sub_series)) {
    if (sub_series$Quarter[i] == 0){
      return(NA)
    }
    else{
    sub_series$True_Date[i] <- paste0(sub_series$Year[i],"-",ceiling((as.numeric(sub_series$Quarter[i])*(3))-1),"-15")
    sub_series$Julian_Day[i] <- JD(as.POSIXct(sub_series$True_Date[i]))
    sub_series$Julian_Day[i] <- floor(sub_series$Julian_Day[i])
    }
  }
  
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  year_unique <- as.data.frame(unique(sub_series$Julian_Day))
  
  

  
  # Here I create the empty dataframe which will be populated with abundance data
  
  abundance_df <- as.data.frame(matrix(data = NA, 
                                       nrow = length(unique(sub_series$Julian_Day)), 
                                       ncol = nrow(ref_list)))
  
  colnames(abundance_df) <- ref_list[, 1]
  rownames(abundance_df) <- year_unique[,1]
  
  # These loops will populate the empty data frame; it will simply add in the value from the main dataset into the
  # correct cell based on the row/column ID's.
  
  for (i in 1:nrow(abundance_df)) {
    sub_temp <- subset(sub_series, Julian_Day == rownames(abundance_df)[i])

    
    for (j in 1:ncol(abundance_df)) {
      abundance_list <- vector()
      
      for (k in 1:nrow(sub_temp)) {
        
        if ((colnames(abundance_df[j])) == toString(sub_temp[k, 5])){
          abundance_list <- append(abundance_list, sub_temp[k, 6])
          
        }
      }
      abundance_list <- as.numeric(abundance_list)
      average <- ceiling(mean(abundance_list))
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
  
  row.names(abundance_df) <- (as.numeric(row.names(abundance_df)))
  
  # Replace NA with 0 (this is what they signify)
  
  abundance_df[is.na(abundance_df)] = 0
  
  # Now remove rows that have sums of 0, these are rows of years that were actually not present in the original data.
  
  abundance_df <- abundance_df[rowSums(abundance_df[])> 0, ]
  
  
  return(abundance_df)
  
}

julian_matrix_A_function <- function(check_list){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- abundance_julian_function(TimeSeries_ID)
                    
                    return(temp)
                  }) 
    
    if (names(check_list[i]) == "palearctic_ID"){
      palearctic_mat_A <- nam
      names(palearctic_mat_A) <- check_list$palearctic_ID
    }
    if (names(check_list[i]) == "nearctic_ID"){
      nearctic_mat_A <- nam
      names(nearctic_mat_A) <- check_list$nearctic_ID
    }
    if (names(check_list[i]) == "afrotropics_ID"){
      afrotropics_mat_A <- nam
      names(afrotropics_mat_A) <- check_list$afrotropics_ID
    }
    if (names(check_list[i]) == "neotropics_ID"){
      neotropics_mat_A <- nam
      names(neotropics_mat_A) <- check_list$neotropics_ID
    }
    if (names(check_list[i]) == "australasia_ID"){
      australasia_mat_A <- nam
      names(australasia_mat_A) <- check_list$australasia_ID
    }
    
  }
  list_matrix <- list(palearctic_mat_A, nearctic_mat_A, afrotropics_mat_A, neotropics_mat_A, australasia_mat_A)
  names(list_matrix) <- c("palearctic_mat_A","nearctic_mat_A","afrotropics_mat_A","neotropics_mat_A", "australasia_mat_A")
  return(list_matrix)
}


# test

#abundance_julian_function(survey_ID)

#a <- julian_matrix_A_function(ID_list)

#novelty_a <- list_novelty_A_function(julian_matrices_A)

