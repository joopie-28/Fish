### Building matrices based on varying bin widths as opposed to just taking time as is ###

# We are going to test different bin widths: 2 and 4 years.

# This is for abundance

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
  
  # Here I create the empty dataframe which will be populated with abundance data
  abundance_df <- as.data.frame(matrix(data = NA, nrow = (ceiling(((max_year-min_year+1)/bin_width))), ncol = nrow(ref_list)))
  colnames(abundance_df) <- ref_list[, 1]
  
  # Select the correct years as a function of bin width
  if (bin_width != 1){
    name_list <- list()
    name_list <- seq(min_year, max_year, by = bin_width)
    name_list <- name_list[order(as.numeric((name_list)), decreasing = FALSE)]
    rownames(abundance_df) <- name_list
  }
    
  # These loops will populate the empty data frame; it will simply add in the value from the main dataset into the
  # correct cell based on the row/column ID's.
  for (i in 1:nrow(abundance_df)) {
    sub_temp <- subset(sub_series, Year >= as.numeric(rownames(abundance_df)[i]) & Year < (as.numeric(row.names(abundance_df)[i]) + (bin_width)))
    for (j in 1:ncol(abundance_df)) {
      abundance_list <- vector()
      for (k in 1:nrow(sub_temp)) {
        if ((colnames(abundance_df[j])) == toString(sub_temp[k, 5])){
          abundance_list <- append(abundance_list, sub_temp[k, 6])
        
        }
      }
      average <- mean(abundance_list)
      abundance_df[i,j] <- average
    }
  }

  
  abundance_df <- abundance_df[order(as.numeric((row.names(abundance_df))), decreasing = FALSE),]
  
  # Very small samples will be returned as vectors rather than dataframes, so will instruct the function to skip this survey if
  # That happens
  if (typeof(binary_df) == "double"){
    return("Not enough data")
  }
  
  # Change row names to bin (time)
  row.names(abundance_df) <- (2021-as.numeric(row.names(abundance_df)))
  
  # Replace NA with 0 (this is what they signify)
  abundance_df[is.na(abundance_df)] = 0
  

  return(abundance_df)
  
}

list_matrix_A_bins_function <- function(check_list, bin_width){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- abundance_aggregate_function(TimeSeries_ID, bin_width)
                    if(class(temp)=="matrix"){
                      # remove bins within sites with fewer species than cut off
                      temp <- temp[rowSums(temp) > rich.cutoff,]
                    }
                    return(temp)
                  }) 
    
    if (names(check_list[i]) == "palearctic_ID"){
      palearctic_mat_B <- nam
      names(palearctic_mat_B) <- check_list$palearctic_ID
    }
    if (names(check_list[i]) == "nearctic_ID"){
      nearctic_mat_B <- nam
      names(nearctic_mat_B) <- check_list$nearctic_ID
    }
    if (names(check_list[i]) == "afrotropics_ID"){
      afrotropics_mat_B <- nam
      names(afrotropics_mat_B) <- check_list$afrotropics_ID
    }
    if (names(check_list[i]) == "neotropics_ID"){
      neotropics_mat_B <- nam
      names(neotropics_mat_B) <- check_list$neotropics_ID
    }
    if (names(check_list[i]) == "australasia_ID"){
      australasia_mat_B <- nam
      names(australasia_mat_B) <- check_list$australasia_ID
    }
    
  }
  list_matrix <- list(palearctic_mat_B, nearctic_mat_B, afrotropics_mat_B, neotropics_mat_B, australasia_mat_B)
  names(list_matrix) <- c("palearctic_mat_B","nearctic_mat_B","afrotropics_mat_B","neotropics_mat_B", "australasia_mat_B")
  return(list_matrix)
}

matrix_list_A_bin2 <- list_matrix_B_bins_function(ID_list, 2)
 

novelty_list_A_bin2 <- list_novelty_function(matrix_list_A_bin2)

novelty_analysis_output_A_bin2 <- novel.probability(novelty_list_A_bin2)



venn_plot_function(novelty_analysis_output_A_bin2)
