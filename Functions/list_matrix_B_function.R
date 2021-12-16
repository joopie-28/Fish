### List_matrix_B_function ###

# This function creates a list of presence/absence (binary) matrixes, using 1-year time bins, taking a list of Timeseries ID's as input.
# The output is in the form of a list of binary matrices, with the TimeSeries ID as 'keys'. It requires the 
# "binary_converter_function" to be loaded in as well.

list_matrix_B_function <- function(check_list){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- binary_converter_function(TimeSeries_ID)
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