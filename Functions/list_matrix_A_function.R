### List_matrix_A_function ###

# This function creates a list of abundance matrixes, using 1-year time bins, taking a list of Timeseries ID's as input.
# The output is in the form of a list of abundance matrices, with the TimeSeries ID as 'keys'. It requires the 
# "abundance_identifier_function" to be loaded in as well.

list_matrix_A_function <- function(check_list){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- abundance_identifier_function(TimeSeries_ID)
                    if(class(temp)=="matrix"){
                      # remove bins within sites with fewer species than cut off
                      temp <- temp[rowSums(temp) > rich.cutoff,]
                    }
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
