### list_matrix_A_bins_function ###

# This function takes in a list of TimeSeries ID's and a bin width and returns a list of abundance matrices according 
# to bin width. 

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
