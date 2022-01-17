### list_matrix_B_bins_function ###

# This function is essentially the same as the other matrix_list functions, but it also takes a bin width (it uses the binary
# aggregate function). It takes in a list of TimeSeries ID's and a bin width to return a list of presence/absence matrices.


list_matrix_B_bins_function <- function(check_list, bin_width){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- binary_aggregate_function(TimeSeries_ID, bin_width)
                  
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
