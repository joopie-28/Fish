### List_novelty_A_function ###

# This function takes a list of abundance matrices as input and returns a list of novelty detection matrices. Ensure 
# that the function "identify.novel.gam" is loaded.


list_novelty_A_function <- function(matrix_list){
  
  for (i in 1:length(matrix_list)) {
    
    nam <- lapply(names(matrix_list[[i]][]), 
                  function(ID){
                    print(ID)
                    site.sp.mat <- (matrix_list[[i]][ID])
                    site.sp.mat <- site.sp.mat[[ID]]
                    #site.sp.mat <- site.sp.mat[[ID]]
                    # This line had to be added because there was one timeseries with 0 change over 20 years....
                    if (is.na(site.sp.mat)){
                      return(NA)
                    }
                    else{
                      if(typeof(site.sp.mat) == "character"){
                        return(NA)
                        
                      }
                      else{
                        if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
                          
                          temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
                                                     gam.max.k = -1)
                          # Remove first 5 bins
                          temp <- temp[-c(1:5),]
                          return(temp)
                        }
                        else {
                          return(NA)
                        }
                      }
                    }  
                  })
    
    # Select correct data, delete NA's
    
    if (names(matrix_list[i]) == "palearctic_mat_A"){
      palearctic_novelty_A <- nam
      names(palearctic_novelty_A) <- ID_list$palearctic_ID
      palearctic_novelty_A <- palearctic_novelty_A[!sapply(palearctic_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "nearctic_mat_A"){
      nearctic_novelty_A <- nam
      names(nearctic_novelty_A) <- ID_list$nearctic_ID
      nearctic_novelty_A <- nearctic_novelty_A[!sapply(nearctic_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "afrotropics_mat_A"){
      afrotropics_novelty_A <- nam
      names(afrotropics_novelty_A) <- ID_list$afrotropics_ID
      afrotropics_novelty_A <- afrotropics_novelty_A[!sapply(afrotropics_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "neotropics_mat_A"){
      neotropics_novelty_A <- nam
      names(neotropics_novelty_A) <- ID_list$neotropics_ID
      neotropics_novelty_A <- neotropics_novelty_A[!sapply(neotropics_novelty_A, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "australasia_mat_A"){
      australasia_novelty_A <- nam
      names(australasia_novelty_A) <- ID_list$australasia_ID
      australasia_novelty_A <- australasia_novelty_A[!sapply(australasia_novelty_A, function(x) all(is.na(x)))]
    }
    
    # Add a list of Time Series which will come in handy when evaluating the results.
    
    
  }
  list_novelty <- list(palearctic_novelty_A, nearctic_novelty_A, afrotropics_novelty_A, neotropics_novelty_A, australasia_novelty_A)
  names(list_novelty) <- c("palearctic_novelty_A","nearctic_novelty_A","afrotropics_novelty_A","neotropics_novelty_A", "australasia_novelty_A")
  return(list_novelty)
  
}
