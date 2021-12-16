### List_novelty_B_function ###

# This function takes a list of presence/absence matrices and returns a list of novelty 
# classification matrices. 

list_novelty_B_function <- function(matrix_list){
  
  for (i in 1:length(matrix_list)) {
    
    nam <- lapply(names(matrix_list[[i]][]), 
                  function(ID){
                    print(ID)
                    site.sp.mat <- (matrix_list[[i]][ID])
                    site.sp.mat <- site.sp.mat[[ID]]
                    #site.sp.mat <- site.sp.mat[[ID]]
                    # This line had to be added because there was one timeseries with 0 change over 20 years....
                    if (ID == "G7555"){
                      return(NA)
                    }
                    else{
                      if(typeof(site.sp.mat) == "character"){
                        return(NA)
                      }
                      else{
                        if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
                          
                          temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
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
    
    if (names(matrix_list[i]) == "palearctic_mat_B"){
      palearctic_novelty_B <- nam
      names(palearctic_novelty_B) <- ID_list$palearctic_ID
      palearctic_novelty_B <- palearctic_novelty_B[!sapply(palearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "nearctic_mat_B"){
      nearctic_novelty_B <- nam
      names(nearctic_novelty_B) <- ID_list$nearctic_ID
      nearctic_novelty_B <- nearctic_novelty_B[!sapply(nearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "afrotropics_mat_B"){
      afrotropics_novelty_B <- nam
      names(afrotropics_novelty_B) <- ID_list$afrotropics_ID
      afrotropics_novelty_B <- afrotropics_novelty_B[!sapply(afrotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "neotropics_mat_B"){
      neotropics_novelty_B <- nam
      names(neotropics_novelty_B) <- ID_list$neotropics_ID
      neotropics_novelty_B <- neotropics_novelty_B[!sapply(neotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "australasia_mat_B"){
      australasia_novelty_B <- nam
      names(australasia_novelty_B) <- ID_list$australasia_ID
      australasia_novelty_B <- australasia_novelty_B[!sapply(australasia_novelty_B, function(x) all(is.na(x)))]
    }
    
    # Add a list of TimeSeries which will come in handy when evaluating the results.
    
    
  }
  list_novelty <- list(palearctic_novelty_B, nearctic_novelty_B, afrotropics_novelty_B, neotropics_novelty_B, australasia_novelty_B)
  names(list_novelty) <- c("palearctic_novelty_B","nearctic_novelty_B","afrotropics_novelty_B","neotropics_novelty_B", "australasia_novelty_B")
  return(list_novelty)
  
}