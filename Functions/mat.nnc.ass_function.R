# This function calculates and assigns the contribution of native and non-native species 
# in a community 

mat.nnc.ass <- function(matrix_list){

  for(i in 1:length(matrix_list)) {
    print(i)
    
    if (!is.na(matrix_list[i])){
      if (ncol(matrix_list[[i]]) < 5){
        matrix_list[[i]] <- NA
      }
      else{
        if(nrow(matrix_list[[i]]) < 10){
          matrix_list[[i]] <- NA
        }
      }
    }
  
    if (is.na(matrix_list[[i]])){
      matrix_list[[i]] <- NA
    }
    
    else{


    
      matrix_invasive <- as.data.frame(matrix_list[[i]])
      matrix_invasive$NNC <- 0
      matrix_invasive$NNC_increase <- 0
      matrix_invasive$NAC <- 0
      matrix_invasive$NAC_increase <- 0
      matrix_invasive$INC <- 0
      matrix_invasive$INC_increase <- 0
      matrix_invasive$bins <- as.numeric(rownames(matrix_invasive))
    
      for (j in 1:nrow(matrix_invasive)) {
        NNC <- 0
        NAC <- 0
        INC <- 0
        for (k in 1:ncol(matrix_invasive)) {
        
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "xotic")){
            NNC <- NNC + matrix_invasive[j,k]
          
          }
        
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "ative")){
            NAC <- NAC + matrix_invasive[j,k]
          }
          
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "nvader")){
            INC <- INC + matrix_invasive[j,k]
          }
        
        }
        matrix_invasive$NNC[j] <- NNC
        matrix_invasive$NAC[j] <- NAC
        matrix_invasive$INC[j] <- INC
       
        if (j >= 2){
          
          # We can not work with percentages here so need to think of a way that conveys
          # magnitude of change.
          
          matrix_invasive$NNC_increase[j] <- (matrix_invasive$NNC[j] - matrix_invasive$NNC[j-1])
          matrix_invasive$NAC_increase[j] <- (matrix_invasive$NAC[j] - matrix_invasive$NAC[j-1])
          matrix_invasive$INC_increase[j] <- (matrix_invasive$INC[j] - matrix_invasive$INC[j-1])
        }
      }
      
   
      
      matrix_list[[i]] <- matrix_invasive
    }
  }
  return(matrix_list)

}
