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
  
    if (any(is.na(matrix_list[[i]]))){
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

mat.nnc.ass.V2 <- function(matrix_list){
  
  for(i in 1:length(matrix_list)) {
      
      matrix_invasive <- as.data.frame(matrix_list[[i]])
      matrix_invasive$NNC <- 0
      matrix_invasive$NNC_increase <- 0
      matrix_invasive$NNC_spec <-0
      matrix_invasive$NNC_spec_increase <-0
      matrix_invasive$NAC <- 0
      matrix_invasive$NAC_increase <- 0
      matrix_invasive$NAC_spec_increase <-0
      matrix_invasive$INC <- 0
      matrix_invasive$INC_increase <- 0
      matrix_invasive$INC_spec_increase <-0
      matrix_invasive$bins <- as.numeric(rownames(matrix_invasive))
      
      for (j in 1:nrow(matrix_invasive)) {
        NNC <- 0
        NNC_spec <- 0
        NAC <- 0
        NAC_spec <- 0
        INC <- 0
        INC_spec <- 0
        for (k in 1:ncol(matrix_invasive)) {
          
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "xotic")){
            NNC <- NNC + matrix_invasive[j,k]
            if(matrix_invasive[j,k] != 0){
              NNC_spec <- NNC_spec +1
            }
          }
          
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "ative")){
            NAC <- NAC + matrix_invasive[j,k]
            if(matrix_invasive[j,k] != 0){
              NAC_spec <- NAC_spec +1
            }
          }
          
          if(stri_detect_fixed(colnames(matrix_invasive)[k], "nvader")){
            INC <- INC + matrix_invasive[j,k]
            if(matrix_invasive[j,k] != 0){
              INC_spec <- INC_spec +1
            }
          }
          
        }
        matrix_invasive$NNC[j] <- NNC
        matrix_invasive$NAC[j] <- NAC
        matrix_invasive$INC[j] <- INC
        matrix_invasive$NNC_spec[j] <- NNC_spec
        matrix_invasive$NAC_spec[j] <- NAC_spec
        matrix_invasive$INC_spec[j] <- INC_spec
        
        if (j >= 2){
          
          # We can not work with percentages here so need to think of a way that conveys
          # magnitude of change.
          
          matrix_invasive$NNC_increase[j] <- (matrix_invasive$NNC[j] - matrix_invasive$NNC[j-1])
          matrix_invasive$NAC_increase[j] <- (matrix_invasive$NAC[j] - matrix_invasive$NAC[j-1])
          matrix_invasive$INC_increase[j] <- (matrix_invasive$INC[j] - matrix_invasive$INC[j-1])
          matrix_invasive$NNC_spec_increase[j] <- (matrix_invasive$NNC_spec[j] - matrix_invasive$NNC_spec[j-1])
          matrix_invasive$NAC_spec_increase[j] <- (matrix_invasive$NAC_spec[j] - matrix_invasive$NAC_spec[j-1])
          matrix_invasive$INC_spec_increase[j] <- (matrix_invasive$INC_spec[j] - matrix_invasive$INC_spec[j-1])
        }
      }
      
      
      
      matrix_list[[i]] <- matrix_invasive
    
  }
  return(matrix_list)
  
}



mat.nnc.ass.V3 <- function(matrix_list){
  
  for(i in 1:length(matrix_list)) {
    
    matrix_invasive <- as.data.frame(matrix_list[[i]])
    matrix_invasive$NNC <- 0
    matrix_invasive$NNC_increase <- 0
    matrix_invasive$NNC_spec <-0
    matrix_invasive$NNC_spec_increase <-0
    matrix_invasive$NAC <- 0
    matrix_invasive$NAC_increase <- 0
    matrix_invasive$NAC_spec_increase <- 0
    matrix_invasive$bins <- as.numeric(rownames(matrix_invasive))
    
    for (j in 1:nrow(matrix_invasive)) {
      NNC <- 0
      NNC_spec <- 0
      NAC <- 0
      NAC_spec <- 0
      INC <- 0
      INC_spec <- 0
      for (k in 1:ncol(matrix_invasive)) {
        
        if(stri_detect_fixed(colnames(matrix_invasive)[k], "xotic")){
          NNC <- NNC + matrix_invasive[j,k]
          if(matrix_invasive[j,k] != 0){
            NNC_spec <- NNC_spec +1
          }
        }
        
        if(stri_detect_fixed(colnames(matrix_invasive)[k], "ative")){
          NAC <- NAC + matrix_invasive[j,k]
          if(matrix_invasive[j,k] != 0){
            NAC_spec <- NAC_spec +1
          }
        }
        
        
      }
      matrix_invasive$NNC[j] <- NNC
      matrix_invasive$NAC[j] <- NAC
  
      matrix_invasive$NNC_spec[j] <- NNC_spec
      matrix_invasive$NAC_spec[j] <- NAC_spec

      
      if (j >= 2){
        
        # We can not work with percentages here so need to think of a way that conveys
        # magnitude of change.
        
        matrix_invasive$NNC_increase[j] <- (matrix_invasive$NNC[j] - matrix_invasive$NNC[j-1])
        matrix_invasive$NAC_increase[j] <- (matrix_invasive$NAC[j] - matrix_invasive$NAC[j-1])
        matrix_invasive$NNC_spec_increase[j] <- (matrix_invasive$NNC_spec[j] - matrix_invasive$NNC_spec[j-1])
        matrix_invasive$NAC_spec_increase[j] <- (matrix_invasive$NAC_spec[j] - matrix_invasive$NAC_spec[j-1])
      }
    }
    
    
    
    matrix_list[[i]] <- matrix_invasive
    
  }
  return(matrix_list)
  
}


