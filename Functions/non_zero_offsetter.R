# function for offsetting rows with a lot of 0's and RA of 1 (monospecific)

non_zero_offsetter <- function(matrix){
  
  for (row.index in 1:nrow(matrix)){
    
    row = matrix[row.index,]
    
    if (any(row == 1)){
      
      perfect.abundance <- which(row == 1)
      
      row[perfect.abundance] <- row[perfect.abundance] - 0.001*(length(row)-1)
      
      for(i in 1:length(row)){
        
        if (i != perfect.abundance){
          row[i] = row[i] + 0.001
        }
      }
    }
    matrix[row.index, ] <- row 
  }
  return(matrix)
}
