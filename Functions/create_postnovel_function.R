# Create post-novel groupings
classify.compositions <- function(matrix, pers.data){
  
  # Extract novelty clustering information 
  matrix$novel.groups <- "Non-Novel"
  begin<- which(matrix$bins == pers.data$begin)
  end<-which(matrix$bins == pers.data$end)
  
  # Assign groups 
  if(end < nrow(matrix)){
    matrix$novel.groups[begin:(end-1)] <- 'Novel'
  } else{
    matrix$novel.groups[begin:(end)] <- 'Novel'
  }
  
  prev.group <- "Non-Novel"
  for(i in 2:nrow(matrix)){
    prev.group[i] <- matrix$novel.groups[i-1]
  }
  matrix$prev.group <- prev.group
  
  return(matrix)
}
