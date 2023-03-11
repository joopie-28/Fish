#### Quantifying the post-novel trajectory ####

# Add transition groups to the main data frame
full.novel.mat.season$novel.groups <- "Non-Novel"
full.novel.mat.season$prev.group <- "Non-Novel"

for (i in 1:length(novelty.pers)){
  site <- names(novelty.pers[[i]])
  print(i)
  print(site)
  matrix<-full.novel.mat.season[which(full.novel.mat.season$site == site),]
  nov.data <- novelty.pers[[i]][[site]]
  matrix.expanded<-classify.compositions(matrix, pers.data=nov.data)
  
  full.novel.mat.season[which(full.novel.mat.season$site == site),] <- matrix.expanded
}

# pre-novel - novel - post-1,post-2,post-3,post-4, etc

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

for (i in 1:length(novelty.pers)){
  site <- names(novelty.pers[[i]])
  print(i)
  print(site)
  matrix<-full.novel.mat.season[which(full.novel.mat.season$site == site),]
  nov.data <- novelty.pers[[i]][[site]]
  matrix.expanded<-classify.compositions(matrix, pers.data=nov.data)
  
  return.mat <- nov.matrices[[site]]
  
  return.mat$novel_group <- c(rep('Non-Novel',5), matrix.expanded$novel.groups)
  
}