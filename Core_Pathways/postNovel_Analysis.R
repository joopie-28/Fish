#### Quantifying the post-novel trajectory ####
create.postnovel.groups <- function(ID){
  test.mat <- cbind(nov.matrices[[ID]], 'cat' = c(rep('back', 5), nov.output[[ID]]$cat))
  
  # Add transition groups to the main data frame
  test.mat$novel.groups <- "Non-Novel"
  
  flattened.novelty.pers <- (flatten(novelty.pers))
  
  nov.data = flattened.novelty.pers[[which(sapply(names(flattened.novelty.pers), 
                                                  function(index) is.element(names(nov.matrices[ID]), index)))[1]]]
  if(nov.data$consistency == F){return(NULL)}
  begin<- which(rownames(test.mat) == nov.data$begin)
  end<-which(rownames(test.mat) == nov.data$end)
  
  # Assign groups 
  if(end <= nrow(test.mat)){
    test.mat$novel.groups[begin:(end)] <- 'Novel'
  } else{
    test.mat$novel.groups[begin:(end-1)] <- 'Novel'
  }
  if(test.mat[nrow(test.mat),'novel.groups' ]=='Novel'){ return(NULL)}
  # Also add previous group
  prev.group <- "Non-Novel"
  for(i in 2:nrow(test.mat)){
    prev.group[i] <- test.mat$novel.groups[i-1]
  }
  
  test.mat$prev.group <- prev.group
  
  trans.index <- which(test.mat$novel.groups == 'Non-Novel' & test.mat$prev.group == 'Novel')
  nov.range <- trans.index:nrow(test.mat)
  for(i in 1:length(trans.index:nrow(test.mat))){
    test.mat$novel.groups[nov.range[i]] <- paste0("Pn+",i)
  }
  return(test.mat)
}

grouped.matrices <- lapply(names(nov.matrices), function(index){
  print(index)
  return(create.postnovel.groups(index))
})


# Now calculate the distance metrics

calculate.distance.per.group <- function(input.mat){
  if(is.null(input.mat)){return(NULL)}
  classes <- unique(input.mat$novel.groups)
  groups_df <- as.data.frame(matrix(nrow = length(classes), 
                                    ncol = ncol(input.mat|> dplyr::select(-c(cat, 
                                                                             novel.groups, 
                                                                             prev.group)))+1))
  groups_df[,1] = classes
  for(i in 1:length(classes)){
    centroid <-(sapply((subset(input.mat, novel.groups == classes[i]) |>
                          dplyr::select(-c(cat, novel.groups, prev.group))),mean))
    
    groups_df[i,2:ncol(groups_df)] <- centroid
    
    rownames(groups_df) <- groups_df[[1]]
    
  }
  
  dis.matrix <- as.data.frame(as.matrix(vegdist(groups_df[,-1], method = 'bray')))
  
  new.matrix <- dis.matrix[,1:2]
  new.matrix = round(new.matrix, digits = 2)
  new.matrix$dNdP <- ifelse(new.matrix[,2]/new.matrix[,1] == Inf, 0, 
                            new.matrix[,2]/new.matrix[,1])
  new.matrix$cat <- rownames(new.matrix)
  rownames(new.matrix) <- NULL
  return(new.matrix)
}

DisTrajectories <- rbindlist(lapply(grouped.matrices, function(mat){
  return(calculate.distance.per.group(mat))
})) |>
  mutate(Closest_State = ifelse(Novel > `Non-Novel`, 'Novel', 'Non_Novel')) 

for(i in 1:nrow(DisTrajectories )){
  DisTrajectories $cat.num[i] <- strsplit(DisTrajectories $cat, split = '[+]')[[i]][2]
}

temp <- DisTrajectories[as.numeric(DisTrajectories$cat.num)<8, ]





