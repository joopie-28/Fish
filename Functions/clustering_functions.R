### Clustering algorithm determining persistence length ####
## Based on type I SIMPROF (Clarke et al. 2008)

nov.cluster.id.V3 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  number.names <- as.numeric(rownames(matrix))
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  novel_frame <- novel_frame[1,]
  
  # Find year and row data
  novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <<- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance ="czekanowski", 
            method.cluster = "average", alpha=.05), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.01) })
  
  # Initialize parameters for the length calculations
  # Initialize a class variable
  length <- 0
  Class <- "NONE"
  run <- FALSE
  # Plot dendogram result
  par(mfrow = c(1,1))
  
  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <<- TRUE
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(clust.data, hang=-1,
         main = "Blip Dendrogram", xlab = "")
  })
  
  # This will activate the BLIP module and assign the correct 
  # category
  if(run){
    length <- 1
    start <- novel_bin
    end <- novel_bin
    Class <- "BLIP"
  }
  
  
  # Start of the length calculation module
  if(length == 0){
    
    length.output<-length.calculator(test,novel_bin)
    
    length <- length.output$length 
    start <- length.output$start
    end <- length.output$end
  }
  # Initialize a class variable
  
  
  # Allocate blips needs to be done!!!
  if(start == end){
    
    Class <- "BLIP"
    
    
  }
  
  # Test for full persistence
  if(Class != "BLIP"){
    
    Class <- "Persister"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  
  
  
  # Return data
  return.data <- list(c(list("ID" = ID, "begin" =start, 
                             "end" =end, "length"= length, 
                             "Class" = Class, "clusters" =no.clus)))
  names(return.data) <- ID
  return(return.data)
  
}

nov.cluster.id.V4 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  #number.names <- as.numeric(rownames(matrix[-c(1:5),]))
  number.names <- as.numeric(rownames(matrix))
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  # remove first (last) 5 bins as to not skew data (novelty framework ignores these)
  #matrix <- matrix[-c((nrow(matrix)-4):nrow(matrix)),]
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  if(nrow(novel_frame) > 1){
    novel_frame <- novel_frame[2,]
  }else{
    novel_frame <- novel_frame[1,]
  }
  # Find year and row data
  novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <<- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance =vegdist, 
            method.cluster = "average", alpha=0.05), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.05) })
  
  # Initialize parameters for the length calculations
  # Initialize a class variable
  length <- 0
  Class <- "NONE"
  run <- FALSE
  # Plot dendogram result
  par(mfrow = c(1,1))
  
  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <<- TRUE
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(clust.data, hang=-1,
         main = "Blip Dendrogram", xlab = "")
  })
  
  # This will activate the BLIP module and assign the correct 
  # category
  if(run){
    length <- 1
    start <- novel_bin
    end <- novel_bin
    length.bins <- 1
    Class <- "BLIP"
  }
  
  
  # Start of the length calculation module
  if(length == 0){
    
    length.output<-length.calculator.V2(test,novel_bin)
    
    length <- length.output$length 
    start <- length.output$start
    end <- length.output$end
    
    bin.start <- which(rev(number.names) == start)
    bin.end <- which(rev(number.names) == end)
    length.bins <- bin.start-bin.end
    
    
  }
  # Initialize a class variable
  
  
  # Allocate blips needs to be done!!!
  if(start == end | length.bins == 1){
    
    Class <- "BLIP"
    
    
  }
  
  # Test for full persistence
  if(Class != "BLIP"){
    
    Class <- "Persister"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  
  
  
  # Return data
  return.data <- list(c(list("ID" = ID, "begin" =start, 
                             "end" =end, "length"= length, 
                             "Class" = Class, "clusters" =no.clus,
                             "length.bins" = length.bins)))
  names(return.data) <- ID
  return(return.data)
  
}


length.calculator.V2 <- function(test, novel_bin){
  
  # Extract all dendrogram labels
  names.vec <- unlist(list(test$significantclusters))
  
  # Construct a new named vector with communities assigned to significant clusters
  new.vec <- NULL
  for(i in 1:length(names.vec)){
    for(j in 1:length(test$significantclusters)){
      if (names.vec[i] %in% test$significantclusters[[j]]){
        new.vec[i] <- j
        names(new.vec)[i] <- names.vec[i]
      }
    }
  }
  
  # Ensure community labels are sorted in descending order.
  # Use to reorder the cluster vector.
  sorted.index <- as.character(sort(as.numeric(names(new.vec)), decreasing  =T))
  new.vec <- new.vec[sorted.index]
  
  # Extract the vector index for the novel community
  novel_bin_index <- which(names(new.vec) == novel_bin)
  
  # Extract the novel group
  group.novel<-new.vec[[novel_bin_index]]
  
  
  # Loop through the vector of clusters to find the start and end of the
  # novel composition
  start <- novel_bin
  end <- NULL
  
  for(i in 1:length(new.vec)){
    
    # This reads as: If we bump into a non-novel cluster after having
    # already passed the novel community, then the novel state has ended.
    if(new.vec[i] %!in% group.novel & as.numeric(names(new.vec)[i]) < start){
      end <- as.numeric(names(new.vec)[i])
      break
    }
    
    # If we dont meet these requirements, the novel community doesn't end
    if(is.null(end)){
      end <- as.numeric(names(new.vec)[length(new.vec)])
    }
    
  }
  
  # Correct for novel communities at very end
  
  if(novel_bin == as.numeric(names(new.vec[length(new.vec)]))){
    end <- novel_bin
  }
  
  
  # Persistence length is simply the start of the novel state until we 
  # enter something significantly different.
  length <- start-end
  
  return(list("length" = length,
              "start" = start,
              "end" = end,
              'clustervec' = new.vec))
}


nov.cluster.id.V5 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][2]
  matrix <- matrix[[1]]
  #number.names <- as.numeric(rownames(matrix[-c(1:5),]))
  number.names <- as.numeric(rownames(matrix))
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  # remove first (last) 5 bins as to not skew data (novelty framework ignores these)
  #matrix <- matrix[-c((nrow(matrix)-4):nrow(matrix)),]
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  if(nrow(novel_frame) > 1){
    iter.list <- c(1,2)
    print('double event')
  }else{
    iter.list <- 1
  }
 
  
  return.dat <- lapply(iter.list, function(x){
  # Find year and row data
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  novel_frame <- novel_frame[x,]    
  
  novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
  
  # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
  matrix.temp <- matrix
  rownames(matrix.temp) <- rev(number.names)
  
  # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
  # columns after removing the first 5 rows. This is rare and is addressed by using
  # euclidean distance for those cases instead
  
  test <<- tryCatch(
    simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
            num.simulated = 999, method.distance =vegdist, 
            method.cluster = "average", alpha=0.05), 
    error=function(e) {
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = TRUE,
              num.simulated = 999, method.distance ="euclidean", 
              method.cluster = "average", alpha=.05) })
  
  # Initialize parameters for the length calculations
  # Initialize a class variable
  length <- 0
  Class <- "NONE"
  run <- FALSE
  # Plot dendogram result
  par(mfrow = c(1,1))
  
  # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
  # test. These are immediate blips.
  tryCatch(simprof.plot(test), error=function(e) {
    print("Blip Detected")
    run <<- TRUE
    
    # We will also plot the dendrogram, without SIMPROF coloration.
    # Just a nice visualisation
    dev.off()
    clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
    par(mar = c(3,3,3,3))
    plot(clust.data, hang=-1,
         main = "Blip Dendrogram", xlab = "")
  })
  
  # This will activate the BLIP module and assign the correct 
  # category
  if(run){
    length <- 1
    start <- novel_bin
    end <- novel_bin
    length.bins <- 1
    Class <- "BLIP"
  }
  
  
  # Start of the length calculation module
  if(length == 0){
    
    length.output<-length.calculator.V2(test,novel_bin)
    
    length <- length.output$length 
    start <- length.output$start
    end <- length.output$end
    
    bin.start <- which(rev(number.names) == start)
    bin.end <- which(rev(number.names) == end)
    length.bins <- bin.start-bin.end
    
    
  }
  # Initialize a class variable
  
  
  # Allocate blips needs to be done!!!
  if(start == end | length.bins == 1){
    
    Class <- "BLIP"
    
    
  }
  
  # Test for full persistence
  if(Class != "BLIP"){
    
    Class <- "Persister"
    
  }
  
  if(novel_bin == number.names[length(number.names)]){
    Class <- "END"
    length <- 0
  }
  
  
  
  
  # Return data
  return.data <- list(c(list("ID" = ID, "begin" =start, 
                             "end" =end, "length"= length, 
                             "Class" = Class, "clusters" =no.clus,
                             "length.bins" = length.bins)))
  names(return.data) <- ID
  return(return.data)
  })
  return(return.dat)

}


nov.cluster.id.V6 <- function(matrix){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][1]
  full.ID <- names(matrix)
  matrix <- matrix[[1]]
  #number.names <- as.numeric(rownames(matrix[-c(1:5),]))
  number.names <- as.numeric(rownames(matrix))
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  # remove first (last) 5 bins as to not skew data (novelty framework ignores these)
  #matrix <- matrix[-c((nrow(matrix)-4):nrow(matrix)),]
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  if(nrow(novel_frame) > 1){
    iter.list <- c(1,2)
    print('double event')
  }else{
    iter.list <- 1
  }
  
  
  return.dat <- lapply(iter.list, function(x){
    # Find year and row data
    novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
    novel_frame <- novel_frame[x,]    
    
    novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
    
    novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
    
    # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
    matrix.temp <- matrix
    rownames(matrix.temp) <- rev(number.names)
    
    # Sometimes using bray-curtis/czekanowski can lead to an error where there are 0
    # columns after removing the first 5 rows. This is rare and is addressed by using
    # euclidean distance for those cases instead
    
    test <<- tryCatch(
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = T,
              num.simulated = 999, method.distance =vegdist, 
              method.cluster = "average", alpha=0.05), 
      error=function(e) {})
    
    # Initialize parameters for the length calculations
    # Initialize a class variable
    length <- 0
    Class <- "NONE"
    run <- FALSE
    # Plot dendogram result
    par(mfrow = c(1,1))
    
    # Use a TryCatch expression, as some structures will fail the SIMPROF hypothesis
    # test. These are immediate blips.
    tryCatch(simprof.plot(test), error=function(e) {
      print("Blip Detected")
      run <<- TRUE
      
      # We will also plot the dendrogram, without SIMPROF coloration.
      # Just a nice visualisation
      dev.off()
      clust.data<-(hclust(vegdist(matrix.temp), method = "average"))
      par(mar = c(3,3,3,3))
      plot(clust.data, hang=-1,
           main = "Blip Dendrogram", xlab = "")
    })
    
    # This will activate the BLIP module and assign the correct 
    # category
    if(run){
      length <- 1
      start <- novel_bin
      end <- novel_bin
      length.bins <- 1
      Class <- "BLIP"
    }
    
    
    # Start of the length calculation module
    if(length == 0){
      
      length.output<-length.calculator.V2(test,novel_bin)
      
      length <- length.output$length 
      start <- length.output$start
      end <- length.output$end
      
      bin.start <- which(rev(number.names) == start)
      bin.end <- which(rev(number.names) == end)
      length.bins <- bin.start-bin.end
      
      
    }
    # Initialize a class variable
    
    
    # Allocate blips needs to be done!!!
    if(start == end | length.bins == 1){
      
      Class <- "BLIP"
      
      
    }
    
    # Test for full persistence
    if(Class != "BLIP"){
      
      Class <- "Persister"
      
    }
    
    if(novel_bin == number.names[length(number.names)]){
      Class <- "END"
      length <- 0
    }
    
    
    
    
    # Return data
    return.data <- list(c(list("ID" = ID, "begin" =start, 
                               "end" =end, "length"= length, 
                               "Class" = Class,
                               "length.bins" = length.bins)))
    names(return.data) <- full.ID
    return(return.data)
  })
  return(return.dat)
  
}


# V.7 includes a consistency test between NDF and SIMPROF

nov.cluster.id.V7 <- function(matrix, method_clus, 
                              alpha_clust){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][1]
  full.ID <- names(matrix)
  matrix <- matrix[[1]]
  
  number.names <- as.numeric(rownames(matrix))
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  # Find the average dissimilarities and matrix shapes
  mean_cum <- median(label_frame$raw.min.dist, na.rm =T)
  mean_inst <- median(label_frame$seq.dist, na.rm =T)
  nspp <- ncol(matrix)
  nyears <- nrow(matrix)
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  
  
  # Obtain novel index from the data
  
  novel_frame <- matrix %>% 
    dplyr::filter(str_detect(rownames(matrix), "novel"))
  
  if(nrow(novel_frame) > 1){
    iter.list <- c(1,2)
    print('double event')
  }else{
    iter.list <- 1
  }
  
  
  return.dat <- lapply(iter.list, function(x){
    # Find year and row data
    novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
    novel_frame <- novel_frame[x,]    
    
    novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
    
    novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
    
    # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
    matrix.temp <- matrix
    rownames(matrix.temp) <- rev(number.names)
    
    # Sometimes using bray-curtis/czekanowski can lead to an error where there are a 
    # large number of zeroes in the data. We address this by ofsetting zeroes by 
    # a small number.
    
    test <<- tryCatch(
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = T,
              num.simulated = 999, method.distance =vegdist, 
              method.cluster = method_clus, alpha=alpha_clust), 
      error=function(e) {
        print('ofsetting zeroes by small number')
        
        tryCatch(simprof(data = non_zero_offsetter(matrix.temp), num.expected = 1000, undef.zero = T,
                         num.simulated = 999, method.distance =vegdist, 
                         method.cluster = method_clus, alpha=alpha_clust),
                 error=function(e){
                   NULL
                 })
      })
    
    # Initialize parameters for the length calculations
    # Initialize a class variable
    
    Class <- "NONE"
    calc_length = T
    # Plot dendogram result
    par(mfrow = c(1,1))
    
    # If there is no support for significant clusters. We can not plot the dendrogram or 
    # do any meaningful length calculations.
    
    return.data = tryCatch(simprof.plot(test),
                           error = function(e){
                             print(paste0("No non-random structure in data for alpha = ", alpha_clust))
                             data_error <- list(c(list('consistency' = FALSE,
                                                       'median_seq_dis' = mean_inst,
                                                       'median_cum_dis' = mean_cum,
                                                       'nspp' =nspp,
                                                       'nbins'= nyears)))
                             calc_length <<- F
                             names(data_error) <- full.ID
                             return(data_error)
                           })
    
    
    if(calc_length){
      # Start of the length calculation module
      
      
      length.output<-length.calculator.V2(test,novel_bin)
      
      length <- length.output$length 
      start <- length.output$start
      end <- length.output$end
      
      bin.start <- which(rev(number.names) == start)
      bin.end <- which(rev(number.names) == end)
      length.bins <- bin.start-bin.end
      
      
      
      # Initialize a class variable
      
      
      # Allocate blips needs to be done!!!
      if(start == end | length.bins == 1){
        
        Class <- "BLIP"
        length <- 0
        
        
      }
      
      # Test for full persistence
      if(Class != "BLIP"){
        
        Class <- "Persister"
        
      }
      
      if(novel_bin == number.names[length(number.names)]){
        Class <- "END"
        length <- 0
      }
      
      
      ### Consistency module - is the novel community the start of a new cluster
      clusters <- test$significantclusters
      
      # Find the novel cluster and confirm it is the head of its own cluster.
      consistency <- FALSE
      for (list_index in 1:length(clusters)){
        if (novel_bin %in% as.numeric(clusters[[list_index]]) & 
            novel_bin == max(as.numeric(clusters[[list_index]]))){
          consistency <- TRUE
        }
      }
      
      ### Module for confirming to which cluster the post-novel communities 
      # belong, but only possible if the novel comp ends before the ts.
      
      if(end != number.names[length(number.names)]){
        cluster_vector = length.output$clustervec
        postNov = cluster_vector[which(names(cluster_vector) == end)]
        preNov_clus = cluster_vector[1:(which(names(cluster_vector) == novel_bin)-1)]
        if(postNov %in% preNov_clus){
          FirstPostNov = "Pre_Novel_State"
        }else{
          FirstPostNov = "New_exploratory_State"
        }
      }else{
        FirstPostNov = NA
      }
      
      # Return data in df
      return.data <- list(c(list("ID" = ID, "begin" =start, 
                                 "end" =end, "length"= length, 
                                 "Class" = Class,
                                 "length.bins" = length.bins,
                                 'median_seq_dis' = mean_inst,
                                 'median_cum_dis' = mean_cum,
                                 'nspp' =nspp,
                                 'nbins'= nyears,
                                 'consistency' = consistency,
                                 'FirstPostNov' = FirstPostNov)))
      names(return.data) <- full.ID
    }
    return(return.data)
  })
  return(return.dat)
  
}



















