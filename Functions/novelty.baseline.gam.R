#########################################################
#### Novelty Transience and change point analysis #######
#########################################################

####################
### Results ########
####################

# Analysis of Similarity is a test used to assess two groups are truly different, 
# specifically tailored to community matrix data.

# We first need to assign pre and post novelty groups.

novelty.trajectory.plotter <- function(ID){
  
  matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
  
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[ID]]
    
  }
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                    alpha = 0.05,
                                    metric = "bray",
                                    plot = FALSE, 
                                    site = ID,
                                    plot.data = FALSE,
                                    gam.max.k = -1)

  if(all(label_frame$cat[5:nrow(label_frame)] != "novel")){
    return("No novelty")
  }
  
  # Find the novel community index

  for (i in 1:5) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
  
  }

  # Quick loop to assign names

  for (i in 6:dim(matrix)[1]) {
    for (j in 6:dim(label_frame)[1]) {
    
      if ((rownames(matrix)[i]) == (label_frame$bins[j])){
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", rownames(matrix)[i])
      }
    
    } 
 } 


  # Obtain novel index

  novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
  novel_frame <- novel_frame[1,]
  
  if(nrow(novel_frame) == 0){
  
    novel_index <- nrow(matrix)+1
  } 
  if(nrow(novel_frame) > 0){
  
    novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])

    novel_index <- which(label_frame$bins == novel_bin)
  }
  
  # Now assign two groups based on that index

  # Also return normal names
  
  rownames(matrix) <- label_frame$bins
  
  matrix$stage <- NULL
  
  for (i in 1:nrow(matrix)) {
  
    if(i < novel_index){
      matrix$stage[i] <- "Pre_Novel"
    }
    if (i >= novel_index) {
      matrix$stage[i] <- "Post_Novel"
    }
  
  }

  if(length(unique(matrix[-(novel_index),]$stage)) > 1 ){
    # Run the analysis on post and pre only, exclude the novel community.
    similarity_test <- anosim(matrix[-(novel_index),-(ncol(matrix))], grouping = matrix[-(novel_index),]$stage, permutations = 1000, distance = "bray")
  } else{
    similarity_test <- "No Novelty"
  }
  
  matrix <- matrix[,-(ncol(matrix))]
  
  
  plot <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                    alpha = 0.05,
                                    metric = "bray",
                                    plot = TRUE, 
                                    site = ID,
                                    plot.data = FALSE,
                                    gam.max.k = -1)
  
  
  output <- list(label_frame, similarity_test)
  
  names(output) <- c("Novelty parameters", "ANOSIM")
  
  return(output)
}

transient_novelty_GAM_function <- function(ID, metric){
  par(mfrow=c(1,1), mar = c(0,0,0,0))
  set.k = -1
  
  matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
  
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
    
  }
  
  # NOVEL DISSIMILARITY #####
  # Obtain distance matrix for all sites in timeseries
  site.dist <- as.data.frame(as.matrix(vegdist(matrix,
                                               method=metric)))
  
  # Lets label!
  
  if (nrow(matrix) >= 10) {
    if (ncol(matrix) >=5){
      
      label_frame <- identify.novel.gam(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = metric,
                                        plot = FALSE, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
      
    }
    else{
      return(NA)
    }
  } 
  
  # Going to remove the first 5 of the label frame because
  # these communities are unreliable. The effect of this is
  # then that we do not consider novelty in these first 5.
  # Assign the first 5 a "background" status, solves for the
  # wrongful recognition of novelty in the first 5.
  
  bins <- as.numeric(rownames(matrix))
  
  for (i in 1:5) {
    rownames(site.dist)[i] <- paste0("back-", rownames(site.dist)[i])
    colnames(site.dist)[i] <- paste0("back-", colnames(site.dist)[i])
    
  }
  
  # Quick loop to assign names
  
  for (i in 6:dim(site.dist)[1]) {
    for (j in 6:dim(label_frame)[1]) {
      
      if ((rownames(site.dist)[i]) == (label_frame$bins[j])){
        rownames(site.dist)[i] <- paste0(label_frame$cat[j], "-", rownames(site.dist)[i])
      }
      
      if((colnames(site.dist)[i]) == (label_frame$bins[j])){
        colnames(site.dist)[i] <- paste0(label_frame$cat[j], "-", colnames(site.dist)[i])
      }
    } 
  } 
  
  
  # Calculate difference between novel and selected states
  
  novel_frame <- site.dist %>% filter(str_detect(rownames(site.dist), "novel"))
  
  if(nrow(novel_frame) == 0){
    return("No novelty")
  }
  
  novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
  
  novel_index <- which(label_frame$bins == novel_bin)
  
  
  
  
  # Let's try from the first point
  
  novel.dist <- (c(NA,site.dist[1,-1]))
  
  novel.dist <- as.numeric(novel.dist)
  
  # Compare to all pre-novel states
  
  novel.dist <- sapply(1:dim(site.dist)[1], function(n){
    if(n==1){return(NA)}
    min(site.dist[n,1:(novel_index-1)][-n])
  })
  
  # Remove 0 and 1s for beta regression
  
  novel.dist.tr <- (novel.dist * (length(novel.dist)-1) + 0.5) / length(novel.dist)
  
  # model localised trend over time and extract residuals to get dissimilarity 
  # compared to local mean.
  
  if(var(novel.dist, na.rm=TRUE)==0){return(NULL)}
  
  if(!is.na(set.k)){
    novel.gam <- gam(novel.dist.tr ~ s(bins, bs="cr", k= set.k), 
                     family=betar(),
                     method="REML")
  } else{
    novel.gam <- gam(novel.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # COMPARING OBSERVED TO EXPECTED NOVEL BASELINE DISIMILARITY ####
  
  # This process calculates the p-value of the observed disimilarity score
  # being part of the expected distribution at the point in the time-series.
  
  # convert mu & phi beta parameters to A & B to use in qbeta.
  # mu = the mean predicted dissimilarity at each time point based on the 
  #      additive model.
  # phi = a dispersion parameter, which is set globally across the model.
  #       phi is unfortunately hidden deep in the gam model as a text string.
  
  novel.mu <- c(novel.gam$fitted.values)
  phi <- as.numeric(substr(novel.gam$family$family,
                           regexpr("\\(", novel.gam$family$family)+1,
                           nchar(novel.gam$family$family)-1))
  
  # shape parameters used in qbeta.
  A = novel.mu * phi
  B = phi - A
  
  # predict 5% and 95% prediction intervals from beta distribution parameters. 
  # We use 95% rather than 97.5% because this is a one-tailed probability test.
  # We are not concerned about communities that are MORE similar than predictions.
  # This is done for each bin along the time-series.
  novel.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               novel.p = pbeta(novel.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # Plot the new analysis
  pre_mean <- mean(novel.dist[2:(novel_index-1)])
  post_mean <- mean(novel.dist[(novel_index+1):length(novel.dist)])
  
  plot(novel.dist[-1] ~ bins[-1], type="n", axes=FALSE,
       ylim=c(min(novel.dist[-1], na.rm=TRUE), 
              max(novel.dist[-1], na.rm=TRUE)+0.1), ylab = "", xlab = "")
  polygon(x=c(bins[-1], rev(bins[-1])),
          y=c(novel.p[,1], rev(novel.p[,2])), 
          col="grey80", border=NA)
  lines(plogis(predict(novel.gam)) ~ bins[-1], col="grey50", lwd=2)
  lines(novel.dist[-1] ~ bins[-1])
  abline(v= novel_bin, lty = 2, lwd = 2, col = "orange")
  
  
  segments(min(bins), post_mean, novel_bin, post_mean, col = "red", lty = 2)
  
  segments(novel_bin, pre_mean, max(bins[-1]), pre_mean, col = "red", lty = 2)
  
  axis(side=1)
  mtext(side=1, "Time series age", line=2)
  
  axis(side=2)
  mtext(side=2, text = "Pre-Novelty\ndissimilarity", line=3, las=0)
  

  
  return(data.frame(bins, novel.dist.tr))
  
  
  
}

# Convert to a listing function

novelty.trajectory.lists <- function(check_list){
  nam <- lapply(check_list, 
                function(TimeSeries_ID){
                  print(TimeSeries_ID)
                  
                  temp <- novelty.trajectory.plotter(TimeSeries_ID)
                  
                  return(temp)
                })
  return(nam)
}

cut.off.generator <- function(anosim.lists, number = -1){
  
  anosim.plots.list <- anosim.lists
  
  dis.vector <- vector(length = length(anosim.plots.list))
  
  # Retrieve all dissimilarities
  for (i in 1:length(anosim.plots.list)){
    if(typeof(anosim.plots.list[[i]]$ANOSIM) == "character") {
      # Set to 0 if we cannot perform novelty test
      dis <- 0
      
    }
    else{
      dis <- anosim.plots.list[[i]]$ANOSIM$statistic
      
    }
    dis.vector[i] <- dis
    
  }
  
  novel_indices <- which(dis.vector > number)
  
  true.novelty.output <- anosim.plots.list[novel_indices]
  
  # Add names
  
  for (i in 1:length(true.novelty.output)){
    names(true.novelty.output)[i] <- paste0(true.novelty.output[[i]]$`Novelty parameters`$site[1], ",", 
                                            round(true.novelty.output[[i]]$ANOSIM$statistic, digits = 2))
  }
  
  return(true.novelty.output)
}

# Finds optimal R and returns persistence of novel community
novel.length.checker <- function(novel.output, cut.off){
  
  grand_output <- list()
  run <- FALSE
  for (i in 1:length(novel.output)){
    
    # Extract R statistic and Time Series ID from list names
    statistic <- as.numeric(strsplit(names(novel.output[i]), "-")[[1]][2])
    significance <- novel.output[[i]]$ANOSIM$signif
    
    # Save the pre-optimization R statistic
    statistic.0 <- statistic
    significance.0 <- significance
    
    ID <- strsplit(names(novel.output[i]), "-")[[1]][1]
    print(ID)
    # Now use a modified form of trajectory.plotter to investigate 
    # R statistics of varying configurations.
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
    
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[ID]]
      
    }
    
    
    label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                          alpha = 0.05,
                                          metric = "bray",
                                          plot = FALSE, 
                                          site = ID,
                                          plot.data = FALSE,
                                          gam.max.k = -1)
    
    if(all(label_frame$cat[5:nrow(label_frame)] != "novel")){
      return("No novelty")
    }
    
    # Find the novel community index
    
    for (i in 1:5) {
      rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
      
    }
    
    # Quick loop to assign names
    
    for (i in 6:dim(matrix)[1]) {
      for (j in 6:dim(label_frame)[1]) {
        
        if ((rownames(matrix)[i]) == (label_frame$bins[j])){
          rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", rownames(matrix)[i])
        }
        
      } 
    } 
    
    
    # Obtain novel index
    
    novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
    novel_frame <- novel_frame[1,]
    
    if(nrow(novel_frame) == 0){
      
      novel_index <- nrow(matrix)+1
    } 
    if(nrow(novel_frame) > 0){
      
      novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
      
      novel_index <- which(label_frame$bins == novel_bin)
    }
    
    # Now assign two groups based on that index
    
    # Also return normal names
    
    rownames(matrix) <- label_frame$bins
    
    matrix$stage <- NULL
    
    for (i in 1:nrow(matrix)) {
      
      if(i < novel_index){
        matrix$stage[i] <- "Pre_Novel"
      }
      if (i >= novel_index) {
        matrix$stage[i] <- "Post_Novel"
      }
    }
    
    # Just make sure we store the original matrix
    matrix_1 <- matrix
    
    # Exclude post-novel communities 1 by 1
  
    for (k in nrow(matrix):(novel_index + 1)){
      
      # Just make sure we store the original matrix 
    
      matrix <- matrix[-k, ]
    
      if(length(unique(matrix[-(novel_index),]$stage)) > 1 ){
      
        
       # Run the analysis on post and pre only, exclude the novel community.
        similarity_test <- anosim(matrix[-(novel_index),-(ncol(matrix))], grouping = matrix[-(novel_index),]$stage, permutations = 1000, distance = "bray")
        
  
      # Continuously compare the similarity, once we have reached optimum we stop the loop
        if(similarity_test$statistic <= statistic + 0.05 | similarity_test$signif > 0.05) {
          run <- TRUE
          break
        }else{
          
          statistic <- similarity_test$statistic
          significance <- similarity_test$signif
          
        }
      }
      else{
        # In the case that novelty is far to the end
        novel_length <- as.numeric(rownames(matrix_1)[novel_index]) - as.numeric(rownames(matrix_1)[nrow(matrix_1)])
        output <- list(ID, statistic.0, statistic, novel_length, significance.0 ,significance)
        names(output) <- c("ID", "old_R", "new_R", "length", "old_significance", "new_significance")
      }
    }
    
    # Calculate length of novelty
    if (run){
      optimal_index <- k 
      novel_length <- as.numeric(rownames(matrix_1)[novel_index]) - as.numeric(rownames(matrix_1)[optimal_index])
      significance <- significance
      matrix <- matrix_1
    
      output <- list(ID, statistic.0, statistic, novel_length, significance.0 ,significance)
      names(output) <- c("ID", "old_R", "new_R", "length", "old_significance", "new_significance")
    }
    run <- FALSE
    grand_output <- c(grand_output, list(output))
  }
  for (i in 1:length(grand_output)){
    names(grand_output)[i] <- grand_output[[i]]$ID
  }
  
  # Filter out what we want based on cut.off value
  
  for (i in 1:length(grand_output)){
    
    if(grand_output[[i]]$new_R < cut.off){
      grand_output[[i]] <- NA
    }
    
  }
  
  grand.output.fil <- grand_output[is.na(grand_output) == FALSE]
  
  return(grand.output.fil)
}






#### Novel Length Algo as per 13-06. This is the more complete version
### that involves decay matrices etc. Deemed a bit too complicated for
### the nov.comm paper!

novel.length.depr <- function(novel.output, cut.off){
  
  
  
  grand_output <- lapply(1:length(novel.output), function(i){  
    
    matrix <- novel.output[[i]]
    
    # Create a label frame by running Novel Detection Framework
    label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                          alpha = 0.05,
                                          metric = "bray",
                                          plot =T, 
                                          site = names(novel.output)[i],
                                          plot.data = FALSE,
                                          gam.max.k = -1)
    
    # Print the name so we know what we're on
    name <- names(novel.output)[i]
    print(name)
    
    # Module for proper orientation, can be conditional
    
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
          rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", rownames(matrix)[i])
        }
        
      } 
    } 
    
    
    # Obtain novel index from the data
    
    novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
    novel_frame <- novel_frame[1,]
    
    
    if(nrow(novel_frame) > 0){
      
      novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
      
      novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
      
      
    }
    
    
    # Quick fail safe for novelty at the very end of a time series,
    # The output is simply 1! We also have no idea what the post-trajectory will be
    
    if(novel_index == 1){
      return(list(length = 1, message = "Novelty at very end"))
    }
    
    # Now assign two groups based on that index
    
    matrix$stage <- NULL
    
    for (i in 1:nrow(matrix)) {
      
      if(i <= novel_index){
        matrix$stage[i] <- "Post_Novel"
      }
      if (i > novel_index) {
        matrix$stage[i] <- "Pre_Novel"
      }
    }
    
    #### Primary Module ####
    
    # Set up some vectors for our ANOSIM iterations
    r.vector <- NULL
    sig.vector <- NULL
    timebin <- NULL
    
    # Assign a temporary matrix for the loop
    matrix.temp <- matrix
    
    # Loop over multiple configurations to find the most suitable
    print("Scanning compositional space")
    for (k in 1:(novel_index)){
      message <- NULL
      
      if(length(unique(matrix[-(novel_index),]$stage)) > 1 ){
        
        # Run the analysis on post and pre only, exclude the novel community.
        similarity_test <- anosim(matrix.temp[,-(ncol(matrix.temp))], 
                                  grouping = matrix.temp$stage, 
                                  permutations = 1000, 
                                  distance = "bray")
        
        # Store relevant results
        r.vector[k] <- (similarity_test$statistic)
        sig.vector[k] <- (similarity_test$signif)
        timebin[k] <- nrow(subset(matrix.temp, stage == "Post_Novel"))
        
        # This "restores" the matrix so that we can iterate over it again
        matrix.temp<-matrix
        matrix.temp[1:(novel_index-k), ]$stage <- "Pre_Novel"
      }
      
    }
    # Create a tidy data frame for our iteration results
    ANOSIM.df <- data.frame("R" = r.vector, "Sig.p"=sig.vector, "Set"=timebin)
    
    # Set Significance cut off a bit lower because our data is limited..
    best.R <- ifelse(any(ANOSIM.df$Sig.p < 0.15), max(subset(ANOSIM.df, Sig.p <= 0.15)$R), 0)
    # If we can not get any significance at all, we set length to 0 and call them 'oscillators'
    length <- ifelse(best.R > 0, ANOSIM.df[which(ANOSIM.df$R == best.R), "Set"], 0)
    
    
    
    
    
    
    #### Decayer Module ####
    
    # I would like to add in an extra module that then breaks up the post-novel states
    # to test if they are moving away from the origin. Helps differentiate structured
    # movement back and forth from an exploration of novel compositional space.
    
    # However, we can only do this if the number of post-novel communities is substantial
    # at least 3 or so..
    
    print("Inspecting decay trajectory")
    
    # Extract the 'decay matrix', the communities post novelty
    # we will scan these to see if there's distinct clusters
    decay.matrix <- subset(matrix, stage == "Post_Novel")
    
    # Our criterium for running the decay analysis,
    # ANOSIM wont't work otherwise and getting
    # significant results is hard
    if(nrow(decay.matrix) >4){
      
      # Set up new result vectors
      decay.R <- NULL
      decay.Sig <- NULL
      decay.length <- NULL
      
      # Run the loop, similar to initial procedure
      for(i in (nrow(decay.matrix)-1):1){
        
        decay.matrix$stage <- "novel"
        decay.matrix$stage[c(1:i)] <- "different"
        
        
        # Find the allocation with max distance between groups
        
        test <- anosim(decay.matrix[,-(ncol(decay.matrix))], 
                       grouping = decay.matrix$stage, 
                       permutations = 1000, distance = "bray")
        
        decay.R[i] <- test$statistic
        decay.Sig[i] <- test$signif
        decay.length[i] <- nrow(subset(decay.matrix, stage == "novel"))
        
      }
      # Store these neatly
      DECAY.df <- data.frame("R" = decay.R, "Sig.p"=decay.Sig, "Set"=decay.length)
      decay.best <- DECAY.df[which.max(DECAY.df$R), "Set"] 
      decay.R <- DECAY.df[which.max(DECAY.df$R), "R"] 
      
      #### Pre and Post Comparison Module ####
      # One final test: are the communities surrounding the novel emergence significantly different?
      
      # Initialize a new matrix
      diff.matrix <- matrix
      
      # Remove the novel community such that comparison is just between pre and post
      diff.matrix <- diff.matrix[-(novel_index), ]
      
      # Loop through each module, iteratively removing communities as to find
      # the true pattern and max R. 
      
      max.iteration <- nrow(subset(matrix, stage == "Post_Novel"))
      
      diff.R <- NULL
      diff.Sig <- NULL
      diff.length <- NULL
      # Removing states to see if we can find significant differences
      for (i in 1:(max.iteration-2) ){
        
        # Remove the first row each time unitl we hit the target
        diff.matrix <- diff.matrix[-1, ]
        print(diff.matrix)
        
        diff.test <- anosim(diff.matrix[,-(ncol(diff.matrix))], 
                            grouping = diff.matrix$stage, 
                            permutations = 1000, distance = "bray")
        
        diff.R[i] <- diff.test$statistic
        diff.Sig[i] <- diff.test$signif
        diff.length[i] <- nrow(subset(diff.matrix, stage == "Post_Novel"))
      }
      
      diff.df <- data.frame("R" = diff.R, "Sig.p"=diff.Sig, "Length"=diff.length)
      
      # If there's no difference, it can only be a blip.
      
      sig.df <- subset(diff.df, "Sig.p" <= 0.05)
      
      if(nrow(sig.df) == 0 | max(sig.df$R) < 0.75){
        message <- "BLIP.D"
      }
      
      # We have enough to analyse decay so we will activate the decayer option
      run.decayer <- TRUE
      
      
    }else{
      
      # Provide some empty values in the case we cannot estimate decay df
      DECAY.df <- NA
      decay.best <- NA
      decay.R <- NA
      run.decayer <- FALSE
    }
    # Short persistence if decay and full agree, and our threshold is validated
    if(best.R > 0.75){
      message <- "SHORT.PERSISTENCE"
      
      # Blips are special cases of short persisters
      if(length == 1){
        message <- "BLIP"
        
        # Correct for time series with novelty at the very end
        if(novel_index == nrow(matrix)){
          message <- "SINGLE_END"
        }
      }
    }
    
    # Full persistence will have a high R and a length equal to the 
    # remainder of the series
    if(best.R > 0.75 & length == novel_index){
      message <- "FULL.PERSISTENCE"
      
    }
    
    # Finally, we need something for the edge cases, those that 
    # can not really be put in a box. only running if we have a decay estimate
    if(run.decayer){
      if(decay.R > 0.75 & best.R < 0.75){
        message <- "NOVEL.DECAYER"
        # Use decay length as this is more reasonable
        # in these cases
        if(DECAY.df[which.max(DECAY.df$R), "Sig.p"]<0.05 ){
          length <- decay.best
        }
        
        
      }
    }
    
    
    # Placeholder, I want to extend this module by comparing states
    if(length == 0){
      message <- "OSCILLATOR"
      # Set to 1 because it always last for min 1 year!
      length <- 1
    }
    
    output <- list(name ,ANOSIM.df, best.R, length, novel_index, DECAY.df, message, diff.R)
    names(output) <- c("Simulation" ,"data", "R", "Length", "Nov.emergence", "Decay.stats", "Class", "Pre-post")
    
    
    return(output)
    #}
  })
  return(grand_output)
  
}






# Sensitivity analysis of ANOSIM

ANOSIM.sensitivity <- function(ID){
  
  matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
  
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[ID]]
    
  }
  if(is.null(matrix)){
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[ID]]
    
  }
  
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot = FALSE, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  
  
  
  
  
  
  
  # Randomly assign a novel index. however, we will try to keep this a little in the middle so
  # that the comparison makes sense.
  
  rows <- nrow(matrix)
  middle <- ceiling(rows/2)
  novel_index <- sample(seq(middle-3,middle+3 ), size = 1)

  # Now assign two groups based on that index
  
  # Also return normal names
  
  rownames(matrix) <- label_frame$bins
  
  matrix$stage <- NULL
  
  for (i in 1:nrow(matrix)) {
    
    if(i < novel_index){
      matrix$stage[i] <- "Pre_Novel"
    }
    if (i >= novel_index) {
      matrix$stage[i] <- "Post_Novel"
    }
    
  }
  
  if(length(unique(matrix[-(novel_index),]$stage)) > 1 ){
    # Run the analysis on post and pre only, exclude the novel community.
    similarity_test <- anosim(matrix[-(novel_index),-(ncol(matrix))], grouping = matrix[-(novel_index),]$stage, permutations = 1000, distance = "bray")
  } else{
    similarity_test <- "No Novelty"
  }
  
  matrix <- matrix[,-(ncol(matrix))]
  
  
  plot <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                 alpha = 0.05,
                                 metric = "bray",
                                 plot = TRUE, 
                                 site = ID,
                                 plot.data = FALSE,
                                 gam.max.k = -1)
  
  
  output <- list(label_frame, similarity_test)
  
  names(output) <- c("Novelty parameters", "ANOSIM")
  
  return(output)
}

ANOSIM.sensitivity.lists <- function(check_list){
  nam <- lapply(check_list, 
                function(TimeSeries_ID){
                  print(TimeSeries_ID)
                  
                  temp <- ANOSIM.sensitivity(TimeSeries_ID)
                  
                  return(temp)
                })
  
  output <- cut.off.generator(nam)
  
  
  return(output)
}

ANOSIM.checker <- function(novel.output, cut.off){
  
  grand_output <- list()
  run <- FALSE
  for (i in 1:length(novel.output)){
    
    # Extract R statistic and Time Series ID from list names
    statistic <- as.numeric(strsplit(names(novel.output[i]), "-")[[1]][2])
    significance <- novel.output[[i]]$ANOSIM$signif
    
    # Save the pre-optimization R statistic
    significance.0 <- significance
    statistic.0 <- statistic
    
    ID <- strsplit(names(novel.output[i]), "-")[[1]][1]
    print(ID)
    # Now use a modified form of trajectory.plotter to investigate 
    # R statistics of varying configurations.
    
    matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[ID]]
    
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$nearctic_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$afrotropics_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$neotropics_mat_A[[ID]]
      
    }
    if(is.null(matrix)){
      
      matrix <- Fish_Communities_A$BioRealm_Matrices_A$australasia_mat_A[[ID]]
      
    }
    
    
    label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                          alpha = 0.05,
                                          metric = "bray",
                                          plot = FALSE, 
                                          site = ID,
                                          plot.data = FALSE,
                                          gam.max.k = -1)
    
    
    
    # Randomly assign a novel index. however, we will try to keep this a little in the middle so
    # that the comparison makes sense.
    
    rows <- nrow(matrix)
    middle <- ceiling(rows/2)
    novel_index <- sample(seq(middle-3,middle+3 ), size = 1)
    
    # Now assign two groups based on that index

    
    # Also return normal names
    
    rownames(matrix) <- label_frame$bins
    
    matrix$stage <- NULL
    
    for (i in 1:nrow(matrix)) {
      
      if(i < novel_index){
        matrix$stage[i] <- "Pre_Novel"
      }
      if (i >= novel_index) {
        matrix$stage[i] <- "Post_Novel"
      }
    }
    
    # Just make sure we store the original matrix
    matrix_1 <- matrix
    
    # Exclude post-novel communities 1 by 1
    
    for (k in nrow(matrix):(novel_index + 1)){
      
      # Just make sure we store the original matrix 
      
      matrix <- matrix[-k, ]
      
      if(length(unique(matrix[-(novel_index),]$stage)) > 1 ){
        
        
        # Run the analysis on post and pre only, exclude the novel community.
        similarity_test <- anosim(matrix[-(novel_index),-(ncol(matrix))], grouping = matrix[-(novel_index),]$stage, permutations = 1000, distance = "bray")
        
        # Continuously compare the similarity, once we have reached optimum we stop the loop
        if(similarity_test$statistic <= statistic + 0.05 | similarity_test$signif > 0.05) {
          run <- TRUE
          break
        }else{
          
          statistic <- similarity_test$statistic
          significance <- similarity_test$signif
        }
      }
      else{
        # In the case that novelty is far to the end
        novel_length <- as.numeric(rownames(matrix_1)[novel_index]) - as.numeric(rownames(matrix_1)[nrow(matrix_1)])
        output <- list(ID, statistic.0, statistic, novel_length, significance.0 ,significance)
        names(output) <- c("ID", "old_R", "new_R", "length", "old_significance", "new_significance")
      }
    }
    
    # Calculate length of novelty
    if (run){
      optimal_index <- k 
      novel_length <- as.numeric(rownames(matrix_1)[novel_index]) - as.numeric(rownames(matrix_1)[optimal_index])
      significance <- significance
      matrix <- matrix_1
      
      output <- list(ID, statistic.0, statistic, novel_length, significance.0 ,significance)
      names(output) <- c("ID", "old_R", "new_R", "length", "old_significance", "new_significance")
    }
    run <- FALSE
    grand_output <- c(grand_output, list(output))
  }
  for (i in 1:length(grand_output)){
    names(grand_output)[i] <- grand_output[[i]]$ID
  }
  
  # Filter out what we want based on cut.off value
  
  for (i in 1:length(grand_output)){
    
    if(grand_output[[i]]$new_R < cut.off){
      grand_output[[i]] <- NA
    }
    
  }
  
  grand.output.fil <- grand_output[is.na(grand_output) == FALSE]

  return(grand.output.fil)
}



