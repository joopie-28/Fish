### ANOSIM PERSISTENCE FRAMEWORK ####

# J.M. Sassen 14-06-2022

# Script contains all necessary functions to run the anosim detection algorithm
# as well as functions for simulating matrices

# Testing the ANOSIM persistence method on
# a few simulated time series scenario's

# First create the simulated scenario's;
# these represent cases where the correct
# answer is clear. Use these two simple
# functions to simulate any scenario
# you would like!


# Simulate a one-off novelty emergence followed by
# an immediate return to 'normality'

matrix.simulator <- function(rows, sim.type){
  
  # Simulate a starting condition 
  # Poission picking to set initial state
  vec <- rpois(10, c(runif(10, 0, 500)))
  
  iniState <- vec/sum(vec)
  
  # Specify matrix size (fixed at 10 species)
  matrix.sim <- matrix(nrow = rows, 
                       ncol = 10)
  
  # Different simulation types!
  
  # Stable-state
  if (sim.type == "stable"){
    for (i in 1:nrow(matrix.sim)){
      news <- sapply(iniState, function(x){
        
        temp <- runif(1, x-0.02, x + 0.02)
        if(temp < 0){
          temp <- 0
        }
        return(temp)
      })
      matrix.sim[i,] <- news/sum(news)
    }
  }
  
  # One-off blips
  if(sim.type == "blip"){
    
    loc <- sample(c(6:(rows-5)), 1)
    
    for (i in 1:nrow(matrix.sim)){
      
      # set a random number for the nov community
      # not too far too each end.
      
      if(i > loc | i < loc){
        var <- 0.02
        
      }
      
      else{
        var <- 0.9
        
      }
      
      news <- sapply(iniState, function(x){
        
        temp <- runif(1, x-var, x + var)
        
        if (temp < 0){
          temp <- 0
        }
        
        return(temp)
        
      })
      
      
      matrix.sim[i,] <- news/sum(news)
      
    }
    
  }
  
  # Short persistence 
  if(sim.type == "short.pers"){
    
    loc <- sample(c(6:(rows-5)), 1)
    
    # Make it harder for the algorithm
    # by changing the persistence length
    # every iteration.
    
    pers.length <- round(runif(1, 2, 10), digits = 0)
    print(pers.length)
    for (i in 1:nrow(matrix.sim)){
      
      # set a random number for the nov community
      # not too far too each end.
      
      if(i < loc | i >= loc + pers.length){
        origin <- iniState
        var <- 0.02
      }
      
      if(i == loc){
        origin <- iniState
        var = 0.9
      }
      
      if(i %in% c((loc+1):(loc+(pers.length-1)))){
        origin <- matrix.sim[loc,]
        var <- 0.02
      }
      
      news <- sapply(origin, function(x){
        
        temp <- runif(1, x-var, x + var)
        
        if (temp < 0){
          temp <- 0
        }
        
        return(temp)
        
      })
      
      
      matrix.sim[i,] <- news/sum(news)
      
    }
    
  }
  
  # Full-persistence (no reversion)
  if(sim.type == "full.pers"){
    
    # This sets the location for emergence,
    # randomly
    loc <- sample(c(6:(rows-5)), 1)
    
    # Make it harder for the algorithm
    # by changing the persistence length
    # every iteration.
    
    for (i in 1:nrow(matrix.sim)){
      
      # Here, there is a low variability state prior to novelty
      if(i < loc){
        origin <- iniState
        var <- 0.02
      }
      
      if(i == loc){
        origin <- iniState
        var = 0.9
      }
      # Stays close to novelty towards the end
      if(i %in% c((loc+1):nrow(matrix.sim))){
        origin <- matrix.sim[loc,]
        var <- 0.02
      }
      
      news <- sapply(origin, function(x){
        
        temp <- runif(1, x-var, x + var)
        
        if (temp < 0){
          temp <- 0
        }
        
        return(temp)
        
      })
      
      
      matrix.sim[i,] <- news/sum(news)
      
    }
    
  }
  
  # Slow Decay (this is the most challenging one,
  # calculating length will be somewhat arbitrary!)
  if(sim.type == "slow.decay"){
    
    loc <- sample(c(6:(rows-5)), 1)
    print(paste0("loc", loc))
    for (i in nrow(matrix.sim):1){
      
      
      # Now we depend on the previous state
      # and slowly move away. Set var a little
      # bit higher to allow for real shifts
      if(i < loc){
        origin <- matrix.sim[i+1,]
        var <- 0.1
      }
      # set a random number for the nov community
      # not too far too each end.
      
      if(i > loc){
        origin <- iniState
        var <- 0.02
      }
      
      # Novel emergence
      if(i == loc){
        origin <- iniState
        var = 0.9
      }
      
      news <- sapply(origin, function(x){
        
        temp <- runif(1, x-var, x + var)
        
        if (temp < 0){
          temp <- 0
        }
        
        return(temp)
        
      })
      
      
      matrix.sim[i,] <- news/sum(news)
      
    }
    
  }
  
  
  
  # Convert to proper format
  matrix.df <- as.data.frame(matrix.sim)
  
  # Visualize the simulated scenario
  # and check simulation success
  
  test <- identify.novel.gam.MDS(matrix.df,site = "sim1", metric = "bray", plot =T, plot.data = FALSE, gam.max.k = -1, alpha = 0.05)
  if(!any(test$cat[6:rows] == "novel")){
    # To ensure novelty is actually present!
    print("Oops, recalibrating")
    matrix.df <- NULL
  }
  
  return(matrix.df)
}

# Depends on the ^simulator, allows you to select
# simulation type and number of iterations
generate.matrices <- function(sim.type, quant){
  matrices <- list()
  i <- 1
  while (i <= quant){
    
    temp <- matrix.simulator(30, sim.type)
    
    if (is.null(temp)){
      i <- i
    }
    else{
      matrices[[i]] <- temp
      i <- i + 1
    }
    
    
  }
  
  for (i in 1:length(matrices)){
    
    names(matrices)[i] <- paste0("sim", i)
  }
  
  return(matrices)
}


# The ultimate ANOSIM algorithm!

# Finds optimal R and returns persistence length of novel community
# Works perfectly on simulated data

novel.length.algo <- function(novel.output){
  
  # Apply this large chunk to a list of matrices  
  grand_output <- lapply(1:length(novel.output), function(i){  
    
    matrix <- novel.output[[i]]
    
    # Print the name so we know what we're on
    name <- names(novel.output)[i]
    print(name)
    
    # Create a label frame by running Novel Detection Framework
    label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                          alpha = 0.05,
                                          metric = "bray",
                                          plot =F, 
                                          site = name,
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
    novel_bin <- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
    
    novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
    
    
    # Quick fail safe for novelty at the very end of a time series,
    # The output is simply 1! We also have no idea what the post-trajectory will be
    
    if(novel_index == 1){
      return(list("Length_years" = NA, "Class" = "END"))
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
    
    #### Primary Computation Module ####
    
    # Set up some vectors for our ANOSIM iterations
    r.vector <- NULL
    sig.vector <- NULL
    timebin <- NULL
    
    # Assign a temporary matrix for the loop
    matrix.temp <- matrix
    
    # Store the total lenth of the entire time series
    total_timeseries_length <- nrow(matrix.temp)
    
    # Loop over multiple configurations to find the most suitable
    print("Scanning Compositional Space")
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
    best.R <- ifelse(any(ANOSIM.df$Sig.p <= 0.05), max(subset(ANOSIM.df, Sig.p <= 0.05)$R), 0)
    
    # Extract length from associated R statistic, if none suitable, set length to two time bin
    # i.e. the novel plus one preceding! Note that we need to set 1 to 2 because even when
    # novelty is a blip, qw have a length. 
    length <- ifelse(best.R > 0.75, max(ANOSIM.df[which(ANOSIM.df$R == best.R), "Set"],2), 2)
    
    #### Classification module ####
    
    # Convert time in bins to time in years (not every bin has a time step of 1 year)
    state.end.bin <- as.numeric(strsplit(rownames(matrix[(novel_index-length+1), ]), split = "-")[[1]][2])
    actual.length <- novel_bin - state.end.bin
    
    # Short persistence if decay and full agree, and our threshold is validated
    if(best.R > 0.75){
      message <- "SHORT.PERSISTENCE"
      
      # Blips are special cases of short persisters where a temporary exploration of
      # new compositional space is flanked by two similar states.
      if(length == 2){
        message <- "BLIP"
        
      }
      
      # Full persistence will have a high R and a length equal to the 
      # remainder of the series
      if(length == novel_index){
        message <- "FULL.PERSISTENCE"
      }
    }else{
      
      # This is for cases that lack an R stat > 0.75. These cases reflect a lack of difference between
      # pre and post novel communities, which means the novel community is just an extreme outlier in 
      # an ongoing dynamic.
      message <- "BLIP"
      
    }
    
    # Format output effectively
    output <- list(name ,ANOSIM.df, best.R, length, actual.length, novel_bin, novel_index, total_timeseries_length,message)
    names(output) <- c("Simulation" ,"data", "R", "Length_steps", "Length_years", "Year_of_Emergence", "Emergence_index", "Timeseries_length", "Class")
    
    
    return(output)
    
  })
  
  return(grand_output)
  
}
