# # This functions extract information from the novelty results and returns a data frame with summary
# statistics.

novel.comm.analyzer <- function(full.novel.mat){
  
  # Create a frame with only novel communities for further analysis
  nov.comm.all <- subset(full.novel.mat, cat == "novel")
  
  # Isolate non-blips
  nov.comm.no.blip <- subset(full.novel.mat, cat == "novel" & novel.length > 0)
  
  # New state-level df for analyses of the non-blip novel communities through time
  nov.summary <- rbindlist(lapply(nov.comm.no.blip$site, function(ID){
    
    # Identify the novel points using the length we calculated
    print(ID)
    id <- ID # start with a single to see how it works
    
    # Extract length of the novel state
    state.length <- as.double(unique(nov.comm.no.blip[which(nov.comm.no.blip$site == id), "novel.length"]))
    
    # Use that length to extract all the points that make up the novel state
    
    year.start <- as.double(min(nov.comm.no.blip[which(nov.comm.no.blip$site == id), "bins"]))
    
    year.end <- year.start - state.length # Easy way to find the endpoints
    
    # Add a little variable that tracks whether or not the state is "still counting"
    # We compare the nov. comm. end bin to the timeseries end bin.
    
    pres.stat <- ifelse(year.end == min(full.novel.mat[full.novel.mat$site == id, "bins"]), "Still counting", "Over")
    
    
    # Now we extract the points from our df.
    
    years.nov <- c(year.start:year.end) # this will include values with no communities, but that's fine
    
    nov.stat.comm <- subset(full.novel.mat, site == id & bins %in% years.nov) # There are the communities within a novel state!
    
    # There is always a transition period. Sometimes this is a bit longer than one year,
    # I think it would be fair to add this to the novel community, just to count ext and orig.
    temp <- subset(full.novel.mat, site == id)
    trans.gap <- temp$bins[which(temp$cat == "novel")-1] - temp$bins[which(temp$cat == "novel")] - 1
    transition.period <- ifelse(is_empty(trans.gap), 0, trans.gap)
    
    # compute average diversity metrics, we want PER YEAR, NOT PER TIMEJUMP for extinctions and originations
    
    mean.evenness.n <- mean(nov.stat.comm$evenness) # normal mean is ok here
    mean.shannon.n <- mean(nov.stat.comm$shannon.d)
    mean.ext.n <- sum(nov.stat.comm$ext, na.rm = T)/(state.length +transition.period) # more accurate this way
    mean.orig.n <- sum(nov.stat.comm$orig, na.rm = T)/(state.length +transition.period)
    
    # I think we also need to extract the non-novel communities to do comparisons
    
    # Save total years non-novel in case we want it as a covariate
    
    years.non.nov <- max(temp$bins) - min(temp$bins) - state.length - transition.period
    
    # Non-novel comms
    back.stat.comm <- subset(full.novel.mat, site == id & bins %!in% years.nov)
    
    mean.evenness.b <- mean(back.stat.comm$evenness)
    mean.shannon.b <- mean(back.stat.comm$shannon.d)
    mean.ext.b <- sum(back.stat.comm$ext, na.rm = T)/years.non.nov
    mean.orig.b <- sum(back.stat.comm$orig, na.rm = T)/years.non.nov
    
    # Invaders
    mean.invaders.n <- mean(nov.stat.comm$INC)
    mean.invaders.b <- mean(back.stat.comm$INC)
    
    # Return a new df with state-level metrics
    
    state.level.metrics <- data.frame(matrix(data = NA, nrow =2))
    
    state.level.metrics$state <- NA
    state.level.metrics$length <- NA
    state.level.metrics$mean.even <- NA
    state.level.metrics$mean.shannon <- NA
    state.level.metrics$mean.orig <- NA
    state.level.metrics$mean.ext <- NA
    state.level.metrics$site <- id
    state.level.metrics[,1] <- NULL # remove the auto column 
    state.level.metrics$country <- unique(temp$country)
    state.level.metrics$invaders <- NA
    state.level.metrics$begin <- year.start
    state.level.metrics$end <- year.end
    state.level.metrics$status <- pres.stat
    
    # There's one more fun thing we can add, the ratio between non-novel and novel for orig and ext.
    # It will just be entered for the novel community in the df
    
    # We can now fill in the new df based on some distinctions between states,
    # this looks a bit tedious but is actually the quickest way.
    
    state.level.metrics[1, c("state", "length", "mean.even", 
                             "mean.shannon", "mean.orig", "mean.ext", "invaders")] <- c("novel", state.length, mean.evenness.n,
                                                                                        mean.shannon.n, mean.orig.n, mean.ext.n, mean.invaders.n)
    state.level.metrics[2, c("state", "length", "mean.even", 
                             "mean.shannon", "mean.orig", "mean.ext", "invaders")] <- c("back", years.non.nov, mean.evenness.b,
                                                                                        mean.shannon.b, mean.orig.b, mean.ext.b, mean.invaders.b)
    
    
    
    return(state.level.metrics)
    
  }))
  
  # Add a status variable for persistence type AND for if the state ended before the time series
  nov.comm.all$type <- NA
  nov.comm.all$status <- NA
  
  # Loops to fill in status and type data
  for (i in 1:nrow(nov.comm.all)){
    if (nov.comm.all$site[i] %in% nov.comm.no.blip$site){
      nov.comm.all$type[i] <- "persistent"
    } else {nov.comm.all$type[i] <- "blip"}
  }
  for (i in 1:nrow(nov.comm.all)){
    print(i)
    for (j in 1:nrow(nov.summary)){
      if(nov.comm.all$site[i] == nov.summary$site[j]){
        nov.comm.all$status[i] <- nov.summary$status[j]
      }
    }
  }
  
  # We now have a complete data frame with all the information needed to do some evaluations
  
  # How many are blips?
  perc.blip <- nrow(subset(nov.comm.all, type == "blip"))/(nrow(subset(nov.comm.all, type == "persistent")) + nrow(subset(nov.comm.all, type == "blip")))
  perc.no.blip <- 1-perc.blip
  
  # How many persist till the end of the time series?
  no.end.comm <- nrow(subset(nov.comm.all, status == "Still counting"))/(nrow(subset(nov.comm.all, status == "Still counting")) + nrow(subset(nov.comm.all, status == "Over")))
  
  # Collate output
  output.df <- data.frame("Type" = c("Total","Novel", "Blips", "No Blip", "Unended", "Ended"), 
                          "Absolute"= c(nrow(full.novel.mat), nrow(nov.comm.all), (nrow(nov.comm.all) - nrow(nov.comm.no.blip)), nrow(nov.comm.no.blip), nrow(subset(nov.comm.all, status == "Still counting")), nrow(subset(nov.comm.all, status == "Over"))))
  
  return(list("summary"=output.df, state.level.metrics = nov.summary))
}
