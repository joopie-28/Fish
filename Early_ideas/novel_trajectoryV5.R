# Timothy's function for legacy work.
# I've adjusted it to our data.

novel.trajectoryV5 <- function(novel.object,
                               dissim.method,
                               novel.crop = c(-Inf, Inf)){
  
  print("Triangulating observed trajectory")
  
  sites <- names(novel.object)
  
  
  novel.traj.list <- lapply(c("instant", "cumul", "novel"), function(novel.cat){
    print(novel.cat) 
    
    do.call("rbind", lapply(sites, function(s){
      #print(s)
      temp.nov <- novel.object[[s]]$novelty
      temp.ssmat <- novel.object[[s]]$matrix
      print(s)
      # find position of novel communities
      nov.locs <- as.character(temp.nov$bins)[temp.nov[,novel.cat]]
      
      # subset novel communities to just those that occur in certain time band
      nov.locs <- nov.locs[as.numeric(nov.locs) >= novel.crop[1] &
                             as.numeric(nov.locs) <= novel.crop[2]]
      nov.locs <- nov.locs[!is.na(nov.locs)]
      
      if(length(nov.locs) == 0){return(NULL)}
      
      # return data-frame of post-novel trajectories
      do.call("rbind", lapply(nov.locs, function(l){
        
        # rotate ssmat if we need to
        if(as.numeric(as.character(rownames(temp.ssmat)[1])) < 
           as.numeric(as.character(rownames(temp.ssmat)[dim(temp.ssmat)[1]]))){
          
          temp.ssmat <- temp.ssmat[dim(temp.ssmat)[1]:1, ]
          
        }
        
        # only keep going if novel comm is at least 2 time points from time series end
        if(which(rownames(temp.ssmat) == l)+1 >= nrow(temp.ssmat)){return(NULL)}
        
        # subset ssmat from novel time point -1 ("pre-novel")
        sub.ssmat <- temp.ssmat[(which(rownames(temp.ssmat) == l)-1):nrow(temp.ssmat),]
        sub.dist <- as.matrix(vegdist(sub.ssmat, method=dissim.method))
        
        # create data-frame of dissimilarities after novel comm
        dist.df <- data.frame(sub.dist[-(1:2),1:2])
        
        if(ncol(dist.df)==1){return(NULL)}
        
        colnames(dist.df) <- c("dP", "dN")
        dist.df$delta.dP <- diff(sub.dist[-1,1])
        dist.df$delta.dN <- c(NA, diff(sub.dist[-(1:2),2]))
        dist.df$dP1 <- dist.df$dP[1]
        dist.df$dN1 <- dist.df$dN[1]
        dist.df$bin <- as.numeric(rownames(dist.df))
        dist.df$novel.time <- as.numeric(l)
        dist.df$time.since.novel <- abs(dist.df$bin - as.numeric(l))
        dist.df$site <- s
        dist.df$novelID <- paste0(dist.df$site, ".", dist.df$novel.time)
        
        return(dist.df)
      }))
      
    }))
    
  })
  names(novel.traj.list) <- c("instant", "cumul", "novel")
  
  # now create the same object using randomly selected, non-novel comms
  print("Triangulated non-novel trajectory")
  
  rand.traj.df <- do.call("rbind", lapply(sites, function(s){
    #print(s)
    temp.nov <- novel.object[[s]]$novelty
    temp.ssmat <- novel.object[[s]]$matrix
    
    # find position of novel communities
    back.locs <- as.character(temp.nov$bins)[temp.nov$cat == "back"]
    
    # subset novel communities to just those that occur in certain time band
    back.locs <- back.locs[as.numeric(back.locs) >= novel.crop[1] &
                             as.numeric(back.locs) <= novel.crop[2]]
    back.locs <- back.locs[!is.na(back.locs)]
    
    if(length(back.locs) == 0){return(NULL)}
    
    # return data-frame of post-novel trajectories
    do.call("rbind", lapply(back.locs, function(l){
      
      # rotate ssmat if we need to
      if(as.numeric(as.character(rownames(temp.ssmat)[1])) < 
         as.numeric(as.character(rownames(temp.ssmat)[dim(temp.ssmat)[1]]))){
        
        temp.ssmat <- temp.ssmat[dim(temp.ssmat)[1]:1, ]
        
      }
      
      # only keep going if novel comm is at least 2 time points from time series end
      if(which(rownames(temp.ssmat) == l)+1 >= nrow(temp.ssmat)){return(NULL)}
      
      # subset ssmat from novel time point -1 ("pre-novel")
      sub.ssmat <- temp.ssmat[(which(rownames(temp.ssmat) == l)-1):nrow(temp.ssmat),]
      sub.dist <- as.matrix(vegdist(sub.ssmat, method=dissim.method))
      
      # create data-frame of dissimilarities after novel comm
      dist.df <- data.frame(sub.dist[-(1:2),1:2])
      
      if(ncol(dist.df)==1){return(NULL)}
      
      colnames(dist.df) <- c("dP", "dN")
      dist.df$delta.dP <- diff(sub.dist[-1,1])
      dist.df$delta.dN <- c(NA, diff(sub.dist[-(1:2),2]))
      dist.df$dP1 <- dist.df$dP[1]
      dist.df$dN1 <- dist.df$dN[1]
      dist.df$bin <- as.numeric(rownames(dist.df))
      dist.df$novel.time <- as.numeric(l)
      dist.df$time.since.novel <- abs(dist.df$bin - as.numeric(l))
      dist.df$site <- s
      dist.df$novelID <- paste0(dist.df$site, ".", dist.df$novel.time)
      
      return(dist.df)
    }))
    
  }))
  
  return(list(novel.traj.list = novel.traj.list,
              rand.traj.list = rand.traj.df))  
  
}

# The inputs are a 'fish_communities' list, whether you want to filter your data or not
# and the R statistic you want to filter at (between 0 and 1). You can also choose to
# plot the data in a density map.

nov.len.traj.function <- function(fish.list, filter = T, R, plot){
  
# Extract the ssmats
  raw.ssmats <- do.call(c, fish.list$BioRealm_Matrices_A)

  indices.s <- lapply(names(raw.ssmats), function(index){
  
    out <- str_split(index, pattern = "A.")[[1]][2]
  
    return(out)
  })

  names(raw.ssmats) <- indices.s

# Extract the novelty parameter frame
  raw.nov <- do.call(c, fish.list$BioRealm_Novelty_A)

  indices.n <- lapply(names(raw.nov), function(index){
  
    out.n <- str_split(index, pattern = "A.")[[1]][2]
  
    return(out.n)
  })

  names(raw.nov) <- indices.n
  
  
# Create new object in format for novel trajectories  
  for (i in 1:length(raw.nov)){
    print(i)
  
    for (j in 1:length(raw.ssmats)){
      if (names(raw.nov)[[i]] == names(raw.ssmats)[[j]]){
    
        raw.nov[[i]] <- list("novelty" = raw.nov[[i]], "matrix" = raw.ssmats[[j]])
    
      }
    }
  
  }

# This is an optional filtering section, which allows us to look at novel communities 
# with different R's
  
# uses the novelty length filter function and true novelty output object
# True novelty output is a list of novel communities and their respective R values.
# Excluding those novel communities that occur at the end of their timeseries.
  
  
  if(filter){
    
    print(paste0("Filtering your data based on your R statistic of ", R))
  
    filtered.novelty <- novel.length.checker(true.novelty.output, cut.off = R)
    site.list <- list()
    for (i in 1:length(filtered.novelty)){
      # Create new site list based on filtered values
      site.list <- c(site.list, filtered.novelty[[i]]$ID)
    }
    
    # Now remove values based on presence in filtered list and replace with NA's
    
    for (i in 1:length(raw.nov)){
      
      if(names(raw.nov)[[i]] %notin% site.list){
        
        raw.nov[[i]] <- NA
      }
    }
    
    # Lastly, filter out NA's
    
    raw.nov <- raw.nov[!sapply(raw.nov, anyNA)]
    
  }
  
  output <- list("trajectories" = novel.trajectoryV5(raw.nov, dissim.method = "bray"), "lengths" = filtered.novelty)
  
  if(plot){
    
    p1 <- ggplot(data = output$trajectories$novel.traj.list$novel, aes(x=dN, y=dP)) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
      scale_fill_distiller(palette= "Spectral", direction=-1) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) + ggtitle(paste0("R = ", R)) +
      labs(title = paste0("R = ", R), x = "Distance from Novel State", y = "Distance from Last Pre-Novel State") +
      theme(
        legend.position='topright'
      )
  }
  print(p1)
  return(output)
}

len.traj.full <- nov.len.traj.function(Fish_Communities_A, R = 0.99, plot = T)

 
# Make a GIF!

gif.list <- lapply(c(1:20), FUN = function(x){
  
  test <- subset(len.traj.full$trajectories$novel.traj.list$novel, time.since.novel == x)
  
  
  img.list <-  ggplot(data = test, aes(x=dN, y=dP)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F) +
    scale_fill_distiller(palette="Spectral", direction=-1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position='topright'
    ) + ggtitle(paste0("Time after novelty = ", x))
  
  return(img.list)
})
img.list <- list()

for (i in 1:20){
  print(i)
  ggsave(plot = gif.list[[i]], filename = paste0(i,".png"))
  img.list <- c(img.list, paste0(i, ".png"))
}

img.read <- lapply(img.list, image_read)
img.joined <- image_join(img.read)
animated.img <- image_animate(img.joined, fps = 0.5)
image_write(image = animated.img,
            path = "novelty.gif")

# If dP is < dN, this implies a return to pre-novel conditions








