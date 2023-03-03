# Filters time series by the lag of their nvoel community to the next community of any type.

filter.by.lag <- function(full.novel.mat.season, n_lag){
  full.novel.mat.season$lag_to_next <- NA
  for (site.ID in unique(full.novel.mat.season$site)){
    print(site.ID)
    temp <- subset(full.novel.mat.season, site == site.ID)
    indices <- which(full.novel.mat.season$site == site.ID)
    for (i in 1:nrow(temp)){
      temp$lag_to_next[i] <- abs(temp$bins[i+1] - temp$bins[i])
      
    }
    if(any(temp$cat == 'novel')){
      nov.index <- which(temp$cat == 'novel')
      for (j in nov.index){
        print(j)
        if(temp$lag_to_next[j] > n_lag & temp$novel.class[j] != 'END'){
          temp$lag_to_next <- rep(1000, nrow(temp))
        }
      }
    }
    full.novel.mat.season$lag_to_next[indices] <- temp$lag_to_next
  }
  return(full.novel.mat.season)
}

# Remove time series where the novel community emerges due to lag
filter.by.lag.novel <- function(full.novel.mat.season, n_lag){
  full.novel.mat.season$lag_to_next <- NA
  for (site.ID in unique(full.novel.mat.season$site)){
    print(site.ID)
    temp <- subset(full.novel.mat.season, site == site.ID)
    indices <- which(full.novel.mat.season$site == site.ID)
    for (i in 1:nrow(temp)){
      temp$lag_to_next[i] <- abs(as.numeric(temp$bins[i+1]) - as.numeric(temp$bins[i]))
      
    }
    
    if(any(temp$cat == 'novel')){
      nov.index <- which(temp$cat == 'novel')
      
        prev.lag<-temp[nov.index, 'bin_lag']
        if(any(prev.lag > n_lag)){
          full.novel.mat.season$lag_to_next[indices] <- rep(1000, nrow(temp))
      }
    } 
  }
  return(full.novel.mat.season)
}

