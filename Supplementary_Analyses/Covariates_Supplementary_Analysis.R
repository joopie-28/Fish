##################################################################
# Novel Communities of Freshwater Fishes: a global investigation #
# into their emergence and persistence over the last 50 years	####
##################################################################

################ Supplementary analyses ######################

novel_LagPosition_full <- function(input_list){
  
  # Keep track of iteration
  i <- 1
  
  nov.mat <- lapply(input_list , function(mat){
    
    # Simple message for the viewer
    print(paste0(i, " of ", length(input_list)))
    
    # provide a site name
    ID <- names(input_list)[i]
    
    # Run novelty detection framework
    temp <- identify.novel.gam(site.sp.mat = mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
                               gam.max.k = -1)
    
    temp$position <- 1:nrow(temp)
    
    i <<- i + 1
    
    
    return(temp)
  })
  
  return(nov.mat)
}
novel_LagPosition_filtered <- function(input_list){
  
  # Keep track of iteration
  i <- 1
  
  nov.mat <- lapply(input_list , function(mat){
    
    # Simple message for the viewer
    print(paste0(i, " of ", length(input_list)))
    
    # provide a site name
    ID <- names(input_list)[i]
    
    # Run novelty detection framework
    temp <- identify.novel.gam(site.sp.mat = mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
                               gam.max.k = -1)
    # Remove first 5 bins
    temp <- temp[-c(1:5),]
    # Rename position 1:n, after filtering
    temp$position <- 1:nrow(temp)
    
    i <<- i + 1
    
    
    return(temp)
  })
  
  return(nov.mat)
}
LagConsistencyFilter_Supp <- function(fullPos_Lag_df,nlag){
  
  for(site in fullPos_Lag_df$site){
    print(site)
    # Subset the data by site
    sub <- fullPos_Lag_df[which(fullPos_Lag_df$site == site),]
    nov.lag <- sub$bin.lag[which(sub$novel == 1)]
    
    
    if(all(sub$novel == 0, na.rm = T)){
      next
    }
    # If true, keep the data. We also remove those time series with lag of +nlag!
    if(any(nov.lag < nlag)){
      next
    }
    else{
      # If false, the novel community is unreliable and we remove it and the time series
      print('Removing these rows')
      fullPos_Lag_df <- fullPos_Lag_df[-which(fullPos_Lag_df$site == site),]
    }
  }
  return(fullPos_Lag_df)
}

# Full dataset
fullPos_Lag <- novel_LagPosition_full(matrix_list_seasonality)

# First 5 removed
First5_Pos_Lag <- novel_LagPosition_filtered(matrix_list_seasonality)

# Full dataset to df
fullPos_Lag_df <- rbindlist(fullPos_Lag) |>
  dplyr::select(site, novel, cumul, instant, bin.lag, position)|>
  separate(site, c('Site', 'Quarter'), remove = F) 

# Truncated dataset, where first 5 bins removed and lag > 3 removed
First5_Pos_NoLag_df  <- rbindlist(First5_Pos_Lag) |>
  dplyr::select(site, novel, cumul, instant, bin.lag, position) |>
  separate(site, c('Site', 'Quarter'), remove = F) |>
  LagConsistencyFilter_Supp(nlag = 4)

# Inspect relationships in full data
mod.s1 <-glmer(novel ~ bin.lag+position+ (1|Quarter/Site), 
      data = fullPos_Lag_df, family= 'binomial') 
summary(mod.s1)

# Observe renewed relationship after removal of first 5 bins
mod.s2 <-glmer(novel ~ bin.lag+position+(1|Quarter/Site), 
               data = First5_Pos_NoLag_df, family= 'binomial') 
summary(mod.s2) ## lag non-significant, position remains problematic










summary(glmer(novel ~ NNC*Latitude+ppd_pk_uav+ (1|Quarter/site_ID), family = 'binomial', data=FullEnvFrame))
