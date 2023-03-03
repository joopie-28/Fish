# Functions for creating species-by-community matrices incorporating seasonality


create.seasonal.ssmat <- function(n.ID){
  
  # Isolate time series data from overall survey data
  df<-subset(Survey_Data, TimeSeriesID == n.ID )
  
  # Create data frames for each timeseries-quarter pair
  output<-lapply(unique(df$Quarter), function(quarter){
    
    sub.1 <- subset(df, Quarter == quarter)
    
    sub.2<-as.data.frame(tapply(sub.1$Abundance,
                                list(sub.1$Year, sub.1$Species),
                                sum, na.rm = T))
    
    # Convert absences (NA) to 0's
    sub.2[is.na(sub.2)]<-0
    
    # Convert to relative abundance
    ss.mat <- sub.2/rowSums(sub.2)
    
    # Convert rownames to time in past
    rownames(ss.mat) <- 2021-as.numeric(rownames(ss.mat))
    
    # Filter time series with less than 10 timesteps and at least 2 species
    if(nrow(ss.mat) > 9 & ncol(ss.mat) > 1){
      return(ss.mat)
    }
    else{
      return(NULL)
    }
    
  })
  # Assign merger names for later use
  names(output) <- paste0(n.ID,".", unique(df$Quarter))
  return(output)
}

list_matrix_seasonality_function <- function(ID_list){
  ids <- unlist(ID_list)
  
  # Create community-by-species matrices for each season-timeseries pair
  mats<-lapply(ids, function(x){
    print(x)
    create.seasonal.ssmat(x)
  })
  
  # Homgenise list
  output<-unlist(mats, recursive = F)
  output.final<-Filter(Negate(is.null), output)
  
  return(output.final)
}

novelty.detection.gam <- function(input_list){
  
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
    i <<- i + 1
    
    # Remove first 5 bins
    temp <- temp[-c(1:5),]
    
    
    return(temp)
  })
  
  return(nov.mat)
}

create.mod.frame <- function(novel.output){
  big.df<-rbindlist(nov.output)
  
  temp_df <- data.frame(do.call("rbind", strsplit(as.character(big.df$site), ".",
                                                  fixed = TRUE)))
  big.df[,c("site")] <-NULL
  names(temp_df) <- c("site", "Quarter")
  
  return.frame<-cbind(temp_df, big.df)
  return(return.frame)
}
