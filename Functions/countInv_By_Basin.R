#### Invader presence per basin ####

countInv_By_Basin <- function(time_series_data, basin){
  
  temp <- subset(time_series_data, HydroBasin == basin)
  country <- unique(temp$Country)
  # if more than one country, just take majority
  if(length(country > 1)){
   country <- names(which.max(table(temp$Country)))
  }
  sites <- temp$TimeSeriesID
  basinSub <- subset(Survey_Data, TimeSeriesID %in% sites)
  
  basin_df <- as.data.frame(tapply(basinSub$Abundance,
                                   list(basinSub$Year, basinSub$Species),
                                   sum, na.rm = T))
  
  # Convert absences (NA) to 0's
  basin_df[is.na(basin_df)]<-0
  
  # Convert to relative abundance
  basin.mat <- basin_df/rowSums(basin_df)
  
  # Convert rownames to time in past
  rownames(basin.mat) <- 2021-as.numeric(rownames(basin.mat))
  
  # extract country information to match with invader database using mget()
  
  invader_list <- get0(tolower(paste0(country, '_invaders')))
  # Might want to extend this to include NNC later..

  invader_by_basin <- ifelse(basin.mat[,] > 0, 1, 0) |>
    as.data.frame() |>
    # only do this if we have countries with invaders
    select(colnames(basin.mat)[colnames(basin.mat) %in% invader_list$Species]) |>
    rowSums() |>
    as.data.frame() |>
    mutate(Basin = basin,
           Country = country)
  
  colnames(invader_by_basin) <- c('total_inv', 'basin', 'country')
  invader_by_basin$bins <- rownames(invader_by_basin)
  invader_by_basin$year <- 2021-as.numeric(rownames(invader_by_basin))
  return(invader_by_basin)
}

# we want to model chance of emergence, chance of being persistent, and length of persistence. These are our core topics.

# 1. has been done, and can be done using full.novel.mat and any covariates in that. invbybasin adds the basin level invaders
# 2. can be done by a frame containing nov communities at a time series level. These include environmental variables at a time series level
# but also community-specific characteristics such as evenness, richness..

# explore how these analyses can be redone at different HYBAS levels.

countEx_By_Basin <- function(time_series_data, basin){
  
  temp <- subset(time_series_data, HydroBasin == basin)
  country <- unique(temp$Country)
  # if more than one country, just take majority
  if(length(country > 1)){
    country <- names(which.max(table(temp$Country)))
  }
  sites <- temp$TimeSeriesID
  basinSub <- subset(Survey_Data, TimeSeriesID %in% sites)
  
  basin_df <- as.data.frame(tapply(basinSub$Abundance,
                                   list(basinSub$Year, basinSub$Species),
                                   sum, na.rm = T))
  
  # Convert absences (NA) to 0's
  basin_df[is.na(basin_df)]<-0
  
  # Convert to relative abundance
  basin.mat <- basin_df/rowSums(basin_df)
  
  # Convert rownames to time in past
  rownames(basin.mat) <- 2021-as.numeric(rownames(basin.mat))
  
  # extract country information to match with invader database using mget()
  invader_list <-  get0(tolower(paste0(country, '_country_level'))) |>
    dplyr::filter(Status == 'Native') # native, because we will use %!in% downstream
    
  
  # Ask: any non-native species in the basin?
  
  invader_by_basin <- ifelse(basin.mat[,] > 0, 1, 0) |>
    as.data.frame() |>
    # only do this if we have countries with invaders
    dplyr::select(colnames(basin.mat)[colnames(basin.mat) %!in% invader_list$Species]) |>
    rowSums() |>
    as.data.frame() |>
    mutate(Basin = basin,
           Country = country)
  
  colnames(invader_by_basin) <- c('total_inv', 'basin', 'country')
  invader_by_basin$bins <- rownames(invader_by_basin)
  invader_by_basin$year <- 2021-as.numeric(rownames(invader_by_basin))
  return(invader_by_basin)
}
