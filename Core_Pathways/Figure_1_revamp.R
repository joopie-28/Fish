# Revamping figure 1

# Extract basins

plot.basin.summary <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  
  # Set up plotting background
  palearctic.plot <- ne_countries(scale = "Large", continent = 'Europe', type = 'sovereignty',
                                  returnclass = c("sp", "sf"))

  ymin = 40
  ymax = 70
  xmin= -10
  xmax= 35
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(palearctic.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  
  
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-geo.timeseries.full |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(Country != "CIV")|>
    mutate(binary_novel = ifelse(novel > 0, 1, 0)) |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('Site' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'binary_novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$binary_novel == 0),'observed'], non.nov)
      nov <- c(basin.sub[which(basin.sub$binary_novel == 1),'observed'], nov)
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$binary_novel == 0),'observed']+
                             basin.sub[which(basin.sub$binary_novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$binary_novel == 1),'observed']/basin.sub[which(basin.sub$binary_novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = rbp(10)[as.numeric(cut(plotting.frame$proportion_novel,breaks =10))],
       border=rgb(0,0,0,0),
       add =T)
  
}

pdf(file = "/Users/sassen/Desktop/Figure_test.pdf",
    width = 15,
    height = 12.5)
plot.basin.summary.bycommunity.PAL(full.novel.mat.season, HYBAS_scheme, 8)
plot.basin.summary.bycommunity.NEA(full.novel.mat.season, HYBAS_scheme, 8)
plot.basin.summary.bycommunity.AUS(full.novel.mat.season, HYBAS_scheme, 8)
dev.off()


# Need a better colour ramp
plot.basin.summary.bycommunity.PAL <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  
  # Set up plotting background
  palearctic.plot <- ne_countries(scale = "Large", continent = 'Europe', type = 'sovereignty',
                                  returnclass = c("sp", "sf"))
  
  ymin <<- 40
  ymax <<- 70
  xmin= -10
  xmax= 35
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(palearctic.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = T, bg = 'white', xaxs = 'i', yaxs = 'i')
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))

  
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-full.novel.mat.season |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(country != "CIV") |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('site_ID' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$novel == 0),'observed'], non.nov)
      nov <- c(ifelse(nrow(basin.sub[which(basin.sub$novel == 1),'observed'])==0,0,
                      basin.sub[which(basin.sub$novel == 1),'observed']), nov)
      
      # Continue here... wip
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$novel == 0),'observed']+
                             basin.sub[which(basin.sub$novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$novel == 1),'observed']/basin.sub[which(basin.sub$novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Add colour ramp
  colour.vec<- rbp(5)
  prop.nov.ByBasin <- prop.nov.ByBasin |>
    mutate(colramp = ifelse(proportion_novel == 0, colour.vec[1],
                            ifelse(proportion_novel > 0 & proportion_novel <= 0.05,colour.vec[2],
                                   ifelse(proportion_novel <= 0.10 & proportion_novel > 0.05 ,colour.vec[3],
                                          ifelse(proportion_novel <= 0.15 & proportion_novel > 0.10,colour.vec[4],colour.vec[5])))))
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = plotting.frame$colramp,
       border=rgb(0,0,0,0),
       add =T)

}

plot.basin.summary.bycommunity.NEA <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  
  # Set up plotting background
  nearctic.plot <- ne_countries(scale = "Large", geounit = c('United States of America', 'Canada', "Mexico"), type = 'sovereignty',
                                returnclass = c("sp", "sf"))
  ymin = 25
  ymax = 50
  xmin= -130
  xmax= -55
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(nearctic.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  
  
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_na_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_na_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-full.novel.mat.season |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(country != "CIV") |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('site_ID' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$novel == 0),'observed'], non.nov)
      nov <- c(ifelse(nrow(basin.sub[which(basin.sub$novel == 1),'observed'])==0,0,
                      basin.sub[which(basin.sub$novel == 1),'observed']), nov)
      
      # Continue here... wip
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$novel == 0),'observed']+
                             basin.sub[which(basin.sub$novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$novel == 1),'observed']/basin.sub[which(basin.sub$novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Add colour ramp
  colour.vec<- rbp(5)
  prop.nov.ByBasin <- prop.nov.ByBasin |>
    mutate(colramp = ifelse(proportion_novel == 0, colour.vec[1],
                            ifelse(proportion_novel > 0 & proportion_novel <= 0.05,colour.vec[2],
                                   ifelse(proportion_novel <= 0.10 & proportion_novel > 0.05 ,colour.vec[3],
                                          ifelse(proportion_novel <= 0.15 & proportion_novel > 0.10,colour.vec[4],colour.vec[5])))))
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = plotting.frame$colramp,
       border=rgb(0,0,0,0),
       add =T)
  
}

plot.basin.summary.bycommunity.AUS <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  
  # Set up plotting background
  aus.plot <- ne_countries(scale = "Large", geounit = c('australia', 'new zealand', 'papua new guinea', 'indonesia'), type = 'map_units',
                           returnclass = c("sp", "sf"))
  ymin = -47 
  ymax = -7
  xmin = 110 
  xmax = 160
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(aus.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  
  segments(x0=149, y0=-29, x1 =149 , y1 =-25 )
  segments(x0=149, y0=-29, x1 =155 , y1 =-29 )
  segments(x0=155, y0=-29, x1 =155 , y1 =-25 )
  segments(x0=155, y0=-25, x1 =149 , y1 =-25 )
  segments(x0 = 155, y0 = -29, x1=135.2, y1= -47.6   )
  segments(x0 = 149, y0 = -25, x1=,116, y1= -29   )
  par(fig = c(0.175, 0.54, 0.145, 0.45), new = T, mai= c(0,0,0,0), mar=c(0,0,0,0)) 
  qld_bbox = st_bbox(c(ymin = -29, ymax = -25, xmin =149, xmax = 160))
  plot(crop(ne_countries(scale = 'Large', country = c('Australia'), type = 'sovereignty',
                         returnclass = c("sp", "sf")), qld_bbox), col = 'lightgrey', bg= 'white')
  box(which = 'plot',  lty = 1)
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_au_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_au_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-full.novel.mat.season |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(country != "CIV") |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('site_ID' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$novel == 0),'observed'], non.nov)
      nov <- c(ifelse(nrow(basin.sub[which(basin.sub$novel == 1),'observed'])==0,0,
                      basin.sub[which(basin.sub$novel == 1),'observed']), nov)
      
      # Continue here... wip
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$novel == 0),'observed']+
                             basin.sub[which(basin.sub$novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$novel == 1),'observed']/basin.sub[which(basin.sub$novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Add colour ramp
  colour.vec<- rbp(5)
  prop.nov.ByBasin <- prop.nov.ByBasin |>
    mutate(colramp = ifelse(proportion_novel == 0, colour.vec[1],
                            ifelse(proportion_novel > 0 & proportion_novel <= 0.05,colour.vec[2],
                                   ifelse(proportion_novel <= 0.10 & proportion_novel > 0.05 ,colour.vec[3],
                                          ifelse(proportion_novel <= 0.15 & proportion_novel > 0.10,colour.vec[4],colour.vec[5])))))
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = plotting.frame$colramp,
       border=rgb(0,0,0,0),
       add =T)
  
}




plot.basin.summary.bycommunity.PAL(full.novel.mat.season, HYBAS_scheme, 8)

# NEED TO:
# make sure seasonality is not influencing these numbers
# Colour legend

complete.basin.novelty.plotter.PAL <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  # Define layout of the plot
  layout(matrix(c(1,2), ncol=2), widths=c(2,1), height=4)
  
  # Create the plot of Europe with hydrobasins
  par(mar = c(5, 3, 0, 0))
  
  # Set up plotting background
  palearctic.plot <- ne_countries(scale = "Large", continent = 'Europe', type = 'sovereignty',
                                  returnclass = c("sp", "sf"))
  
  ymin <- 40
  ymax <- 70
  xmin= -10
  xmax= 35
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(palearctic.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = F, bg = 'white', xaxs = 'i', yaxs = 'i')
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_eu_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-full.novel.mat.season |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(country != "CIV") |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('site_ID' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$novel == 0),'observed'], non.nov)
      nov <- c(ifelse(nrow(basin.sub[which(basin.sub$novel == 1),'observed'])==0,0,
                      basin.sub[which(basin.sub$novel == 1),'observed']), nov)
      
      # Continue here... wip
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$novel == 0),'observed']+
                             basin.sub[which(basin.sub$novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$novel == 1),'observed']/basin.sub[which(basin.sub$novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Add colour ramp
  colour.vec<- rbp(5)
  prop.nov.ByBasin <- prop.nov.ByBasin |>
    mutate(colramp = ifelse(proportion_novel == 0, colour.vec[1],
                            ifelse(proportion_novel > 0 & proportion_novel <= 0.05,colour.vec[2],
                                   ifelse(proportion_novel <= 0.10 & proportion_novel > 0.05 ,colour.vec[3],
                                          ifelse(proportion_novel <= 0.15 & proportion_novel > 0.10,colour.vec[4],colour.vec[5])))))
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = plotting.frame$colramp,
       border=rgb(0,0,0,0),
       add =T)
  
  # Next bit
  
  ymin.map = par('usr')[3] 
  ymax.map = par('usr')[4]
  # define the bin width
  bin_width <- 2
  
  df <- full.novel.mat.season |>
    dplyr::filter(BioRealm == 'Palearctic')
  df_nov <- subset(df, novel == 1)
  # create a new column 'lat_group' by dividing latitude by the bin width and rounding down
  df$lat_group <- floor(df$Latitude / bin_width) * bin_width
  df_nov$lat_group <- floor(df_nov$Latitude / bin_width) * bin_width
  
  # group the dataframe by 'lat_group' and count the number of rows
  df_counts <- aggregate(df$Latitude, by = list(df$lat_group), length)
  
  df_counts_nov <- aggregate(df_nov$Latitude, by = list(df_nov$lat_group), length)
  
  # create a new dataframe with latitude range and number of observations
  range_df <- data.frame(lat_range = paste(ifelse(df_counts$Group.1 >= 0, paste0(df_counts$Group.1, 'N'), paste0(abs(df_counts$Group.1), 'S')),
                                           ifelse(df_counts$Group.1 + bin_width - 1 >= 0, paste0(df_counts$Group.1 + bin_width - 1, 'N'), paste0(abs(df_counts$Group.1 + bin_width - 1), 'S')),
                                           sep = '-'),
                         no_obs = df_counts$x,
                         min_lat = df_counts$Group.1) |>
    dplyr::filter(min_lat >= ymin)
  
  range_nov_df <- data.frame(lat_range = paste(ifelse(df_counts_nov$Group.1 >= 0, paste0(df_counts_nov$Group.1, 'N'), paste0(abs(df_counts_nov$Group.1), 'S')),
                                               ifelse(df_counts_nov$Group.1 + bin_width - 1 >= 0, paste0(df_counts_nov$Group.1 + bin_width - 1, 'N'), paste0(abs(df_counts_nov$Group.1 + bin_width - 1), 'S')),
                                               sep = '-'),
                             no_obs = df_counts_nov$x,
                             min_lat = df_counts_nov$Group.1) |>
    dplyr::filter(min_lat >= ymin)
  
  
  # Create the plot of Europe with hydrobasins
  par(mar = c(5, 0, 0, 1))
  
  range_df$propnovel <- range_nov_df$no_obs/range_df$no_obs
  
  plot(NULL, xlim=c(0,max(range_df$propnovel)), 
       ylim=c(ymin.map,ymax.map-1.3), axes= F, xlab = "")
  axis(side=1)
  mtext("Proportion of Novel Communities", 
        side=1,line=3, cex.lab=1)
  
  for(i in 1:(nrow(range_df)-1)){
    print(i)
    polygon(y=c(range_df$min_lat[i], range_df$min_lat[i+1],range_df$min_lat[i+1] ,range_df$min_lat[i] ),
            x= c(0,0,range_df$propnovel[i], range_df$propnovel[i]), col='grey70')
    # polygon(y=c(range_df$min_lat[i], range_df$min_lat[i+1],range_df$min_lat[i+1] ,range_df$min_lat[i] ),
    #      x= c(0,0,range_df$no_obs[i], range_df$no_obs[i]), col = "grey70")
    #polygon(y=c(range_nov_df$min_lat[i], range_nov_df$min_lat[i+1],range_nov_df$min_lat[i+1] ,range_nov_df$min_lat[i] ),
    #       x= c(0,0,range_nov_df$no_obs[i], range_nov_df$no_obs[i]), col = "orange")
    
    
  }
  
}

complete.basin.novelty.plotter.NEA <- function(full.novel.mat.season, HYBAS_scheme, HYBAS_level){
  # Define layout of the plot
  layout(matrix(c(1,2), ncol=2), widths=c(2,1), height=4)
  
  # Create the plot of Europe with hydrobasins
  par(mar = c(5, 3, 0, 0))
  
  # Set up plotting background
  nearctic.plot <- ne_countries(scale = "Large", geounit = c('United States of America', 'Canada', "Mexico"), type = 'sovereignty',
                                returnclass = c("sp", "sf"))
  ymin = 28
  ymax = 50
  xmin= -130
  xmax= -62
  
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(nearctic.plot, realm_bbox), col = 'lightgrey', border = rgb(0,0,0,0),
       axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  
  ## Add hydrobasin statistics
  
  # Extract desired level of hydrobasins
  lay_EU_6 <- sort(ogrListLayers("/Users/sassen/Desktop/Hydrobasins/hybas_na_lev01-12_v1c"))[HYBAS_level]
  
  # Extract geometry of basins
  basin.shapes_EU <- st_read(dsn = "/Users/sassen/Desktop/Hydrobasins/hybas_na_lev01-12_v1c", 
                             layer = lay_EU_6)
  colname_level <- paste0('HYBAS_ID_', HYBAS_level)
  # Compute summary statistics per basin
  
  temp<-full.novel.mat.season |>
    #filter civ for now, slight issue with hybas scheme
    dplyr::filter(country != "CIV") |>
    left_join(HYBAS_scheme[,c('TimeSeriesID',colname_level)], by = c('site_ID' = 'TimeSeriesID')) |>
    group_by_at(vars(starts_with('HYBAS_ID_'),'novel')) |>
    summarise(observed = n()) |>
    st_drop_geometry()
  
  # Compute, for each basin, the proportion of timeseries where novelty occured, 
  # versus those that did not.
  non.nov <- NULL
  nov <- NULL
  prop.nov <- NULL
  
  for (basin in unique(temp[[colname_level]])){
    print(basin)
    basin.sub <- temp[which(temp[1] == basin),]
    # Convoluted, but essentially divides novelty by non-novelty at time series level
    if(nrow(basin.sub) > 1){
      non.nov <- c(basin.sub[which(basin.sub$novel == 0),'observed'], non.nov)
      nov <- c(ifelse(nrow(basin.sub[which(basin.sub$novel == 1),'observed'])==0,0,
                      basin.sub[which(basin.sub$novel == 1),'observed']), nov)
      
      # Continue here... wip
      prop.nov <- c(ifelse(basin.sub[which(basin.sub$novel == 0),'observed']+
                             basin.sub[which(basin.sub$novel == 1),'observed'] > 0,  
                           basin.sub[which(basin.sub$novel == 1),'observed']/basin.sub[which(basin.sub$novel == 0),'observed'], 
                           NA), prop.nov)
    }else{
      # if no novel time series in basin, just add 0
      non.nov <- c(basin.sub$observed, non.nov)
      nov <- c(0, nov)
      prop.nov <- c(0, prop.nov)
    }
    
  }
  # Store in a new df upon loop completion
  prop.nov.ByBasin <- data.frame('basin' =unique(temp[1]),
                                 'proportion_novel' = unlist(prop.nov),
                                 'novel' = unlist(nov),
                                 'non_novel' = unlist(non.nov))
  
  # Create a colour ramp for plotting
  rbp<-colorRampPalette(c('gray60','orange'))
  
  # Add colour ramp
  colour.vec<- rbp(5)
  prop.nov.ByBasin <- prop.nov.ByBasin |>
    mutate(colramp = ifelse(proportion_novel == 0, colour.vec[1],
                            ifelse(proportion_novel > 0 & proportion_novel <= 0.05,colour.vec[2],
                                   ifelse(proportion_novel <= 0.10 & proportion_novel > 0.05 ,colour.vec[3],
                                          ifelse(proportion_novel <= 0.15 & proportion_novel > 0.10,colour.vec[4],colour.vec[5])))))
  # Isolate sampled basins for plotting
  filtered.basin.shapes <- basin.shapes_EU[basin.shapes_EU$HYBAS_ID %in% HYBAS_scheme[,colname_level],]
  
  plotting.frame<-left_join(filtered.basin.shapes, prop.nov.ByBasin, by = c('HYBAS_ID'= colname_level )) |>
    dplyr::filter(!is.na(proportion_novel))
  
  plot(plotting.frame$geometry,
       col = plotting.frame$colramp,
       border=rgb(0,0,0,0),
       add =T)
  
  # Next bit
  
  ymin.map = par('usr')[3] 
  ymax.map = par('usr')[4]
  # define the bin width
  bin_width <- 2
  
  df <- full.novel.mat.season |>
    dplyr::filter(BioRealm == 'Nearctic')
  df_nov <- subset(df, novel == 1)
  # create a new column 'lat_group' by dividing latitude by the bin width and rounding down
  df$lat_group <- floor(df$Latitude / bin_width) * bin_width
  df_nov$lat_group <- floor(df_nov$Latitude / bin_width) * bin_width
  
  # group the dataframe by 'lat_group' and count the number of rows
  df_counts <- aggregate(df$Latitude, by = list(df$lat_group), length)
  
  df_counts_nov <- aggregate(df_nov$Latitude, by = list(df_nov$lat_group), length)
  
  # create a new dataframe with latitude range and number of observations
  range_df <- data.frame(lat_range = paste(ifelse(df_counts$Group.1 >= 0, paste0(df_counts$Group.1, 'N'), paste0(abs(df_counts$Group.1), 'S')),
                                           ifelse(df_counts$Group.1 + bin_width - 1 >= 0, paste0(df_counts$Group.1 + bin_width - 1, 'N'), paste0(abs(df_counts$Group.1 + bin_width - 1), 'S')),
                                           sep = '-'),
                         no_obs = df_counts$x,
                         min_lat = df_counts$Group.1) |>
    dplyr::filter(min_lat >= ymin)
  
  range_nov_df <- data.frame(lat_range = paste(ifelse(df_counts_nov$Group.1 >= 0, paste0(df_counts_nov$Group.1, 'N'), paste0(abs(df_counts_nov$Group.1), 'S')),
                                               ifelse(df_counts_nov$Group.1 + bin_width - 1 >= 0, paste0(df_counts_nov$Group.1 + bin_width - 1, 'N'), paste0(abs(df_counts_nov$Group.1 + bin_width - 1), 'S')),
                                               sep = '-'),
                             no_obs = df_counts_nov$x,
                             min_lat = df_counts_nov$Group.1) |>
    dplyr::filter(min_lat >= ymin)
  
  # Sometimes the nov.df and df are different lengths. we need to account for this.
 range_nov_df <- left_join(range_df, range_nov_df, by= c('lat_range')) |>
    mutate_at(c('no_obs.y'), ~replace_na(.,0)) |>
    mutate(no_obs = no_obs.y,
           min_lat = min_lat.x) |>
    dplyr::select(lat_range, no_obs, min_lat)
  
  # Create the plot of Europe with hydrobasins
  par(mar = c(5, 0, 0, 1))
  
  range_df$propnovel <- range_nov_df$no_obs/range_df$no_obs
  
  plot(NULL, xlim=c(0,max(range_df$propnovel)), 
       ylim=c(ymin.map,ymax.map-1.3), axes= F, xlab = "")
  axis(side=1)
  mtext("Proportion of Novel Communities", 
        side=1,line=3, cex.lab=1)
  
  for(i in 1:(nrow(range_df)-1)){
    print(i)
    polygon(y=c(range_df$min_lat[i], range_df$min_lat[i+1],range_df$min_lat[i+1] ,range_df$min_lat[i] ),
            x= c(0,0,range_df$propnovel[i], range_df$propnovel[i]), col='grey70')
    # polygon(y=c(range_df$min_lat[i], range_df$min_lat[i+1],range_df$min_lat[i+1] ,range_df$min_lat[i] ),
    #      x= c(0,0,range_df$no_obs[i], range_df$no_obs[i]), col = "grey70")
    #polygon(y=c(range_nov_df$min_lat[i], range_nov_df$min_lat[i+1],range_nov_df$min_lat[i+1] ,range_nov_df$min_lat[i] ),
    #       x= c(0,0,range_nov_df$no_obs[i], range_nov_df$no_obs[i]), col = "orange")
    
    
  }
  
}

venn_plot_main <- function(full.novel.mat){
  
  rand.preds <- lapply(1:3, function(n){
    success_cat = c("instant", "cumul", "novel")[n]
    
    full.novel.mat$success = full.novel.mat[, success_cat]
    mod<- glmer(success~1+bin_lag+position + (1|Quarter/site_ID), data = full.novel.mat, family = 'binomial')
    
    pred.df <- as.data.frame(summary(mod)$coefficients[1,])
    pred.df$taxa.rand <- summary(mod)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(mod)$varcor$length[1,1]
    return(pred.df)
  })
  
  names(rand.preds) <- c("instant", "cumul", "novel")
  
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.43,0.65,0.51)
  circle.radius <- c(0.2, 0.22)
  circle.y <- 0.5
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  par(lheight=0.85)
  text(x=circle.cent[3]+0.02, y=circle.y+0.35, 
       adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=1.4)
  
  par(xpd=FALSE)
  # turn I circle to spatial polygons to get overlap
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  N.sp <- raster::intersect(I.sp, C.sp)
  I.sub <- gDifference(I.sp, C.sp)
  C.sub <- gDifference(C.sp, I.sp)
  
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border="black", lwd=1)
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=1.5)
  
  text(x=0.225, y=circle.y+0.2,
       labels=c("Instantaneous\nnovelty (I)    "), 
       col="red", adj=1, cex=1.4)
  text(x=0.82, y=circle.y+0.2,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=1.4)
  
  overall.means <- t(sapply(rand.preds, function(x){
    x <- t(x)
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.073, 0.081, 0.02)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.8)
  
  text(x = c(circle.cent + c(-0.076, 0.075, 0.015)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=1.7)
}

pdf(file = "/Users/sassen/Desktop/Figure_test_EU.pdf",
    width = 15,
    height = 12.5)
complete.basin.novelty.plotter.PAL(full.novel.mat.season, HYBAS_scheme, 8)
dev.off()

pdf(file = "/Users/sassen/Desktop/Figure_test_NA.pdf",
    width = 25,
    height = 8)
complete.basin.novelty.plotter.NEA(full.novel.mat.season, HYBAS_scheme, 7)
dev.off()

pdf(file = "/Users/sassen/Desktop/Figure_test_venn.pdf",
    width = 15,
    height = 15)
venn_plot_main(full.novel.mat.season)
dev.off()

