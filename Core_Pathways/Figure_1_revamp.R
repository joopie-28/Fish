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
  
  plot(crop(palearctic.plot, realm_bbox), col = 'gray60', border = rgb(0,0,0,0),
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
  rbp<-colorRampPalette(c('lightgrey','orange'))
  
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
    height = 18)
plot.basin.summary(full.novel.mat.season, HYBAS_scheme, 8)
dev.off()



