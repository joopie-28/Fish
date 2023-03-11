###### Script for extracting HYDROBASIN ID codes from coordinates #########

create_basin_TS <- function(time_series_data){
  
  basin_list <- lapply(1:12, function(Hydro_Level){
    print(paste0("Extracting Basin Code at level", Hydro_Level))
    # Extract data at level of choice
    lay_EU = sort(ogrListLayers("/Users/sassen/Desktop/hybas_eu_lev01-12_v1c"))[Hydro_Level]
    lay_NA = sort(ogrListLayers("/Users/sassen/Desktop/hybas_na_lev01-12_v1c"))[Hydro_Level]
    lay_AUS = sort(ogrListLayers("/Users/sassen/Desktop/hybas_au_lev01-12_v1c"))[Hydro_Level]
    lay_ASI = sort(ogrListLayers("/Users/sassen/Desktop/hybas_as_lev01-12_v1c"))[Hydro_Level]
    lay_SA = sort(ogrListLayers("/Users/sassen/Desktop/hybas_sa_lev01-12_v1c"))[Hydro_Level] 
    lay_AF = sort(ogrListLayers("/Users/sassen/Desktop/hybas_af_lev01-12_v1c"))[Hydro_Level]
    
    # Combine all regions
    world_Basin_Level <- rbind(st_read(dsn = "/Users/sassen/Desktop/hybas_eu_lev01-12_v1c", 
                                       layer = lay_EU),
                               st_read(dsn = "/Users/sassen/Desktop/hybas_na_lev01-12_v1c", 
                                       layer = lay_NA), 
                               st_read(dsn = "/Users/sassen/Desktop/hybas_au_lev01-12_v1c", 
                                       layer = lay_AUS),
                               st_read(dsn = "/Users/sassen/Desktop/hybas_as_lev01-12_v1c", 
                                       layer = lay_ASI),
                               st_read(dsn = "/Users/sassen/Desktop/hybas_sa_lev01-12_v1c", 
                                       layer = lay_SA),
                               st_read(dsn = "/Users/sassen/Desktop/hybas_af_lev01-12_v1c", 
                                       layer = lay_AF))
    
    # Create a new spatial data frame with hydrobasin ID codes
   
    # planar geometries
    sf_use_s2(FALSE) 
    
    SpatialTimeSeries <- time_series_data |>
      # filter those with non-used time series
      filter(TimeSeriesID %in% full.novel.mat.season$site_ID) |>
      # convert to spatial object
      st_as_sf(coords = c("Longitude", "Latitude"), crs = WG84) |>
      # find intersections with european hydrobasins
      st_intersection(world_Basin_Level) |>
      # select relevant columns
      dplyr::select('TimeSeriesID',
             'HYBAS_ID') |>
      # Tidyverse needs some odd syntax for renaming using paste0()
      rename(!! paste0('HYBAS_ID_', Hydro_Level) := "HYBAS_ID")
    TimeSeries <- st_drop_geometry(SpatialTimeSeries)
    return(TimeSeries)
  })
  
  return(Reduce(left_join,basin_list))
  
}


