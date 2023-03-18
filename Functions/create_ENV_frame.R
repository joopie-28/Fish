# Function for creating a data frame with variables extracted from the HydroAtlas

create_ENV_frame <- function(df, HYBAS_Level, HYBAS_scheme, variables_list){
  
  lay = ogrListLayers("/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb")[HYBAS_Level]
  
  # test
  HydroAtlas_lev <- st_read(dsn = "/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb", 
                             layer = lay)
  
  HydroAtlas_lev$HYBAS_ID <- as.factor(HydroAtlas_lev$HYBAS_ID)
  
  
  HYBAS_sub <- HYBAS_scheme[, c(1, which(colnames(HYBAS_scheme) == paste0('HYBAS_ID_', HYBAS_Level)))]
  HYBAS_sub[, 2] <- as.factor(HYBAS_sub[,2])
  class.frame <- df|> 
    rename('HYBAS_ID' = 'Basin') |> 
    left_join(HYBAS_sub, by = c("Site" = "TimeSeriesID")) |>
    dplyr::select(!"HYBAS_ID") |>
    rename("HYBAS_ID" = paste0('HYBAS_ID_', HYBAS_Level)) |>
    
    left_join(as.data.frame(HydroAtlas_lev[,c('HYBAS_ID',variables_list)]), by = 'HYBAS_ID')

  return(class.frame)
}



