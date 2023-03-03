# Environmental Variables Hydrosheds integration

# These functions extract the environmental variables from
# hydrosheds at a given HYBAS LEVEL.

library(rgdal)
library(sf)
library(tidyverse)
library(lme4)

create_ENV_frame <- function(full.novel.mat.season, HYBAS_Level, HYBAS_scheme){
  
  lay = ogrListLayers("/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb")[HYBAS_Level]
  
  # test
  HydroAtlas_lev <- st_read(dsn = "/Users/sassen/Desktop/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb", 
                             layer = lay)
  
  HydroAtlas_lev$HYBAS_ID <- as.factor(HydroAtlas_lev$HYBAS_ID)
  
  
  HYBAS_sub <- HYBAS_scheme[, c(1, which(colnames(HYBAS_scheme) == paste0('HYBAS_ID_', HYBAS_Level)))]
  HYBAS_sub[, 2] <- as.factor(HYBAS_sub[,2])
  class.frame <- full.novel.mat.season |> 
    rename('HYBAS_ID' = 'basin') |> 
    left_join(HYBAS_sub, by = c("site_ID" = "TimeSeriesID")) |>
    select(!"HYBAS_ID") |>
    rename("HYBAS_ID" = paste0('HYBAS_ID_', HYBAS_Level)) |>
    
    left_join(as.data.frame(HydroAtlas_lev[,c('HYBAS_ID','run_mm_syr',
                                               'dis_m3_pyr', 'riv_tc_ssu',
                                               'dor_pc_pva',
                                               "crp_pc_sse", 'pst_pc_sse',
                                               'pac_pc_sse', 'hft_ix_s09','ppd_pk_sav')]), by = 'HYBAS_ID')

  return(class.frame)
}


class.frame <- create_ENV_frame(full.novel.mat.season, 12, HYBAS_scheme)






summary(glmer(cbind(novel, length-novel)~run_mm_syr +dis_m3_pyr +pst_pc_sse + (1|HYBAS_ID), 
            data = class.frame, family = 'binomial'))

