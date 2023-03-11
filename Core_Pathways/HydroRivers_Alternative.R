#### HydroRivers Alternative Analyses ####


lay = ogrListLayers("/Users/sassen/Desktop/RiverATLAS_Data_v10.gdb/RiverATLAS_v10.gdb")

# test
RiverAtlas <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10.gdb/RiverATLAS_v10.gdb")

lay = ogrListLayers("/Users/sassen/Desktop/HydroRIVERS_v10.gdb/HydroRIVERS_v10.gdb")

# test
RiverID_eu <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")

plot(ne_countries(scale = "Large", country = 'spain', type = 'countries',
             returnclass = c("sp", "sf")), axes = T, xlim = c(-5,0), ylim = c(42.5,44))

plot(RiverID_eu[which(RiverID_eu$HYBAS_L12 %in% EnviroByTS_L12$HYBAS_ID[110:120]),'HYRIV_ID'],
     add = T)
plot(st_buffer(EnviroByTS_L12$geometry[110:120],150), col = 'black', 
     pch = 19, cex = 15, add = T)

# need to buffer
View(st_intersects(st_buffer(test, 1),st_buffer(EnviroByTS_L12$geometry[15], 150)))


enviro_aus <- subset(EnviroByTS_L12, BioRealm == 'Australasia')


test<-RiverID_aus[which(RiverID_aus$HYBAS_L12 %in% enviro_aus$HYBAS_ID),c('HYRIV_ID','geometry')]
                 


RiverID_eu <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")
RiverID_aus <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_au_shp/HydroRIVERS_v10_au_shp")
RiverID_af <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af_shp")
RiverID_sa <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_sa_shp/HydroRIVERS_v10_sa_shp")
RiverID_na <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp")
RiverID_as <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as_shp")
# write the function 

ExtractRiverFromSurvey <- function(EnviroByTS_L12, country_input, RiverID){
  
  # Subset data based on biorealm to reduce computing time
  enviro_sub <- subset(EnviroByTS_L12, Country == country_input)
  RiverID_Sub<-RiverID[which(RiverID$HYBAS_L12 %in% enviro_sub$HYBAS_ID),
                    c('HYRIV_ID','geometry')]
  
  print('Extracting river ID from point data')
  # Find closest line to point
  output <- as.data.frame(geosphere::dist2Line(p=st_coordinates(enviro_sub$geometry), 
                                                 line = as_Spatial(RiverID_Sub$geometry))) 
  
  output$HYRIV_ID <- NA
  output$site <- NA
  for(i in 1:nrow(output)){
    output$HYRIV_ID[i] <- RiverID_Sub$HYRIV_ID[output$ID[i]]
    output$site[i] <- enviro_sub$ID[i]
  }     
  
  output <- output |>
    select(HYRIV_ID, site)
  
 return(output)
}

RiverID_Global <- rbindlist(list(RiverID_eu, RiverID_af, RiverID_as, 
                                 RiverID_aus, RiverID_na, RiverID_sa))

globalRivers_Extracted <- rbindlist(lapply(countries.suf.data, function(country){
  print(country)
  ExtractRiverFromSurvey(EnviroByTS_L12, country, RiverID_Global)
}))

# This was a very intensive computation, due to the size of the HydroRiver Database,
# so we will be saving the obejct in an .RDS file.

saveRDS(globalRivers_Extracted, file = './outputs/globalRivers_Extracted.rds')

####pseudocode
# Then we need to extract environemntal variables
Rivers_ENV <- #read table from qgis
  |> select('relevant columns')
# should be as easy as 

globalRivers_Extracted |>
  left_join(Rivers_ENV, by = 'HYRIV_ID')

# And then, simply join on geo.timeseries.full
                     
                     
                     
