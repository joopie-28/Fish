# Function for matching survey point data with river ID from the hydrorivers database.

# Some survey points are slightly off the closest river, so we match these points with
# the closest river 

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