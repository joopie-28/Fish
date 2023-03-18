#### HydroRivers Alternative Analyses ####

# This script handles all importing and pre-processing of the HydroRivers database
# It allows us to extract data from millions of reaches, streams and rivers around 
# the globe, and match these with our fish surveys.

RiverID_eu <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")
RiverID_aus <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_au_shp/HydroRIVERS_v10_au_shp")
RiverID_af <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af_shp")
RiverID_sa <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_sa_shp/HydroRIVERS_v10_sa_shp")
RiverID_na <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp")
RiverID_as <- st_read(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as_shp")

# Function for extracting River ID's from data

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


# Now we need to import all the environmental variables.

AUSRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_au")
EURiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                            lay = "RiverATLAS_v10_eu")
NARiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                            lay = "RiverATLAS_v10_na")
SANRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                            lay = "RiverATLAS_v10_sa_north")
SASRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                         lay = "RiverATLAS_v10_sa_south")
AFRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                            lay = "RiverATLAS_v10_af")
ASIRiver_Data <- st_read(dsn ="/Users/sassen/Desktop/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp",
                            lay = "RiverATLAS_v10_as")
   
environmental_variables_river = data.frame('Variable' = c("Natural_Discharge_Annual",
                                                    "Land_Surface_Runoff_Annual",
                                                    "Degree_Regulation",
                                                    "River_Area",
                                                    "Elevation",
                                                    "Terrain_slope",
                                                    "Stream_Gradient",
                                                    "Aridity_Index",
                                                    "Cropland_Extent",
                                                    "Pasture_Extent",
                                                    "Protected_Area_Extent",
                                                    "Population_Density",
                                                    "Urban_Extent",
                                                    "Human_Footprint"),
                                           'Code' = c('dis_m3_pyr',
                                                "run_mm_cyr",
                                                "dor_pc_pva",
                                                "ria_ha_csu",
                                                "ele_mt_cav",
                                                "slp_dg_cav",
                                                "sgr_dk_rav",
                                                "ari_ix_cav",
                                                "crp_pc_cse",
                                                "pst_pc_cse",
                                                "pac_pc_cse",
                                                "ppd_pk_cav",
                                                "urb_pc_cse",
                                                "hft_ix_c09")) 

# Create this massive list with all environmental variables for all possible rivers.
data_list <- list(AUSRiver_Data, EURiver_Data, NARiver_Data , 
     SANRiver_Data, SASRiver_Data, AFRiver_Data,
     ASIRiver_Data)

# Filter out rows and variables to create a workable amount of data
for (i in 1:length(data_list)){
  print(i)
  data_list[[i]] <- data_list[[i]] |> 
    select(c(1:14,environmental_variables_river$Code)) |>
    filter(HYRIV_ID %in% globalRivers_Extracted$HYRIV_ID)
}

# Collate into a large dataframe
RiverData_Global <- rbindlist(data_list)
rm(data_list)

# Save it to the environemnt so we do not have to recompute these data.
saveRDS(RiverData_Global, file = './outputs/RiverData_Global.rds')

# Join all data, such that we have environmnetal variables
# linked with a timeseries ID.
MergedENV_by_RivID <- globalRivers_Extracted |>
  left_join(RiverData_Global, by = 'HYRIV_ID') |>
  rename(ID = site)

# Join all environmental data to our time series 
# data frame; at both site and community levels.

# This is the site-level data frame
FullGeoFrame <- geo.timeseries.full |>
  left_join(MergedENV_by_RivID, by = ("ID")) |>
  mutate(binary_novel = ifelse(novel > 0, 1, 0))|>
  separate(ID, c('Site', 'Quarter'), remove = F)

# This is the community-level data frame
FullEnvFrame <- fullNovFrame_complete |>
  rename(ID = site) |>
  left_join(MergedENV_by_RivID, by = c('ID'))

# This frame includes only novel communities and is used
# for persistence length modelling
PersistenceFrame <- FullEnvFrame|>
  filter(novel.class == 'Persister' | novel.class == 'BLIP') |>
  filter(consistency == T) |>
  mutate(binary_pers = ifelse(novel.class == "Persister", 1, 0))


#### End of script ####

plot(ne_countries(scale = "Large", country = 'australia', type = 'countries',
                  returnclass = c("sp", "sf")), axes = T, xlim = c(152,154), 
     ylim = c(-28,-26))
points(FullEnvFrame$Longitude, FullEnvFrame$Latitude, pch = 19, col = 'black',
     cex = 0.4)
plot(FullEnvFrame$geometry, col = as.factor(FullEnvFrame$basin), add = T)

summary(glm(binary_novel~1, data= FullEnvFrame, family = 'binomial'))

test <- FullEnvFrame %>% 
  mutate_at(c('position','bin_lag',"gain", "loss" ,"run_mm_cyr", "INC_increase","total_inv", "DIST_DN_KM","DIST_UP_KM","dis_m3_pyr",'ria_ha_csu',"ele_mt_cav","urb_pc_cse", "crp_pc_cse", "pac_pc_cse", "pst_pc_cse"), 
                           ~(scale(., center =T, scale =T) %>% as.vector)) |>
  mutate(True_Novel = ifelse(novel == 1 & consistency == TRUE, 1, 0))
mod1 = glm(True_Novel~position+bin_lag+gain+loss+evenness+ INC_increase*total_inv+dis_m3_pyr+dis_m3_pyr*urb_pc_cse+ele_mt_cav+urb_pc_cse+DIST_DN_KM,data=(test), family = 'binomial')

plot_model(mod1, vline.color = "red")


#31:38
