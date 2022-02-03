# Find the attributes of the Basins and link them to basin name.

library(rgdal)
library(raster)
library(rgeos)
library(sf)
library(tidyverse)
library(rfishbase)

### Matching basin names and ID using Geolocation, so that we can determine invasives per basin ####

# Load in the shapefile which contains basin names AND coordinate-based shapefiles (Tedesco et Al.)
occurence_shapefile <- read_sf(dsn = "/Users/sassen/Desktop/datatoFigshare/Basin042017_3119.shp")

# We now need to allocate a basin name for every TimeSeries_ID in time_series_data,
# so that we can reliably tag our fish species.

# This function matches available names to HydroBasin ID's
basin_name_code <- basin_name_match_function(occurence_shapefile)

# Let's add a column to the original timeseries data which will contain basin name.

time_series_data$Basin_name <- NA

for (i in 1:nrow(time_series_data)) {
  print(paste0("Working on TimeSeries ", i, " out of ", nrow(time_series_data)))
  for (j in 1:nrow(basin_name_code)) {
    if(time_series_data$HydroBasin[i] == basin_name_code$HydroBasin[j]){
      
      time_series_data$Basin_name[i] <- basin_name_code$Basin_name[j]
    }
  }
}




# Now we have Tedesco et al names for the majority of the basins we're interested in!
# We can now go on and tag the fish species based on the Tedesco database!

# However, the Tedesco database actually is lacking quite a few species in each  basin
# (I don't know why). The only solution I can come up with for this is to make a "backup-list"
# with the exotic/native status of all species WITHIN A COUNTRY, based on the fishbase status.
# So, the procedure for determining fish status will become -> cross-reference fish species in a
# particular basin using the Tedesco database -> if there is no status for the species in that basin
# -> descend to the country level and assign status based on that. I have created separate files
# for each country and stored these in the "exotic_databases_countries" file in the R-project.


# It is nice that we have found novelty. Novelty, in our context, is rapid and substantive 
# change of the identity and abundance of species within a certain locality (i.e. changes in 
# the assemblage structure). We thus know that communities that we have classified as 
# novel using our framework are by default different to previous communities at that 
# locality. So, comparing stability of communities where novelty occurred versus those where it 
# did not seems trivial, because this instability/rapid change is what makes them novel in the first 
# place. More interesting would be looking into instability in terms of which fish species 
# are actually undergoing major population fluctuations. Are they native ones? Or are they exotic? 
# Perhaps "novelty" at the ecological timescale (decades) is simply a breakdown of synchrony (i.e. 
# population fluctuations of different species are positively correlated) ? We could stretch this train
# of thought even further; what is the average trophic level of the population, and is novelty at all 
# associated with changes in this metric? In my mind, novel community states should be associated with 
# a reshuffling of the ecological interactions within that community. Is this reflected in the population
# fluctuations of different fish species (grouped into feeding guilds, trophic levels, etc.)? 

# Plot species timeseries 

df <- NULL
y <- (Fish_Communities_A[["BioRealm_Matrices_A"]][["palearctic_mat_A"]][["G135"]])
for(i in 1:ncol(y)){
temp_df <- data.frame(x=1:nrow(y), y=y[,i], col=rep(i:i, each=nrow(y)))
df <- rbind(df,temp_df)} 

ggplot(df,aes(x=x,y=log(y),group=col,colour=factor(col))) + geom_line() # plot data








