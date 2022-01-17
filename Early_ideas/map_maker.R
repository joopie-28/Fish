library("rnaturalearth")
library("rnaturalearthdata")
library(data.table)

# Creating a dataframe with the timeseries where we found novelty
# This is basically the set-up

community <- GLM_lists_B$GLM_input_B_2

community <- subset(community,  bins <50)

novel_points <- by(community$novel, INDICES = community$site, FUN = sum)

novel_points <- as.list(novel_points)

df = data.frame(matrix(vector(), length(novel_points), 2,
                       dimnames=list(c(), c("TimeSeriesID", "Novelty"))),
                stringsAsFactors=F)

df$TimeSeriesID <- names(novel_points)
df$Novelty <- novel_points

df <- subset(df, Novelty > 0)       

# These are novel based on abundance A2 bins
novelty_subset <- subset(time_series_data, TimeSeriesID %in% df$TimeSeriesID)

# These are not novel based on abundance A2 bins
placeholder <- subset(time_series_data, TimeSeriesID %in% community$site)

non_novelty_subset <- subset(placeholder, !(TimeSeriesID %in% df$TimeSeriesID))


# Create a map with colour-coded timeseries locations!
# There's a couple variations here

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) + xlim(c(-120, 30)) + ylim(c(30, 70)) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Novel community emergence in Europe and the USA", subtitle = "between 1970 and 1980 (bin size = 2 years, data = presence/absence)") +
  geom_point(data = non_novelty_subset, aes(x = non_novelty_subset$Longitude, y = non_novelty_subset$Latitude), size = 2, 
             shape = 23, fill = "darkred") +
  geom_point(data = novelty_subset, aes(x = novelty_subset$Longitude, y = novelty_subset$Latitude), size = 2, 
             shape = 23, fill = "gold") 

ggplot(data = world) + xlim(c(-120, 30)) + ylim(c(30, 70)) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Novel community emergence in Europe and the USA", subtitle = "between 1970 and 1980 (bin size = 2 years, data = presence/absence)") +
  geom_point(data = non_novelty_subset, aes(x = non_novelty_subset$Longitude, y = non_novelty_subset$Latitude), size = 2, 
             shape = 23, fill = "darkred") +
  geom_point(data = novelty_subset, aes(x = novelty_subset$Longitude, y = novelty_subset$Latitude), size = 2, 
             shape = 23, fill = "gold") 


# Some further analyses, a histogram of timeseries and novelty per country

legend <- c("Novel" = "gold", "Total" = "grey")

ggplot(placeholder, aes(x = Country)) + ylab("Number of TimeSeries") +
  geom_histogram(stat = "count", color = "grey", fill = "grey") + 
    geom_histogram(data = subset(time_series_data, TimeSeriesID %in% df$TimeSeriesID),
                 stat = "count", color = "gold", fill = "gold") +
  scale_color_grey() + theme_classic() + labs(title = "TimeSeries country distribution (Abundance)") 

  
  


  

  
  







