# River maps

library(maps)
library(rworldmap) # to plot axes
library(rgdal) # to load the shapefile
library(rworldxtra)

newmap <- getMap(resolution = "high")

# Remember to extract RAR files before using this function
riversData <- readOGR(dsn = "/Users/sassen/Desktop/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp", "HydroRIVERS_v10")

lat1 <- runif(5,  35, 55)
long1 <- runif(5,  0, 50)
lat2 <- runif(5,  35, 55)
long2 <- runif(5,  0, 50)

svg("map_Europe.svg", height=4, width=6)
par(mar=c(3, 3, 2, 2))
waterColor <- 'cadetblue1'
plot(newmap, xlim = c(-5, 9), ylim = c(42, 51), 
     asp = 1,lty=2, lwd=1,
     bg=waterColor, col="#ffffff",)

plot(riversData, col=waterColor, add=T) # plot rivers
map.axes()
points(x=long1, y=lat1, pch= 15, col="black", cex=1.2)
points(x=long2, y=lat2, pch= 17, col="red", cex=1.2)
legend("topright", bg="white", pt.cex = 1.2,
       legend = c('group1', 'group2'), col=c('black', 'red'), pch = c(15, 17))
