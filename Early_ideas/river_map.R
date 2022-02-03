# River maps

library(maps)
library(rworldmap) # to plot axes
library(rgdal) # to load the shapefile
library(rworldxtra)

newmap <- getMap(resolution = "low")

# Remember to extract RAR files before using this function
riversData <- readOGR(dsn = "/Users/sassen/Desktop/HydroRIVERS_v10_shp/HydroRIVERS_v10_shp", "HydroRIVERS_v10")

basin_geo_data <- readOGR(dsn = "/Users/sassen/Desktop/eu_bas_15s_beta", "eu_bas_15s_beta")

lat1 <- runif(5,  35, 55)
long1 <- runif(5,  0, 50)
lat2 <- runif(5,  35, 55)
long2 <- runif(5,  0, 50)

svg("map_Europe.svg", height=4, width=6)
par(mar=c(3, 3, 2, 2))
waterColor <- 'cadetblue1'
plot(newmap, xlim = c(0, 5), ylim = c(40, 49), 
     asp = 1,lty=1, lwd=1,
     bg= "white", col="#ffffff", axes = T)

plot(hydrobasins_12[c(29813:30231), 4], add =TRUE, col = "blue") # 590
plot(hydrobasins_12[c(30231:30413), 4], add = TRUE, col = "red") # 620






plot(riversData, col= "grey" , add=T) # plot basins

plot(flow_data, col = waterColor, add=T) # plot rivers
plot(hydrobasins_12[c(45000:46000), 1], add =TRUE)
dev.off()

map.axes()
points(x=long1, y=lat1, pch= 15, col="black", cex=1.2)
points(x=long2, y=lat2, pch= 17, col="red", cex=1.2)
legend("topright", bg="white", pt.cex = 1.2,
       legend = c('group1', 'group2'), col=c('black', 'red'), pch = c(15, 17))
