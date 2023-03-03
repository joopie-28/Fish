# Plotting functions figure 1.

map.realm.plotter <- function(full.novel.mat.season){
  
  par(mfrow=c(2,2))
  col <- subset(full.novel.mat.season, cat =='novel')
  col.df<-geo.timeseries[which(geo.timeseries$TimeSeriesID %in% col$site_ID),]
  non <- subset(full.novel.mat.season, cat !='novel')
  non.df<-geo.timeseries[which(geo.timeseries$TimeSeriesID %in% non$site_ID),]
  # Import geodata from rnaturalearth package and set cropping parameters
  
  palearctic.plot <- ne_countries(scale = "Large", continent = 'Europe', type = 'sovereignty',
                                  returnclass = c("sp", "sf"))
  ymin = 40
  ymax = 70
  xmin= -12
  xmax= 50
  
  text.offset = c(-1.95, 1.8,7)
  fig.offset = c('y.offset' = 0.5, 'cex1' = 1, 'cex2'=0.8, 'circl.cent1' = 0.77, 'diam.large'=0.18, 'diam.small'=0.17)
  par(mar=c(0,2,1,1))
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(palearctic.plot, realm_bbox), col = 'grey97', axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  text("a)", y=69,x=-9)
  points(x=non.df$Longitude, y=non.df$Latitude, pch = 19, cex=0.75, col = 'black')
  points(x=col.df$Longitude, y=col.df$Latitude, pch = 19, cex=0.75, col = 'orange')
  
  
  nearctic.plot <- ne_countries(scale = "Large", geounit = c('United States of America', 'Canada', "Mexico"), type = 'sovereignty',
                                returnclass = c("sp", "sf"))
  ymin = 25
  ymax = 70
  xmin= -130
  xmax= -55
  
  text.offset = c(-1.4, 2.1,7.6)
  fig.offset = c('y.offset' = 1, 'cex1' = 0.9, 'cex2'=0.7, 'circl.cent1' = 0.75,'diam.large'=0.15, 'diam.small'=0.145)
  
  par(mar=c(0,1,1,0))
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(nearctic.plot, realm_bbox), col = 'grey97', axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  text("b)", y = 69, x=-128)
  
  points(x=non.df$Longitude, y=non.df$Latitude, pch = 19, cex=0.75, col = 'black')
  points(x=col.df$Longitude, y=col.df$Latitude, pch = 19, cex=0.75, col = 'orange')
  
  
  aus.plot <- ne_countries(scale = "Large", geounit = c('australia', 'new zealand', 'papua new guinea', 'indonesia'), type = 'map_units',
                           returnclass = c("sp", "sf"))
  ymin = -47 
  ymax = 0 
  xmin = 110 
  xmax = 180
  
  text.offset = c(-2.5, 3, 7.5)
  fig.offset = c('y.offset' = 1.5, 'cex1' = 1, 'cex2'=0.8, 'circl.cent1' = 0.74, 'diam.large'=0.18, 'diam.small'=0.17)
  par(mar=c(0,2,0,1))
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(aus.plot, realm_bbox), col = 'grey97', axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  text("c)", x=112, y=-2)
  
  segments(x0=149, y0=-29, x1 =149 , y1 =-25 )
  segments(x0=149, y0=-29, x1 =155 , y1 =-29 )
  segments(x0=155, y0=-29, x1 =155 , y1 =-25 )
  segments(x0=155, y0=-25, x1 =149 , y1 =-25 )
  segments(x0 = 155, y0 = -29, x1=135.2, y1= -47.6   )
  segments(x0 = 149, y0 = -25, x1=,116, y1= -29   )
  
  points(x=non.df$Longitude, y=non.df$Latitude, pch = 19, cex=0.75, col = 'black')
  points(x=col.df$Longitude, y=col.df$Latitude, pch = 19, cex=0.75, col = 'orange')
  par(mar=c(0,1,0,2))
  
  venn_plot_main(full.novel.mat.season)
  text("d)", x=.05,y=0.93)
  
  par(fig = c(0.03, 0.2, 0.02, 0.22), new = T, mai= c(0,0,0,0)) 
  qld_bbox = st_bbox(c(ymin = -29, ymax = -25, xmin =149, xmax = 160))
  plot(crop(ne_countries(scale = 'Large', country = c('Australia'), type = 'sovereignty',
                         returnclass = c("sp", "sf")), qld_bbox), col = 'grey97', bg= 'white')
  
  points(x=non.df$Longitude, y=non.df$Latitude, pch = 19, cex=0.75, col = 'black')
  points(x=col.df$Longitude, y=col.df$Latitude, pch = 19, cex=0.75, col = 'orange')
  box(which = 'plot',  lty = 1)
  
  
  
}
venn_plot_main <- function(full.novel.mat){
  
  rand.preds <- lapply(1:3, function(n){
    success_cat = c("instant", "cumul", "novel")[n]
    
    full.novel.mat$success = full.novel.mat[, success_cat]
    mod<- glmer(success~1+bin_lag+position + (1|Quarter/site_ID), data = full.novel.mat, family = 'binomial')
    
    pred.df <- as.data.frame(summary(mod)$coefficients[1,])
    pred.df$taxa.rand <- summary(mod)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(mod)$varcor$length[1,1]
    return(pred.df)
  })
  
  names(rand.preds) <- c("instant", "cumul", "novel")
  
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.43,0.65,0.51)
  circle.radius <- c(0.2, 0.22)
  circle.y <- 0.5
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  par(lheight=0.85)
  text(x=circle.cent[3]+0.02, y=circle.y+0.35, 
       adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=1.4)
  
  par(xpd=FALSE)
  # turn I circle to spatial polygons to get overlap
  I.sp <- Polygon(cbind(I$x, I$y))
  I.sp <- SpatialPolygons(list(Polygons(list(I.sp), ID = "a")))
  
  C.sp <- Polygon(cbind(C$x, C$y))
  C.sp <- SpatialPolygons(list(Polygons(list(C.sp), ID = "a")))
  
  N.sp <- raster::intersect(I.sp, C.sp)
  I.sub <- gDifference(I.sp, C.sp)
  C.sub <- gDifference(C.sp, I.sp)
  
  plot(I.sub, col="red", add=TRUE, border="black", lwd=1)
  plot(C.sub, col=cumul.col, add=TRUE, border="black", lwd=1)
  plot(N.sp, col="orange", add=TRUE, border="black", lwd=1.5)
  
  text(x=0.225, y=circle.y+0.2,
       labels=c("Instantaneous\nnovelty (I)    "), 
       col="red", adj=1, cex=1.4)
  text(x=0.82, y=circle.y+0.2,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=1.4)
  
  overall.means <- t(sapply(rand.preds, function(x){
    x <- t(x)
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.073, 0.081, 0.02)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.8)
  
  text(x = c(circle.cent + c(-0.086, 0.1, 0.015)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=1.7)
}
