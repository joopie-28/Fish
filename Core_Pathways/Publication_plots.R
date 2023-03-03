# Publication plots

# 22-11-2022 
# Script for publication-worthy plots


# All relevant functions
map.realm.plotter <- function(realm){
  # Import geodata from rnaturalearth package and set cropping parameters
  if(realm == 'Australasia'){
    biorealm_background <- ne_countries(scale = "Large", geounit = c('australia', 'new zealand', 'papua new guinea', 'indonesia'), type = 'map_units',
                                        returnclass = c("sp", "sf"))
    ymin = -47 
    ymax = 0 
    xmin = 110 
    xmax = 180
    
    text.offset = c(-2.5, 3, 7.5)
    fig.offset = c('y.offset' = 1.5, 'cex1' = 1, 'cex2'=0.8, 'circl.cent1' = 0.74, 'diam.large'=0.18, 'diam.small'=0.17)
  }
  if(realm == "Nearctic"){
    biorealm_background <- ne_countries(scale = "Large", geounit = c('United States of America', 'Canada', "Mexico"), type = 'sovereignty',
                                        returnclass = c("sp", "sf"))
    ymin = 25
    ymax = 70
    xmin= -130
    xmax= -55
    
    text.offset = c(-1.4, 2.1,7.6)
    fig.offset = c('y.offset' = 1, 'cex1' = 0.9, 'cex2'=0.7, 'circl.cent1' = 0.75,'diam.large'=0.15, 'diam.small'=0.145)
  }
  if(realm == "Palearctic"){
    biorealm_background <- ne_countries(scale = "Large", continent = 'Europe', type = 'sovereignty',
                                        returnclass = c("sp", "sf"))
    ymin = 40
    ymax = 70
    xmin= -12
    xmax= 50
    
    text.offset = c(-1.95, 1.8,7)
    fig.offset = c('y.offset' = 0.5, 'cex1' = 1, 'cex2'=0.8, 'circl.cent1' = 0.77, 'diam.large'=0.18, 'diam.small'=0.17)
    
  }
  
  # Set up the plotting background
  par(mar=c(2,2,2,2))
  realm_bbox = st_bbox(c(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax))
  
  plot(crop(biorealm_background, realm_bbox), col = 'grey97', axes = F, bg = 'white', xaxs = 'i', yaxs = 'i' )
  if(realm == "Palearctic"){
    mtext("b)", side = 3, line=0.5, adj = c(-4,0))
  }
  if(realm == "Nearctic"){
    mtext("d)", side = 3, line=0.5, adj = c(-4,0))
  }
  if(realm == "Australasia"){
    mtext("c)", side = 3, line=0.5, adj = c(-4,0))
  }
  
  
  if(realm=='Australasia'){
    segments(x0=149, y0=-29, x1 =149 , y1 =-25 )
    segments(x0=149, y0=-29, x1 =155 , y1 =-29 )
    segments(x0=155, y0=-29, x1 =155 , y1 =-25 )
    segments(x0=155, y0=-25, x1 =149 , y1 =-25 )
    segments(x0 = 155, y0 = -29, x1=135.2, y1= -47.6   )
    segments(x0 = 149, y0 = -25, x1=,116, y1= -29   )
  }
  col.df <- subset(full.novel.mat, cat =='novel')
  points(x=col.df$long, y=col.df$lat, pch = 19, cex=0.75, col = 'orange')
  non.df <- subset(full.novel.mat, cat !='novel')
  points(x=non.df$long, y=non.df$lat, pch = 19, cex=0.75, col = 'black')
  points(x=col.df$long, y=col.df$lat, pch = 19, cex=0.75, col = 'orange')
  box(which = 'plot',  lty = 1)
  
  # Plot Japan inset for the Palearctic zone
  if(realm == "Palearctic"){
    par(fig = c(0.55, .926, 0.035, 0.5), new = T) 
    jap_bbox = st_bbox(c(ymin = 32, ymax = 50, xmin =130, xmax = 150))
    plot(crop(ne_countries(scale = 'Large', country = c('Japan'), type = 'sovereignty',
                           returnclass = c("sp", "sf")), jap_bbox), col = 'grey97', bg= 'white')
    col.df <- subset(full.novel.mat, cat =='novel')
    non.df <- subset(full.novel.mat, cat !='novel')
    points(x=non.df$long, y=non.df$lat, pch = 19, cex=0.75, col = 'black')
    points(x=col.df$long, y=col.df$lat, pch = 19, cex=0.75, col = 'orange')
    box(which = 'plot',  lty = 1)
  }
  # Plot QLD inset for Aus zone
  if(realm == "Australasia"){
    par(fig = c(0.1, .4, 0.035, 0.4), new = T, mai= c(0,0,0,0)) 
    qld_bbox = st_bbox(c(ymin = -29, ymax = -25, xmin =149, xmax = 160))
    plot(crop(ne_countries(scale = 'Large', country = c('Australia'), type = 'sovereignty',
                           returnclass = c("sp", "sf")), qld_bbox), col = 'grey97', bg= 'white')
    col.df <- subset(full.novel.mat, cat =='novel')
    non.df <- subset(full.novel.mat, cat !='novel')
    points(x=non.df$long, y=non.df$lat, pch = 19, cex=0.75, col = 'black')
    points(x=col.df$long, y=col.df$lat, pch = 19, cex=0.75, col = 'orange')
    box(which = 'plot',  lty = 1)
    
    
  }
  
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
  
  text(x=c(circle.cent + c(-0.073, 0.081, 0.04)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.8)
  
  text(x = c(circle.cent + c(-0.086, 0.1, 0.012)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=1.7)
 
  
  
  
  
  
}

emergence_modeller <- function(realm){
  
  # Australia has a very good record in terms of sampling gaps.
  # The bin-lag variabel thus can not be modelled for this region.
  taxa.prob.models <- lapply(1:3, function(n){
    
    success_cat = c("instant", "cumul", "novel")[n]
    
    full.novel.mat$success = full.novel.mat[, ..success_cat]
    
    if (realm == "Australasia"){
      mod<- glm(success~1, data = subset(full.novel.mat, BioRealm == realm), family = 'binomial')
    }
    else{
      mod<- glm(success~1+bin_lag+position, data = subset(full.novel.mat, BioRealm == realm), family = 'binomial')
    }
    
    
    
    pred.df <- as.data.frame(summary(mod)$coefficients)
    pred.df$taxa.rand <- summary(mod)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(mod)$varcor$length[1,1]
    return(pred.df)
  })
  
  names(taxa.prob.models) <- c("instant", "cumul", "novel")
  return(taxa.prob.models)
}

venn_plot_function.v3 <- function(realm){
  rand.preds <- lapply(1:3, function(n){
    success_cat = c("instant", "cumul", "novel")[n]
    
    full.novel.mat$success = full.novel.mat[, ..success_cat]
    if(realm != "Australasia"){
      mod<- glm(success~1+bin_lag+position, data = subset(full.novel.mat, BioRealm == realm), family = 'binomial')
    }else{
      mod<- glm(success~1, data = subset(full.novel.mat, BioRealm == realm), family = 'binomial')
      
    }
    pred.df <- as.data.frame(summary(mod)$coefficients[1,])
    pred.df$taxa.rand <- summary(mod)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(mod)$varcor$length[1,1]
    return(pred.df)
  })
  
  names(rand.preds) <- c("instant", "cumul", "novel")
  
  overall.means <- t(sapply(rand.preds, function(x){
    x <- t(x)
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.35,0.65,0.54)
  if(overall.means[1,1] > overall.means[2,1]){
    circle.radius <- c(0.3, 0.27)
  }else{
    circle.radius <- c(0.27, 0.3)
  }
  
  circle.y <- 0.5
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  
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
  
  
  text(x=c(circle.cent + c(-0.075, 0.085, -0.04)),
       y= rep(circle.y, 3) + 0.05,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.3)
  
  text(x = c(circle.cent + c(-0.085, 0.1, -0.04)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=1.1)
}

venn_plot_function.v4 <- function(mod.list){
  rand.preds <- lapply(1:3, function(n){
    success_cat = c("instant", "cumul", "novel")[n]
    
    full.novel.mat$success = full.novel.mat[, success_cat]
      mod<- glm(success~bin_lag+position + BioRealm, data = subset(full.novel.mat, BioRealm == realm), family = 'binomial')

     
    
    pred.df <- as.data.frame(summary(mod)$coefficients[1,])
    pred.df$taxa.rand <- summary(mod)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(mod)$varcor$length[1,1]
    return(pred.df)
  })
  
  names(rand.preds) <- c("instant", "cumul", "novel")
  
  overall.means <- t(sapply(rand.preds, function(x){
    x <- t(x)
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.35,0.65,0.54)
  if(overall.means[1,1] > overall.means[2,1]){
    circle.radius <- c(0.3, 0.27)
  }else{
    circle.radius <- c(0.27, 0.3)
  }
  
  circle.y <- 0.5
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  
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
  
  
  text(x=c(circle.cent + c(-0.075, 0.085, -0.04)),
       y= rep(circle.y, 3) + 0.05,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=1.3)
  
  text(x = c(circle.cent + c(-0.085, 0.1, -0.04)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=1.1)
}

### All plots created from here on ####

# Plot 1. SIMPROF-NMDS example
mds.cluster.plotter <- function(matrix){
  
  # Set plotting params
  par(mfrow = c(1,2))
  
  
  
  
  # Extract labels
  label <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                  alpha = 0.05,
                                  metric = "bray",
                                  plot = F, 
                                  site = "NA",
                                  plot.data = FALSE,
                                  gam.max.k = -1)
  
  matrix$category <- label$cat
  matrix$colour <- NA
  
  # Assign colours
  for (i in 1:nrow(matrix)) {
    if(matrix$category[i] == "cumul"){
      matrix$colour[i] <- "skyblue"
    }
    if(matrix$category[i] == "instant"){
      matrix$colour[i] <- "red1"
    }
    if(matrix$category[i] == "novel"){
      matrix$colour[i] <- "orange"
    }
    if(matrix$category[i] == "back"){
      (matrix$colour[i] <- "grey")}
  }
  
  # Run MDS 
  NMDS=metaMDS(matrix[,-c(ncol(matrix), (ncol(matrix)-1))], 
               k=2, trymax = 10000)
  plot(NMDS, type = "n", ylab = "MDS2", xlab = "MDS1")
  usr<-par("usr")
  text(x = median(c(usr[3], usr[1])), 
       y = usr[4],
       labels = "a)",
       adj = c(-0.6, 1.6),
       col = "black")
  
  points(x = NMDS$points[,"MDS1"], 
         y = NMDS$points[, "MDS2"], 
         bg = matrix$colour,
         pch = 21,
         col = "black",
         cex = 1.4)
  
  for (i in 1:(nrow(NMDS$points)-1)){
    
    arrows(x0 = NMDS$points[i,"MDS1"], 
           y0 = NMDS$points[i, "MDS2"], 
           x1 = NMDS$points[i+1,"MDS1"], 
           y1 = NMDS$points[i+1, "MDS2"], length = 0.05, lwd = 1)
  }
  
  print("Clustering")
  
  test <- simprof(data = matrix[,-c(ncol(matrix), (ncol(matrix)-1))], num.expected = 1000,
                  num.simulated = 999, method.distance ="czekanowski", 
                  method.cluster = "average"
                  ,alpha=0.05, undef.zero = T)
  
  temp <- simprof.plot(test, plot = F)
  
  
  dendro.col.df <- data.frame(labels = as.numeric(labels(temp)))
  dendro.col.df$colour <- NA
  dendro.col.df$cat <- NA
  
  for(i in 1:nrow(dendro.col.df)){
    
    index <- which(as.numeric(label$bins) == dendro.col.df$labels[i])
    dendro.col.df$cat[i] <- label$cat[index]
    
    if(dendro.col.df$cat[i] == "cumul"){
      dendro.col.df$colour[i] <- "skyblue"
    }
    if(dendro.col.df$cat[i] == "instant"){
      dendro.col.df$colour[i] <- "red1"
    }
    if(dendro.col.df$cat[i] == "novel"){
      dendro.col.df$colour[i] <- "orange"
    }
    if(dendro.col.df$cat[i] == "back"){
      dendro.col.df$colour[i] <- "grey"}
  }
  
  # Plot the Dendrogram
  
  labels_colors(temp) <- dendro.col.df$colour
  
  labels(temp) <- paste0(dendro.col.df$cat, "-", dendro.col.df$labels)
  
  plot(temp, ylab = "Height")
  usr<-par("usr")
  text(x = median(c(usr[3], usr[1])), 
       y = usr[4],
       labels = "b)",
       adj = c(-1.8, 1.6),
       col = "black")
  
}

pdf(file = "/Users/sassen/Desktop/plot_1.pdf",
    width = 12,
    height = 7)

mds.cluster.plotter(matrices[[1]][["G7620"]])
dev.off()

#### Plot 2A. Global emergence rates Venn ###
pdf(file = "/Users/sassen/Desktop/plot_2a.pdf",
    width = 15,
    height = 12)

venn_plot_main(full.novel.mat.season)

dev.off()

### Plot 2b. Palearctic emergence rates ####

pdf(file = "/Users/sassen/Desktop/plot_2b.pdf",
    width = 15,
    height = 12)

map.realm.plotter('Palearctic')
par(fig = c(0.55, 1, 0.535, 1), new = T) 
venn_plot_function.v3('Palearctic')
dev.off()


### Plot 2c. Nearctic emergence rates ####
pdf(file = "/Users/sassen/Desktop/plot_2c.pdf",
    width = 15,
    height = 12)

map.realm.plotter('Nearctic')
par(fig = c(0.55, 1, 0.535, 1), new = T) 
venn_plot_function.v3('Nearctic')
dev.off()

### Plot 2d. Australasian emergence rates ####
pdf(file = "/Users/sassen/Desktop/plot_2d.pdf",
    width = 15,
    height = 12)

map.realm.plotter('Australasia')
par(fig = c(0.55, 1, 0.535, 1), new = T) 
venn_plot_function.v3('Australasia')
dev.off()


# Plot 3. Trajectory summary stats

figure.3 <- function(full.novel.mat){
  nov.sub<-subset(full.novel.mat.season, novel.class != "NONE"& novel.class != "END")
  nov.sub$novel.class <- as.factor(nov.sub$novel.class)
  nov.tbl <- table(nov.sub$novel.class)/nrow(nov.sub)
  names(nov.tbl) <- c("Blip Event", 'Persistent State')
  barplot(nov.tbl, ylim = c(0,0.75), ylab = 'Percentage of Novel Communities', col ="skyblue4")
  box()
}

pdf(file = "/Users/sassen/Desktop/Novel_communities/Main_text/figures/plot_3.pdf",
    width = 12,
    height = 12)

par(mar=c(5.1,7.5,4.1,7.5))
figure.3(full.novel.mat.season)
dev.off()

# Plot 4.a-b Extinction and origination rates vioplot

figure.4b <- function(data.nov){
  layout(matrix(1:2, nrow=1))
  # Relevel variables for more consistent plot
  data.nov$novel.class <- factor(data.nov$novel.class , levels=c("BLIP",'NONE', "Persister"))
  par(mar=c(4,4,2,0.5))
  vioplot(orig~novel.class, data.nov, col = "green", horizontal = F, axes=F, yaxt='n', xlab="", ylab='')
  axis(side = 2)
  axis(side = 1, at =1:3 , labels=c("Background",'Blip', "Persister"))
  box()
  mtext("Local Origination (Number of Species)", side = 2, line=2)
  mtext("a)", side = 3, line=0.5, adj = c(-4,0))
  
  par(mar=c(4,0.5,2,4))
  vioplot(ext~novel.class, data.nov, col = "red", horizontal = F,yaxt = "n", axes=F, xlab="", ylab='')
  axis(side = 4)
  axis(side = 1, at =1:3 , labels=c("Background",'Blip', "Persister"))
  box()
  mtext("Local Extinction (Number of Species)", side = 4, line=2)
  mtext("b)", side = 3, line=0.5, adj = c(-4,0))
}

pdf(file = "/Users/sassen/Desktop/plot_4ab.pdf",
    width = 22,
    height = 12)
data.nov=subset(full.novel.mat,novel.class != 'END')
figure.4b(data.nov)

dev.off()

# Plot 4c-d. Persistence and demography

pdf(file = "/Users/sassen/Desktop/plot_4cd.pdf",
    width = 22,
    height = 12)

data.nov=subset(full.novel.mat,novel.class != 'END' & novel.class != 'NONE')
demo.persistent.plots(persistent.trajectory.mod)
dev.off()

######################################################
######### Tables. ####################################
######################################################

csv.writer <- function(model, name){
  
  write.csv(as.data.frame(summary(model)$coef), paste0("./outputs/", 
                                                     name, ".csv"))
  print(paste0('CSV ', name, ' written to outputs'))
}

## Table S1 A, B, C,. Global Emergence Rates ####
# Already created in other functions

## Table S2. A B C, Realm-specific rates ##

pal<-emergence_modeller('Palearctic')$novel
nea<-emergence_modeller('Nearctic')$novel
aus<-emergence_modeller('Australasia')$novel
realm.list <- list(pal, nea, aus)
names(realm.list) <- c('pal', 'nea', 'aus')

write.csv(pal, "./outputs/pal-novelty.csv")
write.csv(nea, "./outputs/nea-novelty.csv")
write.csv(aus, "./outputs/aus-novelty.csv")

## Table S3. Timeseries level - emergence rates ##

csv.writer(timeseries.level.rates, 'timeseries_level_rates')

## Table S4. Local demography effects on persistence ##

csv.writer(persistent.trajectory.mod, 'demography_persistence')

# Table S5. Origination and Extinction rates

data.nov=subset(full.novel.mat,novel.class != 'END'  & cat != 'cumul' & cat != 'instant')
ext.mod<-(glmer(ext~novel.class + (1|site), data=data.nov, family = "poisson"))
orig.mod<-(glmer(orig~novel.class + (1|site), data=data.nov, family = "poisson"))

csv.writer(ext.mod, 'extinction_byclass')
csv.writer(orig.mod, 'origination_byclass')




### Supplementary figures and tables for Simon ###

# 1. model persistence type as a function of lag between nc and next bin #


pdf(file = "/Users/sassen/Desktop/Novel_communities/Comments_and_Supps/figures/plot_1.pdf",
    width = 22,
    height = 12)

visreg(lag.mod, 'lag_to_next', scale = "response", partial = F, rug = F, ylim = c(0,1),
       ylab = 'persistence type (0 = blip, 1 = persistent)', xlab = 'lag between novel and succeeding community (years)')
points(novel.class.binary~lag_to_next, data = lag.mod.data, pch=19, col = alpha('black', alpha=0.4), cex=1)

dev.off()

pdf(file = "/Users/sassen/Desktop/Novel_communities/Comments_and_Supps/figures/plot_2.pdf",
    width = 22,
    height = 12)

hist(lag.mod.data$lag_to_next, xlim = c(0,15), xlab  = 'lag between novel and succeeding bin', main = '')

dev.off()
csv.writer(lag.mod, 'class_lag_effects')




# 2. Distribution of novel communities persistence for varying Tau ####

# A line graph with Tau on the x-axis, and two y-axes: left shows number of blips, 
# right shows number of persistent.

pdf(file = "/Users/sassen/Desktop/Novel_communities/Comments_and_Supps/figures/plot_3.pdf",
    width = 18,
    height = 12)


decay.df<- data.frame('Persister' = c(
length(which(full.novel.mat$novel.class == "Persister")),
length(which(full.novel.mat$class.T3 == "Persister")),
length(which(full.novel.mat$class.T4 == "Persister")),
length(which(full.novel.mat$class.T5 == "Persister")))/487,

"Blip" = c(

length(which(full.novel.mat$novel.class == "BLIP")),
length(which(full.novel.mat$class.T3 == "BLIP" |full.novel.mat$class.T3 == "Blip")),
length(which(full.novel.mat$class.T4 == "BLIP"|full.novel.mat$class.T4 == "Blip")),
length(which(full.novel.mat$class.T5 == "BLIP"|full.novel.mat$class.T5 == "Blip")))/487)

decay.df$Tau <- c(2,3,4,5)
plot(decay.df$Persister~decay.df$Tau, type = 'l', col = 'red', lwd = 1.5,
     xlab = 'T', ylab = 'Proportion', ylim = c(0,1))
lines(decay.df$Blip~decay.df$Tau, type = 'l', col = 'blue', lwd = 1.5)
legend('topleft', legend=c("Blip", "Persistent"),
       col=c("blue", "red"), lty=1, cex=0.8)

dev.off()
