# MDS PLOTTER

# used for plotting a three-way anosim. Dperecated

mds.plotter <- function(matrix, col.vec, ylim, xlim){
  
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
  plot(NMDS, type = "n", ylim=ylim, xlim = xlim,
       yaxt = "n", xaxt = "n")
  points(x = NMDS$points[,"MDS1"], 
         y = NMDS$points[, "MDS2"], 
         bg = col.vec,
         pch = 21,
         col = "black",
         cex = 1.4)
  
  for (i in 1:(nrow(NMDS$points)-1)){
    
    arrows(x0 = NMDS$points[i,"MDS1"], 
           y0 = NMDS$points[i, "MDS2"], 
           x1 = NMDS$points[i+1,"MDS1"], 
           y1 = NMDS$points[i+1, "MDS2"], length = 0.05, lwd = 1)
  }
  
}

pdf(file = "/Users/sassen/Desktop/ANOSIM.pdf",
    width = 8,
    height = 10)

# Set parameters
par(mfrow = c(3,3))
par(oma = c(4,4,2,2))

output <- novel.length.algo(matrices[234])

# find the first index of emergence
first <- output[[1]]$Emergence_index
total <- output[[1]]$Timeseries_length

# Plot some good examples
for(i in 1:3){
  
  if(i == 1){
    col.vec <- c(rep("grey", first-1), rep("orange", (total-first+1)))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 2){
    col.vec <- c(rep("grey", first-1), rep("orange", (total-first-5)), rep("grey", 6 ))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 3){
    col.vec <- c(rep("grey", first-1), rep("orange", 1), rep("grey", total-first))
    
    par(mar = c(2,0,0.1,0))
  }
  
  # Create the plot!
  mds.plotter(nov.matrices[[1]], 'orange', 
              ylim = c(-0.5, 2.5), xlim = c(-1, 1.5)) # Blip
  
  axis(side =1, line =0, at= c(-1, 0, 1))
  
  if(i == 1){
    axis(side =2, line =0)
    text(-1, 2.4, "a", font=4)
  }
  if(i == 3){
    axis(side = 4, line = 0)
  }
  
  # Extract R and Length
  set.length <- length(which(col.vec == "orange"))
  index <- which(output[[1]]$data$Set == set.length)
  set.R <- output[[1]]$data[index, "R"]
  text.1 <- paste0("R  = ", round(set.R,digits = 2))
  text.2 <- paste0("Length = ", set.length)
  text(1, 2.4, text.1, font =4, pos=1)
  text(1,2.15, text.2 , font =4,pos=1)
  
}

output <- novel.length.algo(matrices[220])
# find the first index of emergence
first <- output[[1]]$Emergence_index
total <- output[[1]]$Timeseries_length
steps <- output[[1]]$Length_steps

for(i in 1:3){
  
  if(i == 1){
    col.vec <- c(rep("grey", nrow(matrices[[220]])-steps), rep("orange",steps ))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 2){
    col.vec <- c(rep("grey", nrow(matrices[[220]])-steps), rep("orange",steps-1), rep("grey", 1))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 3){
    col.vec <- c(rep("grey", nrow(matrices[[220]])-steps), rep("orange",steps-2), rep("grey", 2))
    
    par(mar = c(2,0,0.1,0))
  }
  
  # Create the plot!
  mds.plotter(matrices[[220]], col.vec, 
              ylim = c(-0.5, 1), xlim = c(-.5, 1.5)) # Full persister
  
  axis(side =1, line =0, at= c(-1, 0, 1))
  
  if(i == 1){
    axis(side =2, line =0)
    text(-.45, 1.25, "b", font=4)
  }
  
  if(i == 3){
    axis(side = 4, line = 0)
  }
  
  # Extract R and Length
  set.length <- length(which(col.vec == "orange"))
  index <- which(output[[1]]$data$Set == set.length)
  set.R <- output[[1]]$data[index, "R"]
  text.1 <- paste0("R  = ", round(set.R,digits = 2))
  text.2 <- paste0("Length = ", set.length)
  text(1, 1.25, text.1, font =4)
  text(1,1, text.2 , font =4)
  
  
}

output <- novel.length.algo(matrices[[1]][["G8824"]])
# find the first index of emergence
first <- output[[1]]$Emergence_index
total <- output[[1]]$Timeseries_length
steps <- output[[1]]$Length_steps

for(i in 1:3){
  
  if(i == 1){
    col.vec <- c(rep("grey", nrow(matrices[[237]])-first), rep("orange",first))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 2){
    col.vec <- c(rep("grey", nrow(matrices[[237]])-first), rep("orange",3), rep("grey", total-first-3))
    
    par(mar = c(2,0,0.1,0))
  }
  if(i == 3){
    col.vec <- c(rep("grey", nrow(matrices[[237]])-first), rep("orange",1), rep("grey", total-first-1))
    
    par(mar = c(2,0,0.1,0))
  }
  
  mds.plotter(matrices[[237]], col.vec, 
              ylim = c(-.75, 1), xlim = c(-.75, 1.5)) # Short persister
  
  axis(side =1, line =0, at= c(-1, 0, 1))
  
  if(i == 1){
    axis(side =2, line =0)
    text(-.7, 1.25, "c", font=4)
  }
  
  if(i == 3){
    axis(side = 4, line = 0)
  }
  
  # Extract R and Length
  set.length <- length(which(col.vec == "orange"))
  index <- which(output[[1]]$data$Set == set.length)
  set.R <- output[[1]]$data[index, "R"]
  text.1 <- paste0("R  = ", round(set.R,digits = 2))
  text.2 <- paste0("Length = ", set.length)
  text(1, 1.25, text.1, font =4)
  text(1,1, text.2 , font =4)
}

mtext("MDS1", side=1, line=2, cex=1, col="black", outer=TRUE)
mtext("MDS2", side=2, line=2, cex=1, col="black", outer=TRUE)

dev.off()
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
