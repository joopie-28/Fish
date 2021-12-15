### Building matrices based on varying bin widths as opposed to just taking time as is ###

# We are going to test different bin widths: 2 and 4 years.

# Our initial analysis will be done for binary, I think it makes sense to use some code that creates a complete binary
# dataframe, before any cut-offs are applied. 

binary_aggregate_function <- function(survey_identifier, bin_width){
  
  # Here I create some reference lists 
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  year_unique <- as.data.frame(unique(sub_series$Year))
  min_year <- min(sub_series$Year)
  max_year <- max(sub_series$Year)
  
  # We are now condensing the data quite significantly, so I'm going to write off smaller surveys
  if (bin_width == 4){
    if (nrow(year_unique) < 24){
      return("Not enough data when bin size is 4 years")
    }
  }
  
  if (bin_width == 2){
    if (nrow(year_unique) < 12){
      return ("Not enough data when bin size is 2 years")
    }
  }
  
  # Here I create the empty dataframe which will be populated with presence/absence values, based on the bin width
  binary_df <- as.data.frame(matrix(data = NA, nrow = (ceiling(((max_year-min_year+1)/bin_width))), ncol = nrow(ref_list)))
  colnames(binary_df) <- ref_list[, 1]
  
  # Select the correct years as a function of bin width
  if (bin_width != 1){
    name_list <- list()
    name_list <- seq(min_year, max_year, by = bin_width)
    name_list <- name_list[order(as.numeric((name_list)), decreasing = FALSE)]
    rownames(binary_df) <- name_list
  }
  
  # Let's populate this data frame. I will use the column and row names to loop through a subset of the data
  # a couple times. This will essentially crosscheck if a species is present in the set or not. 
  for (i in 1:nrow(binary_df)) {
    sub_temp <- subset(sub_series, Year >= as.numeric(rownames(binary_df)[i]) & Year < (as.numeric(row.names(binary_df)[i]) + (bin_width)))
    for (j in 1:ncol(binary_df)) {
      if (colnames(binary_df[j]) %in% sub_temp$Species){
        binary_df[i,j] <- 1
      }
      else{
        binary_df[i,j] <- 0
      }
    }
  }

  # Order dataframe so that novelty framework does not get confused
  binary_df <- binary_df[order(as.numeric((row.names(binary_df))), decreasing = FALSE),]
  
  # Very small samples will be returned as a string of numbers (type = double) rather than data frames, so will instruct the function to skip this survey if
  # That happens. The timeseries is just too short basically.
  if (typeof(binary_df) == "double"){
    return("Not enough data")
  }
  
  # Change row names to bin (time)
  row.names(binary_df) <- (2021-as.numeric(row.names(binary_df)))
  
  # Now remove rows that have sums of 0, these are rows of years that were actually not present in the original data.
  binary_df <- binary_df[rowSums(binary_df[])> 0, ]
  
  # Make sure the order is OK, it will influence the model results (our priamry reference community HAS to be the oldest one for the
  # model to make sense).
  
  return(binary_df)
}

test <- binary_aggregate_function("G524", 2)

list_matrix_B_bins_function <- function(check_list, bin_width){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[[i]][], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- binary_aggregate_function(TimeSeries_ID, bin_width)
                    if(class(temp)=="matrix"){
                      # remove bins within sites with fewer species than cut off
                      temp <- temp[rowSums(temp) > rich.cutoff,]
                    }
                    return(temp)
                  }) 
    
    if (names(check_list[i]) == "palearctic_ID"){
      palearctic_mat_B <- nam
      names(palearctic_mat_B) <- check_list$palearctic_ID
    }
    if (names(check_list[i]) == "nearctic_ID"){
      nearctic_mat_B <- nam
      names(nearctic_mat_B) <- check_list$nearctic_ID
    }
    if (names(check_list[i]) == "afrotropics_ID"){
      afrotropics_mat_B <- nam
      names(afrotropics_mat_B) <- check_list$afrotropics_ID
    }
    if (names(check_list[i]) == "neotropics_ID"){
      neotropics_mat_B <- nam
      names(neotropics_mat_B) <- check_list$neotropics_ID
    }
    if (names(check_list[i]) == "australasia_ID"){
      australasia_mat_B <- nam
      names(australasia_mat_B) <- check_list$australasia_ID
    }
    
  }
  list_matrix <- list(palearctic_mat_B, nearctic_mat_B, afrotropics_mat_B, neotropics_mat_B, australasia_mat_B)
  names(list_matrix) <- c("palearctic_mat_B","nearctic_mat_B","afrotropics_mat_B","neotropics_mat_B", "australasia_mat_B")
  return(list_matrix)
}

matrix_list_B_bin2 <- list_matrix_B_bins_function(ID_list, 3)

# Calculate novelty and return output in a list
list_novelty_function <- function(matrix_list){
  
  for (i in 1:length(matrix_list)) {
    
    nam <- lapply(names(matrix_list[[i]][]), 
                  function(ID){
                    print(ID)
                    site.sp.mat <- (matrix_list[[i]][ID])
                    site.sp.mat <- site.sp.mat[[ID]]
                    #site.sp.mat <- site.sp.mat[[ID]]
                    # This line had to be added because there was one timeseries with 0 change over 20 years....
                    if (ID == "G7555"){
                      return(NA)
                    }
                    else{
                      if(typeof(site.sp.mat) == "character"){
                        return(NA)
                      }
                      else{
                        if (nrow(site.sp.mat) >= 10 & ncol(site.sp.mat) >=5) {
                          
                          temp <- identify.novel.gam(site.sp.mat = site.sp.mat, alpha = 0.05, metric = "jaccard", site = ID, plot = TRUE, plot.data = FALSE,
                                                     gam.max.k = -1)
                          # Remove first 5 bins
                          temp <- temp[-c(1:5),]
                          return(temp)
                        }
                        else {
                          return(NA)
                        }
                      }
                    }  
                  })
    
    # Select correct data, delete NA's
    
    if (names(matrix_list[i]) == "palearctic_mat_B"){
      palearctic_novelty_B <- nam
      names(palearctic_novelty_B) <- ID_list$palearctic_ID
      palearctic_novelty_B <- palearctic_novelty_B[!sapply(palearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "nearctic_mat_B"){
      nearctic_novelty_B <- nam
      names(nearctic_novelty_B) <- ID_list$nearctic_ID
      nearctic_novelty_B <- nearctic_novelty_B[!sapply(nearctic_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "afrotropics_mat_B"){
      afrotropics_novelty_B <- nam
      names(afrotropics_novelty_B) <- ID_list$afrotropics_ID
      afrotropics_novelty_B <- afrotropics_novelty_B[!sapply(afrotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "neotropics_mat_B"){
      neotropics_novelty_B <- nam
      names(neotropics_novelty_B) <- ID_list$neotropics_ID
      neotropics_novelty_B <- neotropics_novelty_B[!sapply(neotropics_novelty_B, function(x) all(is.na(x)))]
    }
    if (names(matrix_list[i]) == "australasia_mat_B"){
      australasia_novelty_B <- nam
      names(australasia_novelty_B) <- ID_list$australasia_ID
      australasia_novelty_B <- australasia_novelty_B[!sapply(australasia_novelty_B, function(x) all(is.na(x)))]
    }
    
    # Add a list of TimeSeries which will come in handy when evaluating the results.
    
    
  }
  list_novelty <- list(palearctic_novelty_B, nearctic_novelty_B, afrotropics_novelty_B, neotropics_novelty_B, australasia_novelty_B)
  names(list_novelty) <- c("palearctic_novelty_B","nearctic_novelty_B","afrotropics_novelty_B","neotropics_novelty_B", "australasia_novelty_B")
  return(list_novelty)
  
}

novelty_list_B_bin3 <- list_novelty_function(matrix_list_B_bin2)

novelty_analysis_output_B_bin3 <- novel.probability(novelty_list_B_bin4)


venn_plot_function <- function(prob.model.list){
  
  model.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$fixed.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
  rand.preds <- lapply(c(4:6), function(n){
    temp.pred <- prob.model.list$random.prob.models[[n]]$pred.df
    return(temp.pred)
  })
  
  require(plotrix)
  require(sp)
  require(raster)
  require(rgeos)
  require(shape)
  
  cumul.col <- rgb(0.373,0.651,0.765)
  
  if(!is.null(dev.list())){
    dev.off(which=dev.list())
    close.screen(all.screens=TRUE)
  }
  
  
  
  
  split.screen(rbind(c(0.01,0.99,0.05,0.99), # vens
                     
                     c(0.18,0.44,0.08,0.4), # smallplots
                     c(0.7,0.96,0.08,0.4), # smallplots
                     c(0.44,0.7,0.08,0.4))) # smallplots
  
  screen(1)
  par(mar=c(0,0,0,0), ps=8, tcl = -0.2, mgp=c(3,0.5,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
       xlab="", ylab="", xaxs="i")
  
  # top circle
  circle.cent<-c(0.43,0.65,0.54)
  circle.radius <- c(0.22, 0.2)
  circle.y <- 0.65
  
  I<-draw.circle(y=circle.y, x=circle.cent[1], radius=circle.radius[1], nv=500, 
                 border="black", col="white", lwd=2)
  C<-draw.circle(y=circle.y, x=circle.cent[2], radius=circle.radius[2], nv=500, 
                 border="black", col="white", lwd=2)
  
  par(xpd=NA)
  segments(x0=circle.cent + c(0, 0, 0.0175), 
           y0=rep(circle.y, 3) + c(0, 0, 0.24),
           x1=circle.cent + c(-0.15,0.2, 0.0175),
           y1=c(0.3, 0.3, 0.3),
           lwd=1)
  
  rect(xleft=0.4, xright=0.6, ybottom=circle.y+0.285,
       ytop=circle.y+0.45, border=NA, col="white")
  
  par(lheight=0.85)
  text(x=circle.cent[3]+0.02, y=circle.y+0.35, 
       adj=0.5, labels="Novel\ncommunity (N)", col="orange",
       cex=0.9)
  
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
       col="red", adj=1, cex=0.9)
  text(x=0.85, y=circle.y+0.2,
       labels=c("Cumulative\n   novelty (C)"), 
       col=cumul.col, adj=0, cex=0.9)
  
  overall.means <- t(sapply(rand.preds, function(x){
    
    fit <- plogis(x[1,"Estimate"])
    upper <- plogis(x[1,"Estimate"] + 1.96 * x[1, "Std. Error"])
    lower <- plogis(x[1,"Estimate"] - 1.96 * x[1, "Std. Error"])
    
    round(c(fit, lower, upper)*100, 1)
  }))
  
  text(x=c(circle.cent + c(-0.075, 0.085, 0.01)),
       y= rep(circle.y, 3) + 0.025,
       labels=paste0(sprintf(overall.means[,1], fmt = '%#.1f'), "%"), 
       font=2, cex=0.9)
  
  text(x = c(circle.cent + c(-0.085, 0.1, 0.01)),
       y= rep(circle.y - 0.025, 3),
       labels=paste0("(", sprintf(overall.means[,2], fmt = '%#.1f'), "-",
                     sprintf(overall.means[,3], fmt = '%#.1f'),"%)"),
       cex=0.8)
  
  text(x=relative.axis.point(0.1, "x"),
       y=relative.axis.point(0.95, "y"),
       labels="(A)", font=2, cex=0.9)
  
  close.screen(1)
  
  sapply(2:4, function(n){
    screen(n)
    
    par(mar=c(0,0,0,0), ps=8, tcl = -0.2, mgp=c(3,0.4,0), las=1)
    
    point.col <- c("red", cumul.col, "orange")[n-1]
    
    pred.df <- as.data.frame(model.preds[[n-1]])
    
    ylims <- c(0, 0.2)
    
    plot(x=NULL, y=NULL, xlab="", ylab="",
         xlim=c(0,5.25), xaxs="i",
         ylim=ylims, axes=FALSE)
    
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=par("usr")[4], 
         col="white", border=NA)
    
    if(n == 2){
      axis(side=2, at=seq(0,0.2,0.05), cex.axis=0.8)
      mtext(side=2, line=1.3, text="Probability", cex=0.8, las=0)
      axis(side=4, at=seq(0,0.2,0.05), labels=NA)
    } else {
      axis(side=2, at=seq(0,0.2,0.05), labels=NA)
    }
    
    if(n == 4){
      axis(side=4, at=seq(0,0.2,0.05), labels=NA, tcl= 0.2)
    }
    
    axis(side=1, at = c(0.5:5), labels=c("AFRO", "AUS", "NEA", "NEO", "PAL"),
         mgp=c(3,-0.1,0), cex.axis=0.8)
    
    pred.df$upper <- plogis(pred.df[,"fit"] + 1.96 * pred.df[, "se.fit"])
    pred.df$lower <- plogis(pred.df[,"fit"] - 1.96 * pred.df[, "se.fit"])
    
    raw.data <- prob.model.list$fixed.prob.models[[n]]$model$data
    raw.data$success <- raw.data[, c("instant", "cumul", "novel")[(n-1)]]
    raw.data$failure <- raw.data[, c("non.instant", "non.cumul", "non.novel")[(n-1)]]
    raw.data$prop <- raw.data$success / (raw.data$success + raw.data$failure)
    
    points(x =jitter((as.numeric(raw.data$taxa)-0.5), amount=0.2),
           y = raw.data$prop, pch = 16, col=point.col, cex=0.4)
    
    segments(x0=0.5:4.5, x1=0.5:4.5,
             y0=pred.df$lower,
             y1=pred.df$upper)
    
    points(y=plogis(pred.df$fit), x = 0.5:4.5, bg="white",
           pch=21, cex=0.75, lwd=1.25)
    
    text(x=relative.axis.point(0.15, "x"),
         y=relative.axis.point(0.925, "y"),
         labels=paste0("(",LETTERS[c(2,4,3)][n-1],")"), font=2, cex=0.9)
    
    box()
    close.screen(n)
  })
}

venn_plot_function(novelty_analysis_output_B_bin2)
