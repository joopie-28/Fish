##### This is the script that includes the pre-processing steps #####

#### Characterizing the drivers of novel freshwater fish communities ####

novelty_progressor_function <- function(matrix,
                                        metric, 
                                        site){
  
  
  if (typeof(matrix) == "character"){
    return(NA)
  }
  
  if (is.na(matrix)){
    return(NA)
  }
  # Compute contribution of natives
  matrix_invasive <- matrix
  matrix_invasive$NNC <- 0
  matrix_invasive$NNC_increase <- 0
  matrix_invasive$bins <- as.numeric(rownames(matrix_invasive))
  
  for (i in 1:nrow(matrix_invasive)) {
    NNC <- 0
    for (j in 1:ncol(matrix_invasive)) {
      if(stri_detect_fixed(colnames(matrix_invasive)[j], "xotic")){
        NNC <- NNC + matrix_invasive[i,j]
      }
      
    }
    matrix_invasive$NNC[i] <- NNC/rowSums(matrix_invasive[i,c(1:(ncol(matrix_invasive)-3))])
  }
  
  for (i in 2:nrow(matrix_invasive)) {
    
    matrix_invasive$NNC_increase[i] <- (matrix_invasive$NNC[i] - matrix_invasive$NNC[i-1])
  } 
  
  
  
  # Need to make a frame with bins and distance
  
  # Obtain distance matrix for all sites in timeseries
  site.dist <- as.data.frame(as.matrix(vegdist(matrix,
                                               method=metric)))
  
  # Lets label!
  
  label_frame <- identify.novel.gam(matrix, 
                                    alpha = 0.05,
                                    metric = metric,
                                    plot = TRUE, 
                                    site = site,
                                    plot.data = FALSE,
                                    gam.max.k = -1)
  
  # Going to remove the first 5 of the label frame because
  # these communities are unreliable. The effect of this is
  # then that we do not consider novelty in these first 5.
  # Assign the first 5 a "background" status, solves for the
  # wrongful recognition of novelty in the first 5.
  
  for (i in 1:5) {
    rownames(site.dist)[i] <- paste0("back-", rownames(site.dist)[i])
    colnames(site.dist)[i] <- paste0("back-", colnames(site.dist)[i])
    
  }
  
  # Quick loop to assign names
  
  for (i in 6:dim(site.dist)[1]) {
    for (j in 6:dim(label_frame)[1]) {
      
      if ((rownames(site.dist)[i]) == (label_frame$bins[j])){
        rownames(site.dist)[i] <- paste0(label_frame$cat[j], "-", rownames(site.dist)[i])
      }
      
      if((colnames(site.dist)[i]) == (label_frame$bins[j])){
        colnames(site.dist)[i] <- paste0(label_frame$cat[j], "-", colnames(site.dist)[i])
      }
    } 
  } 
  
  
  # Just need to control for novel communities at the very end of a timeseries,
  # as there is no following state.
  
  
  
  pre_novelty <- (site.dist) %>% 
    dplyr::select(starts_with("novel")) # This selects the community directly before
  
  
  # Control for no-novelty scenario's
  if(ncol(pre_novelty) == 0){
    return(NA)
  }
  
  index <- which(rownames(pre_novelty) %in% colnames(pre_novelty[1]))
  
  if(nrow(pre_novelty) == index){
    return(NA)
  }
  
  
  novelty_progression <- as.data.frame(pre_novelty[c(index:nrow(pre_novelty)), 1])
  colnames(novelty_progression) <- "Distance"
  
  bin_list <- list()
  
  for (i in 1:nrow(pre_novelty)) {
    a <- str_split(rownames(pre_novelty), pattern = "-")[[i]][2]
    bin_list <- append(bin_list, a)
  }
  
  novelty_progression$bins <- unlist(bin_list)[c(index:length(bin_list))]
  
  novelty_progression$time_after_novelty <- 0
  
  for (i in 2:nrow(novelty_progression)) {
    
    novelty_progression$time_after_novelty[i] <- as.numeric(novelty_progression$bins[i]) - as.numeric(novelty_progression$bins[1])
    
  }
  
  #novelty_progression$Distance <- (novelty_progression$Distance/max(novelty_progression$Distance)*100) # scale distance
  novelty_progression$site <- site
  
  
  return(novelty_progression)
  
}

output <- lapply(big_list[], 
                 function(TimeSeries_ID){
                   print((TimeSeries_ID))
                   matrix <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[[TimeSeries_ID]]
                   temp <- novelty_progressor_function(matrix = matrix,
                                                       metric = "bray",
                                                       site = TimeSeries_ID)
                   
                   return(temp)
                 })


output_2 <- output[!is.na(output)]

distance_frame <- rbindlist(output_2)
distance_frame$time_after_novelty <- distance_frame$time_after_novelty*-1

a <- as.data.frame(as.matrix(as.list(by(distance_frame$Distance, INDICES = distance_frame$time_after_novelty, FUN = mean))))
colnames(a) <- "Mean"
a$bins <- rownames(a)
b <- as.data.frame(as.matrix(as.list(by(distance_frame$Distance, INDICES = distance_frame$time_after_novelty, FUN = std.error))))
colnames(b) <- "Error"
data <- cbind(a,b)
data$Mean <- as.numeric(data$Mean) 
data$Error <- as.numeric(data$Error)

plot(data$Mean ~ data$bins, type = "l")

ses <- data$Mean + outer(data$Error, c(1,-1))
with(data, 
     plot(
       bins, Mean, type="l", ylim=c(0,1.00), ylab = "Mean Dissimilarity from Novel State (%)",
       xlab = "Years after Novel State", main = "Novelty evolution over time, bin width = 1",
       panel.first=polygon(c(bins,rev(bins)), c(ses[,1],rev(ses[,2])),border=NA, col="gold")
     )
)















### GAM of the invasives; correlate with novelty? ###

novelty_progression$Distance <- (novelty_progression$Distance * (length(novelty_progression$Distance)-1) + 0.5) / length(novelty_progression$Distance)

inv_gam <- gam(novelty_progression$Distance ~ s(novelty_progression$time_after_novelty, bs="cr", k= set.k), 
    family=betar(),
    method="REML")

inv.mu <- c(NA, inv_gam$fitted.values)
phi <- as.numeric(substr(inv_gam$family$family,
                         regexpr("\\(", inv_gam$family$family)+1,
                         nchar(inv_gam$family$family)-1))

# shape parameters used in qbeta.
A = inv.mu * phi
B = phi - A

# predict 5% and 95% prediction intervals from beta distribution parameters. 
# We use 95% rather than 97.5% because this is a one-tailed probability test.
# We are not concerned about communities that are MORE similar than predictions.
# This is done for each bin along the time-series.
inv.p <- do.call("rbind", lapply(1:length(A), function(n){
  data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
             upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
             inv.p = pbeta(novelty_progression$Distance[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
}))





par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,6,0.5,9), las=1)

# This plots theta's etc.. 
plot(novelty_progression$Distance ~ novelty_progression$time_after_novelty, type="n",
     ylim=c(max(novelty_progression$Distance, na.rm=TRUE)+0.1, 
            min(novelty_progression$Distance, na.rm=TRUE)),
     axes=FALSE, ylab = "", xlab = ""
     )

polygon(x=c(novelty_progression$time_after_novelty, rev(novelty_progression$time_after_novelty)),
        y=c(inv.p[,1][-1], rev(inv.p[,2][-1])), 
        col="grey80", border=NA)


lines(plogis(predict(inv_gam)) ~ novelty_progression$time_after_novelty, col="grey50", lwd=2)
lines(novelty_progression$Distance ~ novelty_progression$time_after_novelty)

points(y=novelty_progression$time_after_novelty[inv.p$inv.p< 0.05],
       x=novelty_progression$time_after_novelty[inv.p$inv.p < 0.05], pch=21, bg="green", col="darkgreen")
lims <- par("usr")

axis(side=2)
mtext(side=2, text = "Dissimilarity", line=3.5, las=0)

axis(side=1)
mtext(side=1, "Julian Days", line=2)


dev.off()



# Investigating characteristics of French Novel communities PCA

big_list # France
# basin_list 2080021030 # Per basin

basin_list <- as.list((subset(time_series_data, HydroBasin == 2080021030))$TimeSeriesID)



test <- GLM_lists_A$GLM_input_A_1
novels <- subset(test, novel == 1)

basin_total <- subset(test, site %in% basin_list)
novels_france <- subset(test, site %in% basin_list)

novels_true <- subset(novels, site %in% basin_list)



a <- Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A[]

# Remove text from the european communities
for (i in 1:length(a)) {
  
  if(is.character(a[[i]])) {
    a[[i]] <- NA
  }
  
}

a <- a[!is.na(a)]

# Adjust rownames
for (i in 1:length(a)) {
  print(i)
  
  matrix <- a[[i]]
  
  for (j in 1:nrow(matrix)){
    
    rownames(matrix)[j] <- paste0(rownames(matrix)[j], "-", names(a[i]))
    
  }
  
  a[[i]] <- matrix
  
}
  
  
# keep names after merging
names_list <- list()
for (i in 1:length(a)){
  print(i)
  names_list <- append(names_list, rownames(a[[i]]))
}  

# Bind the communities 
b <- as.data.frame(rbindlist(a, fill = TRUE))

# Assign the correct rownames
for (i in 1:nrow(b)){
  print(i)
  rownames(b)[i] <- (names_list)[[i]]
}




novels_france$coordinates <- paste0(novels_france$bins,"-", novels_france$site)
novels_true$coordinates <- paste0(novels_true$bins,"-", novels_true$site)


b$coordinates <- rownames(b)

novel_communities_france <- subset(b, coordinates %in% novels_france$coordinates)

# Get rid of NA's
novel_communities_france[is.na(novel_communities_france)] <- 0
b_2 <- b
b_2[is.na(b_2)] <- 0

# drop coordinates column for PCA
b <- b[,-which(names(b) %in% "coordinates")]

novel_communities_france <- as.data.frame(novel_communities_france)
novel_communities_france <- novel_communities_france[,-which(names(novel_communities_france) %in% "coordinates")]


b_2 <- novel_communities_france


novel_communities_france <- novel_communities_france %>% 
  select_if(colSums(.) != 0)



b_2  <- b_2  %>% 
  select_if(colSums(.) != 0)

for(i in 1:nrow(b_2)){
  print(i)
  if(rownames(b_2)[i] %in% novels_true$coordinates){
    b_2$groups[i] <- "novel"
  }
  else{
    b_2$groups[i] <- "non_novel"
  }
}


novels.pca <- prcomp(novel_communities_france, center = TRUE,scale. = TRUE)
summary(novels.pca)

novels.pca_b2 <- prcomp(b_2[,-which(names(b_2) %in% "groups")], center = TRUE,scale. = TRUE)
summary(novels.pca_b2)




library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(novels.pca_b3, var.axes = FALSE, alpha = 1)

ggscreeplot(novels.pca)

novels.pca_b3 <- prcomp(Fish_Communities_A$BioRealm_Matrices_A$palearctic_mat_A$G7646, center = TRUE,scale. = TRUE)
summary(novels.pca_b2)




