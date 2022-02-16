##### This is the script that includes the pre-processing steps #####

#### Characterizing the drivers of novel freshwater fish communities ####

###########################################################################
### This script includes main analyses as well as data import & tidying ###
### it requires functions to be loaded in from the function file, as well #
### as access to certain databases. #######################################
###########################################################################

# J.M. Sassen 31-01-2022 

##################################################
#### Step 0 . Load in packages and databases #####
##################################################

# Clear environment (if wanted) and set your working directory
rm(list = ls())

# source functions from 'functions' sub-folder and load them
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

package.loader(c("rgdal", "raster", "rgeos", "sf", 
                 "tidyverse", "rfishbase","mgcv",
                 "vegan", "lme4", "nlme", 
                 "DHARMa", "merTools", "shape",
                 "multcomp", "maptools", "sp", 
                 "divDyn", "plotrix", "raster",
                 "rgeos", "fun", "analogue",
                 "brms", "data.table", "data.table", "stringi", "tinytex", "knitr"))

# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv") # RIVFishTime 
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv") # RIVFishTime

invasives_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Occurrence_Table.csv") # Tedesco et Al. Invasives database
colnames(invasives_data) <- c("Basin", "Species", "Status", "TSN_ITIS_Code", "FishBase_code", "Valid_name", "Sighting") # Tedesco et Al. Invasives database

basin_data <- read.csv2("/Users/sassen/Desktop/datatoFigshare/Drainage_Basins_Table.csv") # Tedesco et Al. Basins geodatabase
colnames(basin_data) <- c("Basin", "Country", "BioRealm", "Endorheic", "Lon", "Lat", "Med_Lon", "Med_Lat", "Surface_Area") # Tedesco et Al. Basins geodatabase

occurence_shapefile <- read_sf(dsn = "/Users/sassen/Desktop/datatoFigshare/Basin042017_3119.shp") # Tedesco et AL. Shapefile

france_country_level <- france_species_status # Country level species data (Fishbase)


# Tidy data where necessary

invasives_data$Species <- gsub(invasives_data$Species, 
                               pattern ="\\.", 
                               replacement = " ")

invasives_data$Valid_name <- gsub(invasives_data$Valid_name, 
                                  pattern ="\\.", 
                                  replacement = " ")




##################################################
#### Step 1. Tagging species status  #############
##################################################

# Here, we will use a combination of the Tedesco et Al.
# database and our own custom species tags. 

# First we create a list of timeseries for all the countries we're 
# interested in.

fra.ID.list <- as.list(subset(time_series_data, Country == "FRA")$TimeSeriesID) # France
gbr.ID.list <- as.list(subset(time_series_data, Country == "GBR")$TimeSeriesID) # Great Britain
swe.ID.list <- as.list(subset(time_series_data, Country == "SWE")$TimeSeriesID) # Sweden
usa.ID.list <- as.list(subset(time_series_data, Country == "USA")$TimeSeriesID) # United States
esp.ID.list <- as.list(subset(time_series_data, Country == "ESP")$TimeSeriesID) # Spain
fin.ID.list <- as.list(subset(time_series_data, Country == "FIN")$TimeSeriesID) # Spain

# This function generates a dataframe matching a HydroBasin code to
# the name used in Tedesco et Al. This allows for tagging species
# on the basin level.

basin_name_code <- basin_name_match_function(occurence_shapefile)

# Use these data to add a basin_name column to the original 
# timeseries data

time_series_data$Basin_name <- NA
for (i in 1:nrow(time_series_data)) {
  print(i)
  for (j in 1:nrow(basin_name_code)) {
    if (time_series_data$HydroBasin[i] == basin_name_code$HydroBasin[j]){
      time_series_data$Basin_name[i] <- basin_name_code$Basin_name[j]
      
    }
  }
}

# Use this function to tag all species in a list of TimeSeries

#fra.stat.matrices <- assign.stat.country(fra.ID.list, country = "FRA") # uses Tedesco basin data as well as country data.
#gbr.stat.matrices <- assign.stat.country(gbr.ID.list, country = "GBR")
#swe.stat.matrices <- assign.stat.country(swe.ID.list, country = "SWE")
#usa.stat.matrices <- assign.stat.country(usa.ID.list, country = "USA")
#esp.stat.matrices <- assign.stat.country(esp.ID.list, country = "ESP")
#fin.stat.matrices <- assign.stat.country(fin.ID.list, country = "FIN")

# Using alternate method (no-non-natives, only invaders and natives)
fra.stat.matrices <- assign.stat.country_nn(fra.ID.list, country = "FRA") 
gbr.stat.matrices <- assign.stat.country_nn(gbr.ID.list, country = "GBR")
swe.stat.matrices <- assign.stat.country_nn(swe.ID.list, country = "SWE")
usa.stat.matrices <- assign.stat.country_nn(usa.ID.list, country = "USA")
esp.stat.matrices <- assign.stat.country_nn(esp.ID.list, country = "ESP")
fin.stat.matrices <- assign.stat.country_nn(fin.ID.list, country = "FIN")






###############################################################
#### Step 2. Computing invasive turnover metrics  #############
###############################################################

# Compute contribution of natives

fra.nnc.matrices <- mat.nnc.ass(fra.stat.matrices) # calculates relative contributions of natives and exotics.
gbr.nnc.matrices <- mat.nnc.ass(gbr.stat.matrices)
swe.nnc.matrices <- mat.nnc.ass(swe.stat.matrices)
usa.nnc.matrices <- mat.nnc.ass(usa.stat.matrices)
esp.nnc.matrices <- mat.nnc.ass(esp.stat.matrices)
fin.nnc.matrices <- mat.nnc.ass(fin.stat.matrices)


###############################################################
#### Step 3. Building a frame containing community metrics ####
###############################################################

fra.novel.mat <- inv.frame.builder(fra.nnc.matrices, country = "FRA") # this function returns the data frame for modelling.
gbr.novel.mat <- inv.frame.builder(gbr.nnc.matrices, country = "GBR")
swe.novel.mat <- inv.frame.builder(swe.nnc.matrices, country = "SWE")
usa.novel.mat <- inv.frame.builder(usa.nnc.matrices, country = "USA")
esp.novel.mat <- inv.frame.builder(esp.nnc.matrices, country = "ESP")
fin.novel.mat <- inv.frame.builder(fin.nnc.matrices, country = "FIN")

usa.novel.mat <- subset(usa.novel.mat, site != "G1040")

# Scale the data frame containing all countries (needs to be done after rbind..)
test <- rbind(fra.novel.mat, gbr.novel.mat, swe.novel.mat, usa.novel.mat, esp.novel.mat, fin.novel.mat)

test$NNC <- scale(test$NNC, center = T, scale = T)
test$NNC_increase <- scale(test$NNC_increase, center = T, scale = T)
test$NAC <- scale(test$NAC, center = T, scale = T)
test$NAC_increase <- scale(test$NAC_increase, center = T, scale = T)
test$INC <- scale(test$INC, center = T, scale = T)
test$INC_increase <- scale(test$INC_increase, center = T, scale = T)
test$bin_lag <- scale(as.numeric(test$bin_lag, center = T, scale = T ))
test$position <- scale(test$position, center = T, scale = T)
test$basin <- as.factor(test$basin)
test$site <- as.factor(test$site)
test[is.na(test)] <- 0

absolute_abundance_data_outlier <- test

saveRDS(binary_data, "./outputs/binary_data.rds")
saveRDS(relative_abundance_data, "./outputs/relative_abundance_data.rds")
saveRDS(Fish_Communities_A_ABS, "./outputs/Fish_Communities_A_ABS.rds")
saveRDS(absolute_abundance_data, "./outputs/absolute_abundance_data.rds")
saveRDS(absolute_abundance_data_outlier, "./outputs/absolute_abundance_data_outlier.rds")

#################################################################
#### Step 4. Modelling emergence of novelty by turnover #########
#################################################################

# I'm now thinking it makes more sense to use absolute abundance data

novel_model_ABS <- glmer(novel ~ bin_lag + position +NAC_increase*INC_increase +(1|country/basin),
                       data = absolute_abundance_data, 
                       family = binomial)


instant_model_ABS <- glmer(instant ~ bin_lag + position + NAC_increase*INC_increase+ (1|country/basin),
                           data = absolute_abundance_data, 
                           family = binomial)

cumul_model_ABS <- glmer(cumul ~ bin_lag + position + NAC_increase*INC_increase+ (1|country/basin),
                         data = absolute_abundance_data, 
                         family = binomial)




#################################################################
#### Step 5. Adding the novelty baseline metric ################
#################################################################


# Recalibration of the I.N.G function

identify.novel.gam.baseline <- function(site.sp.mat, alpha,
                               metric, site, plot=TRUE, plot.data=FALSE,
                               gam.max.k = -1){
  
  require(vegan) # community dissimilarity metrics
  require(mgcv) # additive modelling
  require(arm) # invlogit transformation
  
  # check to see if site species matrix rows run oldest -> youngest
  # if so, rotate site-species matrix to get oldest sites first
  if(as.numeric(as.character(rownames(site.sp.mat)[1])) < 
     as.numeric(as.character(rownames(site.sp.mat)[dim(site.sp.mat)[1]]))){
    
    site.sp.mat <- site.sp.mat[dim(site.sp.mat)[1]:1, ]
    
  }
  
  # calculate dissimilarity matrix, splitting those that can be
  # incorporated into vegdist from those that can't
  if(metric %in% c("chord", "SQchord")){
    
    require(analogue)
    
    site.dist <- as.matrix(distance(site.sp.mat,
                                    method=metric))
    
    if(metric=="chord"){site.dist <- site.dist / sqrt(2)}
    if(metric=="SQchord"){site.dist <- site.dist / 2}
    
  }
  
  if(metric == "hellinger"){
    
    hell.mat <- decostand(site.sp.mat, method="hellinger")  
    site.dist <- as.matrix(dist(hell.mat))
    site.dist <- site.dist / sqrt(2)
    
  }
  
  if(!metric %in% c("chord", "SQchord", "hellinger")){
    
    site.dist <- as.matrix(vegdist(site.sp.mat,
                                   method=metric))
  }
  
  # rescale distance matrices that scale beyond 1
  if(metric %in% c("euclidean", "SQchord", "chord")){
    site.dist <- site.dist / max(site.dist)
  }
  
  # check to make sure there are some non-0 or non-1 values
  # so we can actually estimate novelty
  if(var(as.vector(site.dist)) == 0){
    return(NA)
  }
  
  # convert bin names to numbers to use as continuous variable
  bins <- as.numeric(rownames(site.sp.mat))
  
  # This sets the maximum knot number for the additive model. One of the biggest
  # problems with additive models is overfitting. Ideally we want penalized 
  # splines, where knot number is weighted against variance explained. This will
  # provide a good measure of local mean trends in dissimilarity. But if we want
  # this method to run on smaller time-series, setting penalized splines doesn't work
  # very well. So I've set this to set maximum knots to half the number of bins, rounded
  # down, for time-series with between 4 and 10 bins. Less than 4 will not run using
  # splines, and is set to run as a linear regression below.
  set.k <- ifelse(dim(site.sp.mat)[1] > 10,
                  gam.max.k,
                  ifelse(dim(site.sp.mat)[1] > 4,
                         floor(dim(site.sp.mat)[1])/2,
                         NA))
  
  # NOVEL DISSIMILARITY #####
  # Obtain distance matrix for all sites in timeseries
  site.dist <- as.data.frame(as.matrix(vegdist(matrix,
                                               method=metric)))
  
  # Lets label!
  
  if (nrow(matrix) >= 10) {
    if (ncol(matrix) >=5){
      
      label_frame <- identify.novel.gam(site.sp.mat = Fish_Communities_A_ABS$BioRealm_Matrices_A_2$palearctic_mat_A$G8513, 
                                        alpha = 0.05,
                                        metric = metric,
                                        plot = TRUE, 
                                        site = site,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
    }
    else{
      return(NA)
    }
  }
  else{
    return(NA)
  }
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
  
  
  # Calculate difference between all states and the novel one
  
  novel.dist <- site.dist %>% filter(str_detect(rownames(site.dist), "novel"))
  
  novel.dist <- as.numeric(novel.dist)
  
  
  # Remove 0 and 1s for beta regression
  
  novel.dist.tr <- (novel.dist * (length(novel.dist)-1) + 0.5) / length(novel.dist)
  
  # model localised trend over time and extract residuals to get dissimilarity 
  # compared to local mean.
  
  if(var(novel.dist, na.rm=TRUE)==0){return(NULL)}
  
  if(!is.na(set.k)){
    novel.gam <- gam(novel.dist.tr ~ s(bins, bs="cr", k= set.k), 
                   family=betar(),
                   method="REML")
  } else{
    novel.gam <- gam(novel.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # COMPARING OBSERVED TO EXPECTED NOVEL BASELINE DISIMILARITY ####
  
  # This process calculates the p-value of the observed disimilarity score
  # being part of the expected distribution at the point in the time-series.
  
  # convert mu & phi beta parameters to A & B to use in qbeta.
  # mu = the mean predicted dissimilarity at each time point based on the 
  #      additive model.
  # phi = a dispersion parameter, which is set globally across the model.
  #       phi is unfortunately hidden deep in the gam model as a text string.
  
  novel.mu <- c(novel.gam$fitted.values)
  phi <- as.numeric(substr(novel.gam$family$family,
                           regexpr("\\(", novel.gam$family$family)+1,
                           nchar(novel.gam$family$family)-1))
  
  # shape parameters used in qbeta.
  A = novel.mu * phi
  B = phi - A
  
  # predict 5% and 95% prediction intervals from beta distribution parameters. 
  # We use 95% rather than 97.5% because this is a one-tailed probability test.
  # We are not concerned about communities that are MORE similar than predictions.
  # This is done for each bin along the time-series.
  novel.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               novel.p = pbeta(novel.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # SEQUENTIAL DISIMILARITY ####
  
  # calculate differenced sequential dissimilarities between 
  # time t and t-1
  seq.dist <- c(NA, diag(site.dist[-1,-dim(site.dist)[2]]))
  
  # transform to remove 0s and 1s for beta regression
  seq.dist.tr <- (seq.dist * (length(seq.dist)-1) + 0.5) / length(seq.dist)
  
  # model localised trend over time and extract residuals to get dissimilarity 
  # compared to local mean.
  if(var(seq.dist, na.rm=TRUE)==0){return(NULL)}
  
  if(!is.na(set.k)){
    seq.gam <- gam(seq.dist.tr ~ s(bins, bs="cr", k= set.k), 
                   family=betar(),
                   method="REML")
  } else{
    seq.gam <- gam(seq.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # MINIMUM DISIMILARITY ####
  
  # calculate minimum dissimilarity from time 1 to time t (time 1 being
  # earliest time point)
  min.dist <- sapply(1:dim(site.dist)[1], function(n){
    if(n==1){return(NA)}
    min(site.dist[n,1:n][-n])
  })
  
  # transform to remove 0s and 1s for beta regression
  min.dist.tr <- (min.dist * (length(min.dist)-1) + 0.5) / length(min.dist)
  
  if(var(min.dist, na.rm=TRUE)==0){return(NULL)}
  # model localised trend over time 
  if(!is.na(set.k)){
    min.gam <- gam(min.dist.tr ~ s(bins, bs="cr", k= set.k), 
                   family=betar(),
                   method="REML")
  } else{
    min.gam <- gam(min.dist.tr ~ bins, family=betar(), method="REML")
  }
  
  # COMPARING OBSERVED TO EXPECTED SEQUENTIAL DISIMILARITY ####
  
  # This process calculates the p-value of the observed disimilarity score
  # being part of the expected distribution at the point in the time-series.
  
  # convert mu & phi beta parameters to A & B to use in qbeta.
  # mu = the mean predicted dissimilarity at each time point based on the 
  #      additive model.
  # phi = a dispersion parameter, which is set globally across the model.
  #       phi is unfortunately hidden deep in the gam model as a text string.
  
  seq.mu <- c(NA, seq.gam$fitted.values)
  phi <- as.numeric(substr(seq.gam$family$family,
                           regexpr("\\(", seq.gam$family$family)+1,
                           nchar(seq.gam$family$family)-1))
  
  # shape parameters used in qbeta.
  A = seq.mu * phi
  B = phi - A
  
  # predict 5% and 95% prediction intervals from beta distribution parameters. 
  # We use 95% rather than 97.5% because this is a one-tailed probability test.
  # We are not concerned about communities that are MORE similar than predictions.
  # This is done for each bin along the time-series.
  seq.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               seq.p = pbeta(seq.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # EXPECTED OBSERVED TO MINIMUM DISIMILARITY ####
  
  # convert mu & phi beta parameters to A & B to use in qbeta
  min.mu <- c(NA, min.gam$fitted.values)
  phi <- as.numeric(substr(min.gam$family$family,
                           regexpr("\\(", min.gam$family$family)+1,
                           nchar(min.gam$family$family)-1))
  A = min.mu * phi
  B = (1 - min.mu) * phi
  
  # predict 5% and 95% prediction intervals from beta distribution parameters
  min.p <- do.call("rbind", lapply(1:length(A), function(n){
    data.frame(lwr=qbeta(0.05, shape1=A[n], shape2=B[n]),
               upr=qbeta(0.95, shape1=A[n], shape2=B[n]),
               min.p = pbeta(min.dist[n], shape1=A[n], shape2=B[n], lower.tail=FALSE))
  }))
  
  # CREATE RETURN DATA-FRAME ####
  return.data <- data.frame(site = site,
                            bins = rownames(site.sp.mat),
                            bin.lag = c(NA, abs(diff(as.numeric(rownames(site.sp.mat))))),
                            seq.dist = seq.dist,
                            seq.resid = c(NA, resid(seq.gam)),
                            raw.min.dist = min.dist,
                            min.resid = c(NA, resid(min.gam)),
                            seq.p = seq.p$seq.p,
                            min.p = min.p$min.p,
                            instant = seq.p$seq.p <= alpha,
                            cumul = min.p$min.p <= alpha,
                            novel = seq.p$seq.p <= alpha & min.p$min.p <= alpha,
                            seq.exp=seq.mu,
                            min.exp=min.mu,
                            seq.edf = summary(seq.gam)$s.table[1,"edf"],
                            min.edf = summary(min.gam)$s.table[1,"edf"],
                            n = nrow(site.sp.mat))
  
  attr(return.data$seq.exp, "seq.phi") = as.numeric(substr(seq.gam$family$family,
                                                           regexpr("\\(", seq.gam$family$family)+1,
                                                           nchar(seq.gam$family$family)-1))
  
  attr(return.data$min.exp, "min.phi") = as.numeric(substr(min.gam$family$family,
                                                           regexpr("\\(", min.gam$family$family)+1,
                                                           nchar(min.gam$family$family)-1))
  
  return.data$cat = "back"
  return.data[return.data$instant & 
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                !return.data$cumul, "cat"] = "instant"
  return.data[!return.data$instant & 
                !is.na(return.data$instant) &
                !is.na(return.data$cumul) &
                return.data$cumul, "cat"] = "cumul"
  return.data[return.data$instant & 
                !is.na(return.data$cumul) &
                !is.na(return.data$instant) &
                return.data$cumul, "cat"] = "novel"
  
  return.data$cat.bef <- c(NA, return.data$cat[-nrow(return.data)])
  
  # PLOT TIME SERIES ####
  
  # this section of code generates the plot seen when running the function.
  if(plot){
    par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,6,0.5,9), las=1)
    
    plot(seq.dist ~ bins, type="n",
         ylim=c(max(seq.dist, na.rm=TRUE)+0.1, 
                min(seq.dist, na.rm=TRUE)),
         axes=FALSE)
    
    polygon(x=c(bins, rev(bins)),
            y=c(seq.p[,1], rev(seq.p[,2])), 
            col="grey80", border=NA)
    lines(plogis(predict(seq.gam)) ~ bins[-1], col="grey50", lwd=2)
    lines(seq.dist ~ bins)
    points(y=seq.dist[seq.p$seq.p< 0.05],
           x=bins[seq.p$seq.p < 0.05], pch=21, bg="red")
    lims <- par("usr")
    
    sapply(which(return.data$novel), function(x){
      segments(x0 = bins[x],
               x1 = bins[x],
               y0 = seq.dist[x] + (0.05 * (par("usr")[3] - par("usr")[4])),
               y1 = par("usr")[3], col="orange", lwd=1)
    })
    
    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[4], y1=par("usr")[4])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
    axis(side=2)
    mtext(side=2, text = "Sequential\ndissimilarity", line=3.5, las=0)
    
    
    plot(min.dist ~ bins, type="n", axes=FALSE,
         ylim=c(min(min.dist, na.rm=TRUE), 
                max(min.dist, na.rm=TRUE)+0.1))
    polygon(x=c(bins, rev(bins)),
            y=c(min.p[,1], rev(min.p[,2])), 
            col="grey80", border=NA)
    lines(plogis(predict(min.gam)) ~ bins[-1], col="grey50", lwd=2)
    lines(min.dist ~ bins)
    
    points(y=min.dist[min.p$min.p< 0.05],
           x=bins[min.p$min.p < 0.05], pch=21, bg="skyblue")
    
    sapply(which(return.data$novel), function(x){
      segments(x0 = bins[x],
               x1 = bins[x],
               y0 = min.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
               y1 =par("usr")[4], col="orange", lwd=1)
    })
    
    plot(novel.dist ~ bins, type="n", axes=FALSE,
         ylim=c(min(novel.dist, na.rm=TRUE), 
                max(novel.dist, na.rm=TRUE)+0.1))
    polygon(x=c(bins, rev(bins)),
            y=c(novel.p[,1], rev(novel.p[,2])), 
            col="grey80", border=NA)
    lines(plogis(predict(novel.gam)) ~ bins, col="grey50", lwd=2)
    lines(novel.dist ~ bins)
    
    points(y=novel.dist[novel.p$novel.p< 0.05],
           x=bins[novel.p$novel.p < 0.05], pch=21, bg="skyblue")
    
    sapply(which(return.data$novel), function(x){
      segments(x0 = bins[x],
               x1 = bins[x],
               y0 = min.dist[x] + (0.05 * (par("usr")[4] - par("usr")[3])),
               y1 =par("usr")[4], col="orange", lwd=1)
    
    })
      
      
    par(xpd=NA)
    points(y=rep(par("usr")[4], sum(min.p$min.p< 0.05 & seq.p$seq.p < 0.05, na.rm=TRUE)),
           x=bins[min.p$min.p< 0.05 & seq.p$seq.p < 0.05 & !is.na(min.p$min.p)],
           pch=21, bg="orange")
    par(xpd=TRUE)
    
    segments(x0=par("usr")[1], x1=par("usr")[1], y0=par("usr")[3], y1=par("usr")[4])
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[3])
    segments(x0=par("usr")[2], x1=par("usr")[2], y0=par("usr")[3], y1=par("usr")[4])
    
    axis(side=1)
    mtext(side=1, "Time series age", line=2)
    
    axis(side=2)
    mtext(side=2, text = "Minimum\ndissimilarity", line=3.5, las=0)
    
    par(xpd=NA)
    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("", "", ""),
           pch=c(NA,NA,NA), lwd=c(NA,NA,1),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","orange"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)
    
    legend(x=par("usr")[2], y= par("usr")[4],
           legend=c("Instantaneous", "Cumulative", "Novel"),
           pch=c(21,21,21), lwd=c(NA,NA,NA),
           pt.bg=c("red","skyblue","orange"), col=c("black","black","black"),
           yjust=0.5, xjust=0, bty="n", seg.len=1, x.intersp=0.5)
    par(xpd=FALSE)
  }
  
  if(plot.data){return(list(return.data,
                            cbind(seq.p, c(NA, plogis(predict(seq.gam)))),
                            cbind(min.p, c(NA, plogis(predict(min.gam))))))
  } else {return(return.data)}
}

















# Model Diagnostics (dispersion test)
dev.off()

split.screen(c(2,1))
split.screen(c(1,2), screen = 1)

screen(2)

sim.res_1 <- simulateResiduals(novel_model_ABS)
disp.test_1 <- testDispersion(sim.res_1)

screen(3)
sim.res_2 <- simulateResiduals(instant_model_ABS)
disp.test_2 <- testDispersion(sim.res_2)

screen(4)
sim.res_3 <- simulateResiduals(cumul_model_ABS)
disp.test_3 <- testDispersion(sim.res_3)


print(ifelse(disp.test$p.value <=0.05, 
             paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),
             "Dispersal okay :)"))


plot(novel_model_ABS, 1)


# Summarize model output
pred.df_a <- as.data.frame(summary(instant_model)$coefficients)
pred.df_a$taxa.rand <- summary(instant_model)$varcor$taxa[1,1]
pred.df_a$length.rand <- summary(instant_model)$varcor$length[1,1]

# Plot models"INC_increase [all]"
library(sjPlot)

plot_model(novel_model_ABS, type = "pred",  terms=c("INC_increase [all]"), 
           title = "Probability of Novelty Emergence explained by Invader dynamics",
           axis.title = c("Net change in Invader Relative Abundance","Probability of Novel State (%)"),
           pred.type = "fe", colors = "green", show.data = T)

# Nicer to do it in ggPlot
df_novel <- ggpredict(novel_model_ABS, type = "fe",  terms="INC_increase [all]")
df_instant <- ggpredict(instant_model_ABS, type = "fe" ,terms= "INC_increase [all]")
df_cumul <- ggpredict(cumul_model_ABS, type = "fe" ,terms= "INC_increase [all]")

df_novel_nac <- ggpredict(novel_model_ABS, type = "fe",  terms="NAC_increase [all]")
df_instant_nac <- ggpredict(instant_model_ABS, type = "fe",  terms="NAC_increase [all]")
df_cumul_nac <- ggpredict(cumul_model_ABS, type = "fe",  terms="NAC_increase [all]")

# Plot INC holding all else equal
p1 <- ggplot(df_novel, aes(x, predicted)) +
  geom_line(color = "gold", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "yellow") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
          labs(y = "True Novel State (%)", x = "Invader Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(novel_model_ABS)$coefficients)[5,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")
  
p2 <- ggplot(df_instant, aes(x, predicted)) +
  geom_line(color = "red1", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "red3") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Instantaneous Novel State (%)", x = "Invader Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(instant_model_ABS)$coefficients)[5,4], 
                                                                           digits = 5), " *"), family = "Times New Roman")

p3 <- ggplot(df_cumul, aes(x, predicted)) +
  geom_line(color = "steelblue4", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "mediumturquoise") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Cumulative Novel State (%)", x = "Invader Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(cumul_model_ABS)$coefficients)[5,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")


p4 <- ggplot(df_novel_nac, aes(x, predicted)) +
  geom_line(color = "gold", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "yellow") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "True Novel State (%)", x = "Native Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(novel_model_ABS)$coefficients)[4,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")

p5 <- ggplot(df_instant_nac, aes(x, predicted)) +
  geom_line(color = "red1", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "red3") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Instantaneous Novel State (%)", x = "Native Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(instant_model_ABS)$coefficients)[4,4], 
                                                                         digits = 5), ""), family = "Times New Roman")

p6 <- ggplot(df_cumul_nac, aes(x, predicted)) +
  geom_line(color = "steelblue4", lwd = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "mediumturquoise") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line("black"), panel.border = element_blank()) + 
  labs(y = "Cumulative Novel State (%)", x = "Native Abundance Change") + 
  ylim(c(0:1)) + geom_text(x = -7, y = 0.8, label = paste0("p = ", round(as.data.frame(summary(cumul_model_ABS)$coefficients)[4,4], 
                                                                         digits = 5), " ***"), family = "Times New Roman")
dev.off()



grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)





###########################################
###### Step X. HydroBasin Aggregation #####
###########################################


hydrosheds <- read_sf(dsn ="/Users/sassen/Desktop/HydroRIVERS_v10_eu_shp/HydroRIVERS_v10_eu_shp")

hydrobasins_12 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev12_v1c")

hydrobasins_8 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev08_v1c")

hydrobasins_0 <- read_sf(dsn = "/Users/sassen/Desktop/hybas_eu_lev00_v1c")



# Attempt to match geometry

loc_points <- time_series_data[, c("Longitude", "Latitude", "TimeSeriesID")]

# Converts lat-long coordinates to geometry
loc_points_sf <- st_as_sf(loc_points, 
                          coords = c("Longitude", "Latitude"), 
                          crs = st_crs(hydrobasins_12))

# Set ID to character instead of integer
hydrobasins_12$HYBAS_ID <- as.character(hydrobasins_12$HYBAS_ID)


# Essential line, functions will not run without this setting 
sf::sf_use_s2(FALSE)

# Create a matrix containing HydroBasin ID and Tedesco Basin name
# This piece of code will find the intersection of the coordinates
# and the named basin.

loc_points <- loc_points_sf %>% mutate(
  
  intersection = as.integer(st_intersects(geometry, 
                                          hydrobasins_12)),
  
  hybas_id = if_else(is.na(intersection), '', 
                 hydrobasins_12$HYBAS_ID[intersection])
)

# Store timeseriesID - subbasin combo's in a 'by' list
loc_by_subbasin <- by(loc_points$TimeSeriesID, INDICES = loc_points$hybas_id, FUN = list)

# Filter out the French ones

france_ID <- subset(time_series_data, Country == "FRA")$TimeSeriesID

for (i in 1:length(loc_by_subbasin)){
  print(i)
  for(j in 1:length(loc_by_subbasin[[i]])){
    if(loc_by_subbasin[[i]][j] %notin% france_ID){
      loc_by_subbasin[[i]] <- NA
      
    }
  }
}

# Clean up which leaves the French ones

loc_by_subbasin <- loc_by_subbasin[!is.na(loc_by_subbasin)]



totals_inv <- read_csv("/Users/sassen/Desktop/Invaders_country.csv")
totals_inv$totals <- 0
totals_inv$totals[1] <- nrow(france_country_level_nn)
totals_inv$totals[2] <- nrow(usa_country_level_nn)
totals_inv$totals[3] <- nrow(gbr_country_level_nn)
totals_inv$totals[4] <- nrow(spain_country_level_nn)
totals_inv$totals[5] <- nrow(swe_country_level_nn)
totals_inv$totals[6] <- nrow(finland_country_level_nn)


ggplot(totals_inv) +
  geom_bar(aes(x = Country, y = totals), stat = "identity", colour = "green", fill = "green") +
  scale_color_grey() + theme_classic() +
  geom_bar(aes(x = Country, y = Invaders), stat = "identity", colour = "red", fill = "red") + 
  ylab("Unique species per country")

identify.novel.gam(site.sp.mat = Fish_Communities_A_ABS$BioRealm_Matrices_A_2$palearctic_mat_A$G10853,
                   site = "NA", 
                   metric = "bray", 
                   gam.max.k = -1,
                   alpha = 0.05)
  
  


