### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the "unaltered abundance" data pathway. Remember to load in functions, which are stored in separate files.
# This pathway also includes otpions for changing the bin width.
# You will not find individual functions loaded in here; ensure that the "functions" is imported into Rstudio, this
# folder contains all the essential functions.

######### STEP 1 ####################
### SETUP AND DATA TRANSFORMATION ###


# Clear environment (if wanted) and set your working directory
rm(list = ls())

setwd("YOUR OWN CHOICE HERE")

# source functions from 'functions' sub-folder
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

# Packages
install.packages(c("mgcv", "vegan", "lme4", "nlme", 
                   "DHARMa", "merTools", "shape",
                   "multcomp", "maptools", "sp", 
                   "divDyn", "plotrix", "raster",
                   "rgeos", "fun", "analogue",
                   "brms", "data.table"))

# Load data from RivFishTime
time_series_data <- read.csv("inputs/1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("inputs/1873_2_RivFishTIME_SurveyTable.csv")

# Assign individual time series to BioRealm groups.

palearctic_ID <- as.list(subset(time_series_data, BioRealm == "Palearctic")[,3])

nearctic_ID <- as.list(subset(time_series_data, BioRealm == "Nearctic")[,3])

afrotropics_ID <- as.list(subset(time_series_data, BioRealm == "Afrotropics")[,3])

neotropics_ID <- as.list(subset(time_series_data, BioRealm == "Neotropics")[,3])

australasia_ID <- as.list(subset(time_series_data, BioRealm == "Australasia")[,3])

# Create a list to hold these, serves as input for the matrix creator
ID_list <- list(palearctic_ID, 
                nearctic_ID, 
                afrotropics_ID, 
                neotropics_ID, 
                australasia_ID)

names(ID_list) <- c("palearctic_ID", 
                    "nearctic_ID", 
                    "afrotropics_ID", 
                    "neotropics_ID", 
                    "australasia_ID")

rm(palearctic_ID, 
   nearctic_ID, 
   afrotropics_ID, 
   neotropics_ID, 
   australasia_ID)

# Create matrix list of varying bin widths using this function

matrix_list_A_bin1 <- list_matrix_A_bins_function(ID_list, 
                                                  bin_width = 1)

matrix_list_A_bin2 <- list_matrix_A_bins_function(ID_list, 
                                                  bin_width = 2)






###### Incorporate seasonality #################


# Isolate a single time series
ID.list <- as.list(time_series_data$TimeSeriesID)

create.seasonal.ssmat <- function(n.ID){
   
   # Isolate time series data from overall survey data
   df<-subset(Survey_Data, TimeSeriesID == n.ID )
   
   # Create data frames for each timeseries-quarter pair
   output<-lapply(unique(df$Quarter), function(quarter){

      sub.1 <- subset(df, Quarter == quarter)
      
      sub.2<-as.data.frame(tapply(sub.1$Abundance,
                                   list(sub.1$Year, sub.1$Species),
                                   sum, na.rm = T))
      
      # Convert absences (NA) to 0's
      sub.2[is.na(sub.2)]<-0
      
      # Convert to relative abundance
      ss.mat <- sub.2/rowSums(sub.2)
      
      # Convert rownames to time in past
      rownames(ss.mat) <- 2021-as.numeric(rownames(ss.mat))
      
      # Filter time series with less than 10 timesteps and at least 2 species
      if(nrow(ss.mat) > 9 & ncol(ss.mat) > 1){
         return(ss.mat)
      }
      else{
         return(NULL)
      }
      
   })
   # Assign merger names for later use
   names(output) <- paste0(n.ID,".", unique(df$Quarter))
   return(output)
}

list_matrix_seasonality_function <- function(ID_list){
  ids <- unlist(ID_list)

  # Create community-by-species matrices for each season-timeseries pair
  mats<-lapply(ids, function(x){
     print(x)
     create.seasonal.ssmat(x)
  })
  
  # Homgenise list
  output<-unlist(mats, recursive = F)
  output.final<-Filter(Negate(is.null), output)
  
return(output.final)
}

matrix_list_seasonality <- list_matrix_seasonality_function(ID.list)

novelty.detection.gam <- function(input_list){
  
   # Keep track of iteration
   i <- 1
 
   nov.mat <- lapply(input_list , function(mat){
      
      # Simple message for the viewer
      print(paste0(i, " of ", length(input_list)))
      
      # provide a site name
      ID <- names(input_list)[i]
        
      # Run novelty detection framework
      temp <- identify.novel.gam(site.sp.mat = mat, alpha = 0.05, metric = "bray", site = ID, plot = TRUE, plot.data = FALSE,
                                    gam.max.k = -1)
      i <<- i + 1
      
      # Remove first 5 bins
      temp <- temp[-c(1:5),]
           
           
   return(temp)
   })

  return(nov.mat)
}

nov.output<-novelty.detection.gam(matrix_list_seasonality)

create.mod.frame <- function(novel.output){
  big.df<-rbindlist(nov.output)

  temp_df <- data.frame(do.call("rbind", strsplit(as.character(big.df$site), ".",
                                     fixed = TRUE)))
  big.df[,c("site")] <-NULL
  names(temp_df) <- c("site", "Quarter")

  return.frame<-cbind(temp_df, big.df)
  return(return.frame)
}

create.mod.frame(nov.output)
# Calculate novelty and return output in a list. No difference between Binary and Abundance data here except the similarity index used
# (Jaccard for binary, Bray-Curtis for abundance). Just remember to use the correct matrices as input (list_A for abundance, list_B for binary).

novelty_list_A <- list_novelty_A_function(matrix_list_A_bin1)

novelty_list_A_bin2 <- list_novelty_A_function(matrix_list_A_bin2)

# Calculate probabilities and summarize results
# Use the required list summarizing novelty results as input for the analysis function

novelty_analysis_output_A <- novel.probability(novelty_list_A)

novelty_analysis_output_A_bin2 <- novel.probability(novelty_list_A_bin2)

# Create a final master list for abundance results
Fish_Communities_A <- list(ID_list, 
                           matrix_list_A_bin1, 
                           matrix_list_A_bin2, 
                           novelty_list_A, 
                           novelty_list_A_bin2, 
                           novelty_analysis_output_A, 
                           novelty_analysis_output_A_bin2)

names(Fish_Communities_A) <- c("BioRealm_ID", 
                               "BioRealm_Matrices_A", 
                               "BioRealm_Matrices_A_2", 
                               "BioRealm_Novelty_A", 
                               "BioRealm_Novelty_A_2", 
                               "Analysis_outputs_A", 
                               "Analysis_outputs_A_2")

rm(matrix_list_A_bin1, 
   matrix_list_A_bin2, 
   novelty_list_A, 
   novelty_list_A_bin2, 
   novelty_analysis_output_A, 
   novelty_analysis_output_A_bin2)

saveRDS(Fish_Communities_A, "./outputs/Fish_Communities_A.rds")

####### STEP 2 ##############
### ANALYZING THE RESULTS ###

# This function creates a "results frame" which will serve as input into
# our GLMM/GLM's. Automating this would be tedious so this bit is 
# slightly repetitive.

GLM_input_A_1 <- frame_builder_function(Fish_Communities_A, 
                                        bin_width = 1, 
                                        data_type = "A")

GLM_input_A_2 <- frame_builder_function(Fish_Communities_A, 
                                        bin_width = 2,
                                        data_type = "A")

GLM_lists_A <- list(GLM_input_A_1, 
                     GLM_input_A_2)

names(GLM_lists_A) <- c("GLM_input_A_1", 
                         "GLM_input_A_2")

# I have constructed two models to estimate probabilities whilst 
# taking into account artefact (bin lag, length) effects. The first model is a
# random intercept GLMM, treating site (timeseries_ID) as a random
# intercept. The second model excludes this random intercept as 
# in reality it did not explain any variance.

GLM_lists_A$GLM_output_Random_A_1 <- random_effects_GLMM(GLM_lists_A$GLM_input_A_1, bin_width = 1)

GLM_lists_A$GLM_output_Random_A_2 <- random_effects_GLMM(GLM_lists_A$GLM_input_A_2, bin_width = 2)

GLM_lists_A$GLM_output_Fixed_A_1 <- fixed_effects_GLM(GLM_lists_A$GLM_input_A_1, bin_width = 1)

GLM_lists_A$GLM_output_Fixed_A_2 <- fixed_effects_GLM(GLM_lists_A$GLM_input_A_2, bin_width = 2)

# Tidy up the environment
rm(GLM_input_A_1, GLM_input_A_2)

# Save data file 
saveRDS(GLM_lists_A, "./outputs/GLM_lists_A.rds")


#### End of Novelty Detection


# Estimate the transition probabilities

transition_data_A_2 <- estimate.observed.expected(prob.model.list = Fish_Communities_A$Analysis_outputs_A_2, 
                           novel.list = Fish_Communities_A$BioRealm_Novelty_A_2, 
                           dist.draws = 1e6)

transition_data_B_2 <- estimate.observed.expected(prob.model.list = Fish_Communities_B$Analysis_Outputs_B_2, 
                                                  novel.list = Fish_Communities_B$BioRealm_Novelty_B_2, 
                                                  dist.draws = 1e6)



figure2.plot(trans.df = transition_data_A_2,
             plot.name = "test",
             ylims = log(c(0.1,15)))

figure2.plot.fixed <- function(trans.df, plot.name, xlims, ylims){
   
   cumul.col <- rgb(0.373,0.651,0.765)
   
   light.cols <- rbind(c(205,205,205),
                       c(207,228,237),
                       c(255,179,179),
                       c(255,228,179)) / 255
   light.cols <- rgb(light.cols[,1], light.cols[,2], light.cols[,3])
   
   library(shape)
   library(plotrix)
   
   #pdf(date.wrap(paste0("./plots/transition null model (", 
   #         plot.name, 
   #        ")"),
   # ".pdf"), 
   #  height=5.5, width=7.5, useDingbats = FALSE)
   
   framemat<-rbind(c(0.10,0.47,0.525,0.95),
                   c(0.47,0.84,0.525,0.95),
                   c(0.1,0.47,0.125,0.525),
                   c(0.47,0.84,0.125,0.525),
                   c(0.85,1,0.35,0.75))
   
   # framemat<-rbind(c(0.10,0.8,0.135,0.95),
   #                 c(0.8,1,0.25,0.75))
   # 
   split.screen(framemat)
   
   sapply(1:4, function(n){
      
      screen(n)
      par(mar=c(0,0,0,0), ps=10, tcl=-0.25, las=1, mgp=c(3,0.2,0))  
      
      sub.trans <- trans.df[trans.df$taxa == levels(trans.df$taxa)[n],]
      
      plot(x=NULL, y=NULL, xlim = xlims, ylim = ylims,
           axes=FALSE, xlab="", ylab="", yaxs="i")
      
      abline(h=0, lwd=2,col="grey60", lty="31")
      
      # custom.circle(x=log(1), y=log(1), r=0.5, col=rgb(0.5,0.5,0.5,0.15),
      #               screen=framemat[n,], border=NA)
      # 
      # draw.ellipse(x=log(0.05), y=log(0.775), a=0.85, b=0.35,
      #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
      #              angle= 0)
      # 
      # draw.ellipse(x=log(0.0002), y=log(2.8), a=2, b=1,
      #              border=NA, col=rgb(0.5,0.5,0.5,0.15),
      #              angle = 0)
      # 
      if(n == 3){
         mtext(side=1, line=2.25, at = par("usr")[2],
               text="Expected probability of transition\n(calculated from occurrence probabilities)")
      }
      
      if(n == 1){
         mtext(side=2, line=1.85, las=0, at = par("usr")[3],
               text="Observed probability / Expected probability")
      }
      
      # text(x=log(1), y=log(1.05), pos=3, adj=0.5, offset=0.8,
      #      labels="Background\nshifts", col="grey60", font=1)
      # 
      # text(x=log(0.05), y=log(0.75),  adj=0,
      #      labels="Shifts\nto and\nfrom\nnovelty", col="grey60", font=1)
      # 
      # text(x=log(0.0002), y=log(2.8),  adj=0,
      #      labels="Shifts\nbetween\nnovelty\ncategories", col="grey60", font=1)
      
      if(n %in% 3:4){
         axis(side=1, 
              at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
              labels=c(0.0001, 0.001,0.01,0.1,1,10,100))
      } else {
         axis(side=1, 
              at=log(c(0.0001,0.001,0.01,0.1,1,10,100)),
              labels=NA)
      }
      
      axis(side=1, at=log(c(seq(0.0001,0.001,0.0001),
                            seq(0.001,0.01,0.001),
                            seq(0.01,0.1,0.01),
                            seq(0.1,1,0.1),
                            seq(1,10,1),
                            seq(10,100,10))), tcl=-0.125, labels=NA)
      
      y.locs <- c(seq(1, 10, 1),
                  1/seq(1, 10, 0.5))
      
      axis(side=2, at=log(y.locs), tcl=-0.125, labels=NA)
      
      if(n %in% c(1,3)){
         axis(side = 2, 
              at = log(c(0.5,1,2,3,5,10,20)),
              labels = c(0.5,1,2,3,5,10,20), 
              mgp=c(3,0.5,0))
      } else {
         axis(side = 2, 
              at = log(c(0.5,1,2,3,5,10,20)),
              labels = NA, 
              mgp=c(3,0.5,0))
      }
      
      sub.trans <- sub.trans[order(sub.trans$non.zero),]
      
      sapply(1:dim(sub.trans)[1], function(x){
         print(x)
         
         segments(y0=log(sub.trans$ratio.mean[x]),
                  y1=log(sub.trans$ratio.mean[x]),
                  x0=log(sub.trans$exp.upper[x]),
                  x1=log(sub.trans$exp.lower[x]),
                  col=ifelse(sub.trans$non.zero[x],
                             "black", "grey70"))
         
         segments(y0=log(sub.trans$ratio.upper[x]),
                  y1=log(sub.trans$ratio.lower[x]),
                  x0=log(sub.trans$exp.mean[x]),
                  x1=log(sub.trans$exp.mean[x]),
                  col=ifelse(sub.trans$non.zero[x],
                             "black", "grey70"))
         
         if(sub.trans$non.zero[x]){
            
            aft.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.aft)[x]]
            bef.col = c("grey35", cumul.col, "red", "orange")[as.factor(sub.trans$cat.bef)[x]]
            shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
            border.col = "black"
            
         } else {
            
            aft.col = light.cols[as.factor(sub.trans$cat.aft)[x]]
            bef.col = light.cols[as.factor(sub.trans$cat.bef)[x]]
            shape = c("triangle", "square", "circle", "diamond")[sub.trans$taxa[x]]
            border.col = "grey70"
            
         }
         
         arrow.shape(x=log(sub.trans$exp.mean[x]),
                     y=log(sub.trans$ratio.mean[x]),
                     r=0.2, screen=framemat[1,], rads=c(1.5, 0.5),
                     col = aft.col,
                     shape=shape,
                     border=NA)
         
         arrow.shape(x=log(sub.trans$exp.mean[x]),
                     y=log(sub.trans$ratio.mean[x]),
                     r=0.2, screen=framemat[1,],
                     rads=c(0.5, 1.5), lwd=0.5,
                     shape=shape,
                     col = bef.col, 
                     border=border.col, add.arrow=TRUE)
         
         arrow.shape(x=log(sub.trans$exp.mean[x]),
                     y=log(sub.trans$ratio.mean[x]),
                     r=0.2, screen=framemat[1,], rads=c(0, 2),
                     shape=shape, border=border.col, add.arrow=FALSE)
         
         
      })
      
      text(x=relative.axis.point(0.02, "x"),
           y = relative.axis.point(0.935, "y"),
           labels = paste0("(",LETTERS[n],") ",
                           c("Diatoms", "Foraminifera", "Nannoplankton", "Radiolarians")[n]),
           font=2, adj=0)
      box()
      close.screen(n)
   })
   
   screen(5)
   par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
   plot.new()
   
   par(xpd=NA)
   
   taxa.pos <- rev(seq(0.375,0.675,len=4))
   
   text(x=0.575, y=0.75, adj=0.5,
        labels=bquote(bold(underline("Taxa"))), font=2, cex=0.8)
   sapply(taxa.pos, function(y){
      
      arrow.shape(x=0.125,
                  y=y,
                  r=0.2, screen=framemat[5,], rads=c(1.5,0.5),
                  col="white",
                  shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                  border="black", plot=TRUE)
      
      arrow.shape(x=0.125,
                  y=y,
                  r=0.2, screen=framemat[5,], rads=c(0.5, 1.5),
                  col="white",
                  shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                  border="black", plot=TRUE, add.arrow=TRUE, lwd=0.5)
      
      arrow.shape(x=0.125,
                  y=y,
                  r=0.2, screen=framemat[5,], rads=c(0,2),
                  shape=c("triangle", "square", "circle", "diamond")[which(taxa.pos==y)],
                  border="black", plot=TRUE)
      
      
   })
   # text(x=0.4, pos=4, y=taxa.pos,
   #      labels=c("Nanno", "Foram", "Radio", "Diatom"), 
   #      cex=0.8, offset=0.75)
   
   text(x=0.13, pos=4, y=taxa.pos,
        labels=sort(c("Nannoplankton", "Foraminifera", "Radiolarians", "Diatoms")), 
        cex=0.8, offset=0.75)
   
   par(lheight=0.85)
   text(x=0.15, y=1.05, labels="Preceding\ncommunity", adj=0, cex=0.8,
        lheight=0.5)
   
   text(x=0.85, y=0.925, labels="Succeeding\ncommunity", adj=1, cex=0.8,
        lheight=0.5)
   par(lheight=1)
   
   arrow.shape(x=0.55, y=1.2,
               r=0.25, screen=framemat[5,], rads=c(1.5,0.5),
               col="white", shape="circle",
               border="black", plot=TRUE)
   
   arrow.shape(x=0.45, y=1.2,
               r=0.25, screen=framemat[5,], rads=c(0.5,1.5),
               col="white", shape="circle",
               border="black", plot=TRUE, add.arrow = TRUE)
   
   segments(x0=c(0.125, 0.075, 0.875, 0.925),
            x1=c(0.075, 0.075, 0.925, 0.925),
            y0=c(1.05, 1.05, 0.925, 0.925),
            y1=c(1.05, 1.15, 0.925, 1.15))
   
   Arrows(x0=c(0.075, 0.925),
          x1=c(0.3, 0.7),
          y0=c(1.15, 1.15),
          y1=c(1.2, 1.2),
          arr.length = 0.1, arr.width = 0.1,
          arr.type = "triangle")
   
   rect.pos <- rev(seq(-0.275,0.05, len=4))
   
   rect(xleft=0.05, xright=0.2,
        ytop=rect.pos[1] + 0.03, ybottom=rect.pos[1] - 0.03, col="grey35")
   
   rect(xleft=0.05, xright=0.2,
        ytop=rect.pos[2] + 0.03, ybottom=rect.pos[2] - 0.03, col="red")
   
   rect(xleft=0.05, xright=0.2,
        ytop=rect.pos[3] + 0.03, ybottom=rect.pos[3] - 0.03, col=cumul.col)
   
   rect(xleft=0.05, xright=0.2,
        ytop=rect.pos[4] + 0.03, ybottom=rect.pos[4] - 0.03, col="orange")
   
   par(lheight=0.7)
   text(x=0.125, y=rect.pos - c(0, 0.015, 0.015, 0.015), pos=4, offset=0.75, cex=0.8,
        labels=c("Background", "Instantaneous\nnovelty", 
                 "Cumululative\nnovelty", "Novel\ncommunity"))
   par(lheight=1)
   
   text(x=0.575, y=0.215, adj=0.5,
        labels=bquote(bold(underline("Community"))), font=2, cex=0.8)
   text(x=0.575, y=0.15, adj=0.5,
        labels=bquote(bold(underline("classification"))), font=2, cex=0.8)
   
   par(xpd=FALSE)
   close.screen(5)
   
   dev.off()
}


figure2.plot.fixed(trans.df = transition_data_A_1,
             plot.name = "test",
             ylims = log(c(0.1,15)),
             xlims = c(-1,100))
dev.off()






