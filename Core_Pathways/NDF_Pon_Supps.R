###############################################
#### PoN Versus NDF supplementary analyses ####
###############################################

# J.M. Sassen #
# 27-01-2023  #

# Filter matrices where novelty occurred
nov.matrices <- matrix_list_seasonality[full.novel.mat.season$site[which(full.novel.mat.season$cat == 'novel')]]

PoN.supps <- list()


methods <- c('median', 'average', 'centroid','single', 'complete','ward.D')

PoN.supp <- merge_recurse(lapply(1:length(methods), function(i){
  #print(methods[i])
  PoN.supps[[i]] <- rbindlist(lapply(1:length(nov.matrices), function(x){
    
    print(paste0(methods[i],"_", x))
    cluster.method <- methods[i]
    # The start of a scoring approach
    label_frame <- identify.novel.gam(site.sp.mat = nov.matrices[[x]], 
                                          alpha = 0.05,
                                          metric = "bray",
                                          plot = T, 
                                          site = names(nov.matrices[x]),
                                          plot.data = FALSE,
                                          gam.max.k = -1)
    
    simprof_output <-tryCatch(simprof(data = nov.matrices[[x]], num.expected = 1000, undef.zero = T,
                                      num.simulated = 999, method.distance =vegdist,method.transform = "identity", 
                                      method.cluster = cluster.method, alpha=0.05, sample.orientation = "row"),error = function(e){
                                      "error"
                                      })
    if(any(simprof_output == "error")){
      output<- data.frame(paste0(names(nov.matrices[x]), "-", j ),
                                "Zero-Error")
      names(output) <- c('ID', methods[i]) 
      return(output)
    }
    
    novel_years <- as.numeric(label_frame$bins[which(label_frame$cat == 'novel')])
    
    # add a module to account for 'double novelty' events.
    temp.output <- list()
    for(j in 1:length(novel_years)){
      
      novel_year <- novel_years[j]
      
      clusters <- simprof_output$significantclusters
      
      # Find the novel cluster
      consistency <- FALSE
      for (list_index in 1:length(clusters)){
        if (novel_year %in% as.numeric(clusters[[list_index]]) & 
            novel_year == max(as.numeric(clusters[[list_index]]))){
          consistency <- TRUE
        }
      }
      
      # Consistency true, then we have a success
      
      temp.output[[j]] <- data.frame(paste0(names(nov.matrices[x]), "-", j ),
                                     consistency)
      output <- rbindlist(temp.output)
    }
    
    names(output) <- c('ID', methods[i]) 
    return(output)
    
  }))
  

}))

save(PoN.supp, ".\outputs")


test <- unique(PoN.supp)


index.list <- list()
for(i in 1:length((test))){
   index <-(strsplit(test[i], "-")[[1]][1])
   print(ncol(nov.matrices[[index]]))
   index.list <- c(index, index.list)
  }

errors <- list()
mats <- list()
for (index in index.list){
 errors[[index]] <-identify.novel.gam(site.sp.mat = nov.matrices[[2]], 
                     alpha = 0.05,
                     metric = "bray",
                     plot = T, 
                     site = index,
                     plot.data = FALSE,
                     gam.max.k = -1)
 mats[index] <- nov.matrices[index]
}


test<-unique(PoN.supp$ID[which(PoN.supp$median == "Zero-Error")])



#############################
### Solving Discrepancies ###
#############################

#SIMPROF

#Problem 1: Discrepancy between SIMPROF and Novelty Detection Framework due to very dissimilar timepoints early in the timeseries.
#Action: Experiment with the alpha parameter as a function of time series variability, or alternatively drop these time series.

# alpha needs to be a flexible parameter

placeholder <- identify.novel.gam.MDS(site.sp.mat = nov.matrices[[x]], 
                   alpha = 0.05,
                   metric = "bray",
                   plot = T, 
                   site = x,
                   plot.data = FALSE,
                   gam.max.k = -1)

# Find the average dissimilarities
mean_cum <- median(placeholder$raw.min.dist, na.rm =T)
mean_inst <- median(placeholder$seq.dist, na.rm =T)

if(mean_cum < 0.25){
  alpha_flex = 0.1
} else {
  if (mean_cum > 0.5){
    alpha_flex = 0.01
  }
  else{
    alpha_flex = 0.05
  }
}

dev.off()
simprof.plot(simprof(data = non_zero_offsetter(test), num.expected = 1000, undef.zero = F,
        num.simulated = 999, method.distance =vegdist, 
        method.cluster = "median", alpha=0.2, sample.orientation = "row"))

# We will probably just remove these I think.. 

# Need to build the consistency module into the persistence iteration 
# framework - v.7.


nov.cluster.id.V7 <- function(matrix, method_clus, 
                              alpha_clust){
  
  # Identify novelty in time series
  ID <- strsplit(names(matrix), "[.]" )[[1]][1]
  full.ID <- names(matrix)
  matrix <- matrix[[1]]
  
  number.names <- as.numeric(rownames(matrix))
  
  label_frame <- identify.novel.gam.MDS(site.sp.mat = matrix, 
                                        alpha = 0.05,
                                        metric = "bray",
                                        plot =F, 
                                        site = ID,
                                        plot.data = FALSE,
                                        gam.max.k = -1)
  
  # Find the average dissimilarities and matrix shapes
  mean_cum <- median(label_frame$raw.min.dist, na.rm =T)
  mean_inst <- median(label_frame$seq.dist, na.rm =T)
  nspp <- ncol(matrix)
  nyears <- nrow(matrix)
  
  #### Pre-processing Module ####
  
  if(as.numeric(rownames(matrix))[nrow(matrix)] <  as.numeric(rownames(matrix))[1]){
    
    # Flip such that orientation is correct 
    for(i in 1:length(matrix)){
      matrix[,i] <- rev(matrix[,i])
    }
    rownames(matrix) <- rev(rownames(matrix))
  }
  
  # Assign "background" state to first 5  bins
  
  for (i in (nrow(matrix)):(nrow(matrix)-4)) {
    rownames(matrix)[i] <- paste0("back-", rownames(matrix)[i])
    
  }
  
  # Assign actual states to the remaining bins, based on NDF
  
  for (i in 1:(dim(matrix)[1]-5)) {
    for (j in 1:(dim(label_frame)[1])) {
      
      if ((rownames(matrix)[i]) == (label_frame$bins)[j]){
        
        rownames(matrix)[i] <- paste0(label_frame$cat[j], "-", 
                                      rownames(matrix)[i])
      }
      
    } 
  } 
  

  # Obtain novel index from the data
  
  novel_frame <- matrix %>% 
    dplyr::filter(str_detect(rownames(matrix), "novel"))
  
  if(nrow(novel_frame) > 1){
    iter.list <- c(1,2)
    print('double event')
  }else{
    iter.list <- 1
  }
  
  
  return.dat <- lapply(iter.list, function(x){
    # Find year and row data
    novel_frame <- matrix %>% dplyr::filter(str_detect(rownames(matrix), "novel"))
    novel_frame <- novel_frame[x,]    
    
    novel_bin <<- as.numeric(strsplit(rownames(novel_frame), split = "-")[[1]][2])
    
    novel_index <- which(rownames(matrix) == paste0("novel-", novel_bin))
    
    # Run SIMPROF Routine on data matrix to identify multivariate structure of communities
    matrix.temp <- matrix
    rownames(matrix.temp) <- rev(number.names)
    
    # Sometimes using bray-curtis/czekanowski can lead to an error where there are a 
    # large number of zeroes in the data. We address this by ofsetting zeroes by 
    # a small number.
    
    test <<- tryCatch(
      simprof(data = matrix.temp, num.expected = 1000, undef.zero = T,
              num.simulated = 999, method.distance =vegdist, 
              method.cluster = method_clus, alpha=alpha_clust), 
      error=function(e) {
        print('ofsetting zeroes by small number')
        
        tryCatch(simprof(data = non_zero_offsetter(matrix.temp), num.expected = 1000, undef.zero = T,
                num.simulated = 999, method.distance =vegdist, 
                method.cluster = method_clus, alpha=alpha_clust),
                error=function(e){
                  NULL
                })
      })
    
    # Initialize parameters for the length calculations
    # Initialize a class variable
   
    Class <- "NONE"
    calc_length = T
    # Plot dendogram result
    par(mfrow = c(1,1))
    
    # If there is no support for significant clusters. We can not plot the dendrogram or 
    # do any meaningful length calculations.
    
    return.data = tryCatch(simprof.plot(test),
             error = function(e){
               print(paste0("No non-random structure in data for alpha = ", alpha_clust))
               data_error <- list(c(list('consistency' = FALSE,
                                          'median_seq_dis' = mean_inst,
                                          'median_cum_dis' = mean_cum,
                                          'nspp' =nspp,
                                          'nbins'= nyears)))
               calc_length <<- F
               names(data_error) <- full.ID
               return(data_error)
             })
    
   
    if(calc_length){
      # Start of the length calculation module
      
      
      length.output<-length.calculator.V2(test,novel_bin)
      
      length <- length.output$length 
      start <- length.output$start
      end <- length.output$end
      
      bin.start <- which(rev(number.names) == start)
      bin.end <- which(rev(number.names) == end)
      length.bins <- bin.start-bin.end
      
      
      
      # Initialize a class variable
      
      
      # Allocate blips needs to be done!!!
      if(start == end | length.bins == 1){
        
        Class <- "BLIP"
        length <- 0
        
        
      }
      
      # Test for full persistence
      if(Class != "BLIP"){
        
        Class <- "Persister"
        
      }
      
      if(novel_bin == number.names[length(number.names)]){
        Class <- "END"
        length <- 0
      }
      
      
      ### Consistency module - is the novel community the start of a new cluster
      clusters <- test$significantclusters
      
      # Find the novel cluster and confirm it is the head of its own cluster.
      consistency <- FALSE
      for (list_index in 1:length(clusters)){
        if (novel_bin %in% as.numeric(clusters[[list_index]]) & 
            novel_bin == max(as.numeric(clusters[[list_index]]))){
          consistency <- TRUE
        }
      }
      
      ### Module for confirming to which cluster the post-novel communities 
      # belong, but only possible if the novel comp ends before the ts.
      
      if(end != number.names[length(number.names)]){
       cluster_vector = length.output$clustervec
       postNov = cluster_vector[which(names(cluster_vector) == end)]
       preNov_clus = cluster_vector[1:(which(names(cluster_vector) == novel_bin)-1)]
       if(postNov %in% preNov_clus){
         FirstPostNov = "Pre_Novel_State"
       }else{
         FirstPostNov = "New_exploratory_State"
       }
      }else{
        FirstPostNov = NA
      }
      
      # Return data in df
      return.data <- list(c(list("ID" = ID, "begin" =start, 
                                 "end" =end, "length"= length, 
                                 "Class" = Class,
                                 "length.bins" = length.bins,
                                 'median_seq_dis' = mean_inst,
                                 'median_cum_dis' = mean_cum,
                                 'nspp' =nspp,
                                 'nbins'= nyears,
                                 'consistency' = consistency,
                                 'FirstPostNov' = FirstPostNov)))
      names(return.data) <- full.ID
    }
    return(return.data)
  })
  return(return.dat)
  
}


# Test run, need to do sinlgle and complete..

novelty.pers <- do.call(c, 
                        lapply(1:length(nov.matrices), function(m){
                          print(m)
                          nov.cluster.id.V7(nov.matrices[m],
                                            method_clus = 'single',
                                            alpha_clust = 0.05)
                          
                        })) 



# Write the results back into our main data frame.
full.novel.mat.season$novel.length <- NA
full.novel.mat.season$novel.class <- NA
full.novel.mat.season$consistency <- NA
full.novel.mat.season$median_seq_dis <- NA
full.novel.mat.season$median_cum_dis <- NA
full.novel.mat.season$nspp <- NA
full.novel.mat.season$FirstPostNov <- NA

for(i in 1:length(novelty.pers)){
  indices <- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel")
  if(length(indices) > 1){
    indices<- which(full.novel.mat.season$site == names(novelty.pers[[i]][1]) & full.novel.mat.season$cat == "novel" & full.novel.mat.season$bins == novelty.pers[[i]][[1]]$begin)
  }
  
  full.novel.mat.season$consistency[indices] <- novelty.pers[[i]][[1]]$consistency
  if(novelty.pers[[i]][[1]]$consistency){
    full.novel.mat.season$novel.length[indices] <- novelty.pers[[i]][[1]]$length.bins
    full.novel.mat.season$novel.class[indices] <-  novelty.pers[[i]][[1]]$Class
    full.novel.mat.season$median_seq_dis[indices] <- novelty.pers[[i]][[1]]$median_seq_dis
    full.novel.mat.season$median_cum_dis[indices] <- novelty.pers[[i]][[1]]$median_cum_dis
    full.novel.mat.season$nspp[indices] <- novelty.pers[[i]][[1]]$nspp
    full.novel.mat.season$FirstPostNov[indices] <- novelty.pers[[i]][[1]]$FirstPostNov
  }}
  


# Set non-novel lengths to 0 instead of NA
full.novel.mat.season$novel.length[is.na(full.novel.mat.season$novel.length)] <- 0
full.novel.mat.season$novel.class[is.na(full.novel.mat.season$novel.class)] <- "NONE"
full.novel.mat.season$novel.class <- as.factor(full.novel.mat.season$novel.class)

# Select all relevant data for modelling
class.frame <- full.novel.mat.season |>
  filter(consistency == TRUE) |>
  mutate(binary.class = ifelse(novel.class == 'Persister', 1, 0),
         delta_spp = (gain - loss),
         pers.proportion = novel.length/(total.n-position))

summary(glmer(binary.class ~ n.from.end+shannon.d+prop.loss+binary_invaders+(1|Quarter/site_ID), 
              data = class.frame, family = 'binomial')) # delta evenness?

persistence.frame <- class.frame |>
  filter(novel.class == 'Persister') 

# Need to integrate environemnt here
summary(lm(pers.proportion ~ diversity + evenness + prop.gain + binary_invaders, 
             data = persistence.frame)) # needs to beta or poission

# Model the 'next phase', might need a new modelling frame

summary(glm(firspostnov_binary ~ novel.class + position + n.from.end, 
            data = class.frame, family = 'binomial'))




# In what cluster is the n_end + 1 community?











# SOLVED #
#Problem 2: SIMPROF algorithm has issues with time series with a large proportion of 'monospecific rows' (i.e. complete dominance of 1 species).
#Action: Offsett zero values by a small number (0.001) to get the algorithm to run. Compare hclust tree (non-offset) to simprof tree (offset) and make sure the structure is the same.


test<-nov.matrices[["G9713.3"]]
x <- "G9713.3"

# Offset function for zero-errors

non_zero_offsetter <- function(matrix){
  
  for (row.index in 1:nrow(matrix)){
    
    row = matrix[row.index,]
    
    if (any(row == 1)){
      
      perfect.abundance <- which(row == 1)
      
      row[perfect.abundance] <- row[perfect.abundance] - 0.001*(length(row)-1)
      
      for(i in 1:length(row)){
        
        if (i != perfect.abundance){
          row[i] = row[i] + 0.001
        }
      }
    }
    matrix[row.index, ] <- row 
  }
  return(matrix)
}


non_zero_offsetter(test) # works well




##### BLIPS ####

# Anomaly indicators: mono-dominance, evenness.

# Does blip revert to new or previosuly seen cluster?

# just adapat the iteration framework i think 

# find 1. novel index and novel length. 2. find index of first
# post-novel community. 3. find cluster. 4. ompare cluster to previous
# 5. if new, then cluster is something diffeerent than before, if not
# we have return to status quo. 

test


###### Models ####











