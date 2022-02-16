### Trying to figure out if novelty is transient or not ###
### Using vegan functions to compute similarity between ###
### pre and post novelty states ###########################

library(tidyverse)
library(vegan)
library(dplyr)
library(data.table)

#########################
### 1. Pre-processing ###
#########################

# Unlist our old lists so we can combine them into one large metric
# for comparison.

ID_list_unlisted <- unlist(ID_list, 
                           recursive = FALSE,
                           use.names = FALSE)

matrices_unlisted <- unlist(Fish_Communities_A$BioRealm_Matrices_A, 
                            recursive = FALSE, 
                            use.names = FALSE)

names(matrices_unlisted) <- ID_list_unlisted


# Ensure quarters are made numeric. We jsut replace 2/3 for 2.3 as this approximates it quite well.
Survey_Data$Quarter <- gsub(Survey_Data$Quarter, 
                            pattern ="\\/", 
                            replacement = ".")

Survey_Data$Quarter <- as.numeric(Survey_Data$Quarter)

# Replace NA with 0 (this is what they signify)

Survey_Data$Quarter[is.na(Survey_Data$Quarter)] = 0


# Use the Julian abundance function to create community
# matrices.

abundance_julian_function <- function(survey_identifier){
  
  # Convert bin width to julian days.
  
  
  # Here I create some reference lists which help us create the matrices
  sub_series <- subset(Survey_Data, TimeSeriesID == survey_identifier)
  
  # Create a column which signifies Julian days.
  sub_series$True_Date <- NA
  sub_series$Julian_Day <- NA
  
  for (i in 1:nrow(sub_series)) {
    if (sub_series$Quarter[i] == 0){
      return(NA)
    }
    else{
      sub_series$True_Date[i] <- paste0(sub_series$Year[i],"-",ceiling((as.numeric(sub_series$Quarter[i])*(3))-1),"-15")
      sub_series$Julian_Day[i] <- JD(as.POSIXct(sub_series$True_Date[i]))
      sub_series$Julian_Day[i] <- floor(sub_series$Julian_Day[i])
    }
  }
  
  ref_list <- as.data.frame(unique(sub_series$Species))
  surveys_unique <- as.data.frame(unique(sub_series$SurveyID))
  year_unique <- as.data.frame(unique(sub_series$Julian_Day))
  
  
  
  
  # Here I create the empty dataframe which will be populated with abundance data
  
  abundance_df <- as.data.frame(matrix(data = NA, 
                                       nrow = length(unique(sub_series$Julian_Day)), 
                                       ncol = nrow(ref_list)))
  
  colnames(abundance_df) <- ref_list[, 1]
  rownames(abundance_df) <- year_unique[,1]
  
  # These loops will populate the empty data frame; it will simply add in the value from the main dataset into the
  # correct cell based on the row/column ID's.
  
  for (i in 1:nrow(abundance_df)) {
    sub_temp <- subset(sub_series, Julian_Day == rownames(abundance_df)[i])
    
    
    for (j in 1:ncol(abundance_df)) {
      abundance_list <- vector()
      
      for (k in 1:nrow(sub_temp)) {
        
        if ((colnames(abundance_df[j])) == toString(sub_temp[k, 5])){
          abundance_list <- append(abundance_list, sub_temp[k, 6])
          
        }
      }
      abundance_list <- as.numeric(abundance_list)
      average <- ceiling(mean(abundance_list))
      abundance_df[i,j] <- average
    }
  }
  
  # Order dataframe so that novelty framework does not get confused
  
  abundance_df <- abundance_df[order(as.numeric((row.names(abundance_df))), 
                                     decreasing = FALSE),]
  
  # Very small samples will be returned as vectors rather than dataframes, so will instruct the function to skip this survey if
  # That happens Shouldn't be necessary based on initial filter conditionals.
  
  if (typeof(abundance_df) == "double"){
    return("Not enough data")
  }
  
  # Change row names to bin (time)
  
  row.names(abundance_df) <- (as.numeric(row.names(abundance_df)))
  
  # Replace NA with 0 (this is what they signify)
  
  abundance_df[is.na(abundance_df)] = 0
  
  # Now remove rows that have sums of 0, these are rows of years that were actually not present in the original data.
  
  abundance_df <- abundance_df[rowSums(abundance_df[])> 0, ]
  
  
  return(abundance_df)
  
}

julian_compare_A_function <- function(check_list){
  
  for (i in 1:length(check_list)) {
    nam <- lapply(check_list[], 
                  function(TimeSeries_ID){
                    print(TimeSeries_ID)
                    temp <- abundance_julian_function(TimeSeries_ID)
                    if(class(temp)=="matrix"){
                      # remove bins within sites with fewer species than cut off
                      temp <- temp[rowSums(temp) > rich.cutoff,]
                    }
                    
                    return(temp)
                  }) 
  return(nam)
  }
}

# Our list of abundance matrices with timesteps of Julian days
julian_matrices_A <- julian_compare_A_function(ID_list_unlisted)
names(julian_matrices_A) <- ID_list_unlisted



############################
### 2. Run analysis ########
############################

# Write a function that compares similarity for a given
# community after novelty to the one preceding novelty.
dissimilarity_comparer_function <- function(matrix, 
                                            metric, 
                                            site,
                                            post){
  
    if(is.na(matrix)){
      return(NA)
    }
  
    if(!class(matrix)=="character"){

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

      # Now choose what you wish to compare based the position of the novel
      # community

      # Just need to control for novel communities at the very end of a timeseries,
      # as there is no following state.
  
    
      pre_novelty <- ((site.dist)) %>% 
        dplyr::select(starts_with("novel")) # This selects the community directly before
      
      novel.dist <- site.dist %>% filter(str_detect(rownames(site.dist), "novel"))
      
      
      

      # Control for no-novelty scenario's
      if(ncol(pre_novelty) == 0){
        return("No Novelty")
      }
    
      # Decides which post community we will be comparing
      if(((as.numeric(unlist(str_split(colnames(pre_novelty), 
                                     pattern = "-"))[2])) + post + 1) > ncol(site.dist)){
        return("Novelty at very end or not enough post")
      }
 
      # Get the index of the reference communities
      post_ref <- as.numeric(unlist(str_split(colnames(pre_novelty), 
                                             pattern = "-"))[2]) + 2 # in the case of immediately before and after novelty

      distance_comparison <- pre_novelty[post_ref,]

      # If there's two novel communities, let's just focus on 
      # the first one for now (we'll just treat the second one
      # as a new state when we come to that.
    
      if(typeof(distance_comparison) == "list"){
        distance_comparison <- distance_comparison[, ncol(distance_comparison)]
      }
    
      return(distance_comparison)
    }
    else{
      return(NA)
    }

}


# Write a function to execute everything at once
outcome_function <- function(check_list, 
                             post){
  
  checkout <- lapply(check_list[], 
                     function(TimeSeries_ID){
                       print(TimeSeries_ID)
                       temp <- dissimilarity_comparer_function(matrix = julian_matrices_A[[TimeSeries_ID]],
                                                               metric = "bray",
                                                               site = TimeSeries_ID,
                                                               post = post)
                      
                       return(temp)
                     }) 
  
}


# Write a function that returns a dataframe with 
# community distances for t + 1 and -1 from novelty.
transient_estimator <- function(check_list, 
                                post){
  
  # Outcome, let's run some diagnostics
  big_test <- outcome_function(check_list,
                               post = post)
  names(big_test) <- check_list[]

  # Bind into a dataframe
  big_test_2 <- do.call(rbind.data.frame, 
                        big_test)
  
  colnames(big_test_2) <- "Distance"

  # Remove all the useless entries (NA's, no novelty or
  # novelty in the final community)

  big_test_3 <- subset(big_test_2, Distance >= 0)

  new_df=as.data.frame(big_test_3[!grepl("No",big_test_3$Distance),])
  colnames(new_df) <- 'Distance'
  
  return(new_df)

}

# Useless! Use transition models instead.

# Run it

transient_analysis_post1 <- transient_estimator(ID_list_unlisted, 
                                                post = 1)

transient_analysis_post2 <- transient_estimator(ID_list_unlisted, 
                                                post = 2)

transient_analysis_post3 <- transient_estimator(ID_list_unlisted, 
                                                post = 3)

transient_analysis_post4 <- transient_estimator(ID_list_unlisted, 
                                                post = 4)

transient_analysis_post5 <- transient_estimator(ID_list_unlisted, 
                                                post = 5)

temp_list <- list(transient_analysis_post1, 
                     transient_analysis_post2,
                     transient_analysis_post3, 
                     transient_analysis_post4,
                     transient_analysis_post5)

transient_frame <- data.frame(Position = character(0), mean = numeric(0), SDERR = numeric(0),)

transient_frame <- as.data.frame(matrix(0, ncol = 3, nrow = 5 ))
colnames(transient_frame) <- c("position", "mean", "SDERR")

for (i in 1:5) {
  
  transient_frame$position[i] <- i
  transient_frame$mean[i] <- mean(as.numeric(temp_list[[i]]$Distance))
  transient_frame$SDERR[i] <- sd(as.numeric(temp_list[[i]]$Distance))/sqrt(length(temp_list[[i]]$Distance))
}

# Now compared to the actual novel community

transient_analysis_post1_N <- transient_estimator(ID_list_unlisted, 
                                                post = 1)

transient_analysis_post2_N <- transient_estimator(ID_list_unlisted, 
                                                post = 2)

transient_analysis_post3_N <- transient_estimator(ID_list_unlisted, 
                                                post = 3)

transient_analysis_post4_N <- transient_estimator(ID_list_unlisted, 
                                                post = 4)

transient_analysis_post5_N <- transient_estimator(ID_list_unlisted, 
                                                post = 5)

temp_list_N <- list(transient_analysis_post1_N, 
                  transient_analysis_post2_N,
                  transient_analysis_post3_N, 
                  transient_analysis_post4_N,
                  transient_analysis_post5_N)



transient_frame_N <- as.data.frame(matrix(0, ncol = 3, nrow = 5 ))
colnames(transient_frame_N) <- c("position", "mean", "SDERR")

for (i in 1:5) {
  
  transient_frame_N$position[i] <- i
  transient_frame_N$mean[i] <- mean(as.numeric(temp_list_N[[i]]$Distance))
  transient_frame_N$SDERR[i] <- sd(as.numeric(temp_list_N[[i]]$Distance))/sqrt(length(temp_list_N[[i]]$Distance))
}

############################
#### 3. Plot results #######
############################


ggplot(transient_frame, aes(x = position, y = mean)) +
  geom_point(size = 2, col = "red" ) + ylim(0, 0.75) +
  geom_errorbar(aes(ymin = mean - SDERR, ymax = mean + SDERR), size = 0.3) +
  labs(x = "Timepoints after novelty", y = "Community Dissimilarity (Bray-Curtis)", 
       title = "Dissimilarity of post novelty communities relative to pre-novelty state") + 
  theme_classic()
  
ggplot(transient_frame_N, aes(x = position, y = mean)) +
  geom_point(size = 2, col = "red" ) + ylim(0, 0.75) +
  geom_errorbar(aes(ymin = mean - SDERR, ymax = mean + SDERR), size = 0.3) +
  labs(x = "Timepoints after novelty", y = "Community Dissimilarity (Bray-Curtis)", 
       title = "Dissimilarity of post novelty communities relative to novel state") + 
  theme_classic()








