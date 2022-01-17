### Creating a complete model using covariates such as bin lag and timeseries richness ####

# This is a GLMM with site or timeseriesID treated as a random intercept
random_effects_GLMM <- function(novel.freq.df, 
                                test.model=FALSE, 
                                bin_width){

  temp_df <- novel.freq.df
  
  taxa.prob.models <- lapply(1:3, function(n){

    success_cat = c("instant", "cumul", "novel")[n]
    
    temp_df$success = temp_df[, ..success_cat]
    
# If the random intercept does not explain any variance, run as a GLM with fixed effects only instead.
    print("Running random intercept GLMM")
    temp_df$site <- as.factor(temp_df$site)
    temp_df$bin.lag <- scale(temp_df$bin.lag, center = T, scale = T)
    temp_df$years_before <- scale(temp_df$years_before, center = T, scale = T)
    #temp_df$length <- scale(temp_df$length, center = T, scale = T)
    temp_df$richness <- scale(temp_df$richness, center = T, scale = T)
    #temp_df$TimeSeries_Length <- scale(temp_df$TimeSeries_Length, center = T, scale = T)
    temp_df$position <- scale(temp_df$position, center = T, scale = T)
    
    covariate_random_glmm <- glmer(success ~  (bin.lag) + 
                                 #(years_before) +
                                  position +
                                 #(length) +
                                 (1|site) +
                                 (richness),
                               data=temp_df, family=binomial)

  # group predictions
    print("summarizing")
  
    pred.df <- as.data.frame(summary(covariate_random_glmm)$coefficients)
    pred.df$taxa.rand <- summary(covariate_random_glmm)$varcor$taxa[1,1]
    pred.df$length.rand <- summary(covariate_random_glmm)$varcor$length[1,1]
  
    write.csv(pred.df,
            date.wrap(paste0("./outputs/", 
                             success_cat, " ", 
                             "covariate_random_glmm"),
                      ".csv")) 
    
  return(list(model=covariate_random_glmm,
              pred_df = pred.df))

})
  names(taxa.prob.models) <- c("instant", "cumul", "novel")
  return(taxa.prob.models)
  
}



