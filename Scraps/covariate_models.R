### Creating a complete model using covariates such as bin lag and timeseries richness ####

# Create a giant data frame for individual communities

model_input <- frame_builder_function(Fish_Communities_B, 1, "B")


# This is a GLMM with site or timeseriesID treated as a random intercept
random_effects_GLMM <- function(novel.freq.df, 
                                    test.model=FALSE){

  temp_df <- novel.freq.df
  
  taxa.prob.models <- lapply(1:3, function(n){

    success_cat = c("instant", "cumul", "novel")[n]
    
    temp_df$success = temp_df[, ..success_cat]
    
    print("working on model 1")

    covariate_random_glmm <- glmer(success ~  log(bin.lag) + 
                                   log(TimeSeries_Length) + 
                                   (1|site) + 
                                   log(position) + 
                                   log(richness),
                                 data=temp_df, family=binomial)

# If the random intercept does not explain any variance, run as a GLM with fixed effects only instead.
    print("working on model 2")
    covariate_fixed_glmm <- glm(success ~  bin.lag + 
                                 TimeSeries_Length + 
                                 position + 
                                 richness,
                               data=temp_df, family=binomial)


# Remove non-significant factors
    print("working on model 3")
    covariate_dropped_glmm <- glm(success ~  log(bin.lag) + 
                                 log(position) + 
                                 log(richness),
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
  
  return(taxa.prob.models)
  
}




a <- random_effects_GLMM(GLM_input_B_1)

