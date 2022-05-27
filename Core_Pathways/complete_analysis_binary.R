### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the presence/absence data pathway. Remember to load in functions, which are stored in separate files.
# This pathway also includes otpions for changing the bin width.
# You will not find indiviudal functions loaded in here; ensure that the "functions" is imported into Rstudio, this
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

# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")

# Assign individual time series to BioRealm groups

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
                australasia_ID )

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

# Create matrix list using this function
matrix_list_B_bin1 <- list_matrix_B_bins_function(ID_list, 
                                                  bin_width = 1)

# Create matrix_list with varying bin widths
matrix_list_B_bin2 <- list_matrix_B_bins_function(ID_list, 
                                                  bin_width = 2)

# Calculate novelty and return output in a list
novelty_list_B <- list_novelty_B_function(matrix_list_B_bin1)

novelty_list_B_bin2 <- list_novelty_B_function(matrix_list_B_bin2)

# Calculate probabilities and summarize results
# Use the required list summarizing novelty results as input for the analysis function

novelty_analysis_output_B <- novel.probability(novelty_list_B)

novelty_analysis_output_B_bin2 <- novel.probability(novelty_list_B_bin2)

# Create a final master list for Binary results
Fish_Communities_B <- list(ID_list, 
                           matrix_list_B_bin1, 
                           matrix_list_B_bin2, 
                           novelty_list_B, 
                           novelty_list_B_bin2, 
                           novelty_analysis_output_B, 
                           novelty_analysis_output_B_bin2 )

names(Fish_Communities_B) <- c("BioRealm_ID", 
                               "BioRealm_Matrices_B", 
                               "BioRealm_Matrices_B_2",  
                               "BioRealm_Novelty_B", 
                               "BioRealm_Novelty_B_2", 
                               "Analysis_Outputs_B", 
                               "Analysis_Outputs_B_2")

rm(matrix_list_B_bin1, 
   matrix_list_B_bin2, 
   novelty_list_B, 
   novelty_list_B_bin2, 
   novelty_analysis_output_B,
   novelty_analysis_output_B_bin2)

saveRDS(Fish_Communities_B, "./outputs/Fish_Communities_B.rds")

####### STEP 2 ##############
### ANALYZING THE RESULTS ###

# This function creates a "results frame" which will serve as input into
# our GLMM/GLM's. Automating this would be tedious so this bit is 
# slightly repetitive.

GLM_input_B_1 <- frame_builder_function(Fish_Communities_B, 
                                        bin_width = 1, 
                                        data_type = "B")

GLM_input_B_2 <- frame_builder_function(Fish_Communities_B, 
                                        bin_width = 2,
                                        data_type = "B")

GLM_lists_B <- list(GLM_input_B_1, 
                    GLM_input_B_2)

names(GLM_lists_B) <- c("GLM_input_B_1", 
                        "GLM_input_B_2")

# I have constructed two models to estimate probabilities whilst 
# taking into account covariate effects. The first model is a
# random intercept GLMM, treating site (timeseries_ID) as a random
# intercept. The second model excludes this random intercept as 
# in reality it did not explain any variance.

GLM_lists_B$GLM_output_Random_B_1 <- random_effects_GLMM(GLM_lists_B$GLM_input_B_1, bin_width = 1)

GLM_lists_B$GLM_output_Random_B_2 <- random_effects_GLMM(GLM_lists_B$GLM_input_B_2, bin_width = 2)

GLM_lists_B$GLM_output_Fixed_B_1 <- fixed_effects_GLM(GLM_lists_B$GLM_input_B_1, bin_width = 1)

GLM_lists_B$GLM_output_Fixed_B_2 <- fixed_effects_GLM(GLM_lists_B$GLM_input_B_2, bin_width = 2)

# Tidy up the environment
rm(GLM_input_B_1, GLM_input_B_2)

saveRDS(GLM_lists_B, "./outputs/GLM_lists_B.rds")


transition_data_A_1 <- estimate.observed.expected(prob.model.list = Fish_Communities_A$Analysis_outputs_A, 
                                                  novel.list = Fish_Communities_A$BioRealm_Novelty_A, 
                                                  dist.draws = 1e6)
transition_data_A_2 <- estimate.observed.expected(prob.model.list = Fish_Communities_A$Analysis_outputs_A_2, 
                                                  novel.list = Fish_Communities_A$BioRealm_Novelty_A_2, 
                                                  dist.draws = 1e6)



figure2.plot(trans.df = transition_data_A_2,
             plot.name = "test",
             ylims = log(c(0.1,15)))


