### Analysis of global fish communities using the Novelty detection Framework (Pandolfi et al, 2020) ###

# This is the "unaltered abundance" data pathway. Remember to load in functions, which are stored in separate files.
# This pathway also includes otpions for changing the bin width.
# You will not find indiviudal functions loaded in here; ensure that the "functions" is imported into Rstudio, this
# folder contains all the essential functions.

# Clear environment (if wanted) and set your working directory
rm(list = ls())
#setwd("YOUR OWN CHOICE HERE")

# source functions from 'functions' sub-folder
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), 
       source)

# Packages
install.packages(c("mgcv", "vegan", "lme4", "nlme", 
                   "DHARMa", "merTools", "shape",
                   "multcomp", "maptools", "sp", 
                   "divDyn", "plotrix", "raster",
                   "rgeos", "fun", "analogue",
                   "brms"))

# Load data
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")

# Assign individual time series to BioRealm groups (if you have attempted to do the analysis using the presence/absence pathway
# earlier, this data will be loaded in already).

palearctic_ID <- as.list(subset(time_series_data, BioRealm == "Palearctic")[,3])

nearctic_ID <- as.list(subset(time_series_data, BioRealm == "Nearctic")[,3])

afrotropics_ID <- as.list(subset(time_series_data, BioRealm == "Afrotropics")[,3])

neotropics_ID <- as.list(subset(time_series_data, BioRealm == "Neotropics")[,3])

australasia_ID <- as.list(subset(time_series_data, BioRealm == "Australasia")[,3])

# Create a list to hold these, serves as input for the matrix creator
ID_list <- list(palearctic_ID, nearctic_ID, afrotropics_ID, neotropics_ID, australasia_ID )

names(ID_list) <- c("palearctic_ID", "nearctic_ID", "afrotropics_ID", "neotropics_ID", "australasia_ID")
rm(palearctic_ID, nearctic_ID, afrotropics_ID, neotropics_ID, australasia_ID)

# Create matrix list using this function
matrix_list_A <- list_matrix_A_function(ID_list)

# Create matrix_list with varying bin widths
matrix_list_A_bin2 <- list_matrix_B_bins_function(ID_list, 2)

# Calculate novelty and return output in a list. No difference between Binary and Abundance data here except the similarity index used
# (Jaccard for binary, Bray-Curtis for abundance). Just remember to use the correct matrices as input (list_A for abundance, list_B for binary).
novelty_list_A <- list_novelty_A_function(matrix_list_A)

novelty_list_A_bin2 <- list_novelty_A_function(matrix_list_A_bin2)

# Calculate probabilities and summarize results
# Use the required list summarizing novelty results as input for the analysis function

novelty_analysis_output_A <- novel.probability(novelty_list_A)

novelty_analysis_output_A_bin2 <- novel.probability(novelty_list_A_bin2)

# Create a final master list for abundance results
Fish_Communities_A <- list(ID_list, matrix_list_A, matrix_list_A_bin2, novelty_list_A, novelty_list_A_bin2, novelty_analysis_output_A, novelty_analysis_output_A_bin2)
names(Fish_Communities_A) <- c("BioRealm_ID", "BioRealm_Matrices_A", "BioRealm_Matrices_A_2", "BioRealm_Novelty_A", "BioRealm_Novelty_A_2", "Analysis_outputs_A", "Analysis_outputs_A_2")
rm(matrix_list_A, matrix_list_A_bin2, novelty_list_A, novelty_list_A_bin2, novelty_analysis_output_A, novelty_analysis_output_A_bin2)


# Create a Venn plot of model results

venn_plot_function(Fish_Communities_A$Analysis_Outputs_A)
venn_plot_function(Fish_Communities_A$Analysis_Outputs_A_2)




