# This contains the plot_builder function and can be used to look at the abundance of any species in any timeseries over time. 
# It takes a species name and timeseries ID as input and produces a plot. It requires the files "time_series_data" and "Survey Data" 
# to work.

# Clear workspace (if you want) and set directory to access files.
rm(list = ls()) 
setwd("YOUR OWN CHOICE HERE")

# Required packages
library(lattice)
library(gridExtra)

# Load datafiles
time_series_data <- read.csv("1873_2_RivFishTIME_TimeseriesTable.csv")
Survey_Data <- read.csv("1873_2_RivFishTIME_SurveyTable.csv")
survey_labels <- read.csv("Unfiltered_Survey_Labels.csv")
filtered_labels <- read.csv("Filtered_Survey_Labels.csv")


# The plot-builder function which takes two input variables. Output plots abundance in whatever units
# the surveyers used over time in years. Red dots denote a survey, so anywhere where you see a red dot
# but no blue dot, there is an absence of that species (going to try and make this clearer but will do 
# for simple eyeballing).
plot_builder <- function(species_name, timeseries_code){
  
  i <- which(time_series_data$TimeSeriesID == timeseries_code)
  location <- time_series_data[i, 12]
  
  Time_Series_Data <- subset(Survey_Data, TimeSeriesID == timeseries_code)
  species_data <- subset(Time_Series_Data, Species == species_name)
  
  # I'm using 0.225 here as way of differentiating between surveys of the same location done in the same year but in different
  # quarters. I've averaged the ones that happened at the same time, I'm assuming it has something do to with multiple vessels.
  species_data$real_time <- species_data$Year + (as.numeric(species_data$Quarter)*0.225)
  Time_Series_Data$real_time <- Time_Series_Data$Year + (as.numeric(Time_Series_Data$Quarter)*0.225)
  agg_species <- aggregate(species_data$Abundance, by=list(species_data$real_time), FUN=mean)
  
  zero_mat <- matrix(ncol = 1, nrow = length(unique(species_data$real_time)))
  zero_vector <- data.frame(zero_mat)
  zero_vector[, 1] <- mean(species_data$Abundance)
  
  plot_species <- xyplot(agg_species$x ~ agg_species$Group.1, xlab = "Time (years)", ylab = paste("Abundance","(", Time_Series_Data[1, 7], ")") , type=c("p","l"), main=paste(species_name, "in Timeseries", timeseries_code, "(", location, ")"))
  
  plot_species <- update(plot_species, panel = function(...) {
    panel.xyplot(...)
    panel.xyplot(unique(Time_Series_Data$real_time), zero_vector$zero_mat, col = "red", pch = 9)
  })
  
  return(plot_species)
}

# Test for 4 random species in timeseries G935, on of our most complete ones. 
# (Feel free to have a look in the data and play around with the fishy names and timeseries)
# I might make a dataframe/list with species names etc. so we could automate this.
a <- plot_builder("Lates niloticus", "G935")
b <- plot_builder("Labeo senegalensis", "G935")
c <- plot_builder("Petrocephalus bovei", "G935")
d <- plot_builder("Anguilla anguilla", "G45")

# Create multi-plot. There will be a warning message generated, but that is only because the
# number of red points is greater than the blue points, as red denotes a survey, and blue denotes
# a record (some surveys have no record of a species)
grid.arrange(a, b, c, d, nrow = 2, ncol=2)


