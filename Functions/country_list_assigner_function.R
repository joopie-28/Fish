# Small funtion thar extracts time series based on an ID

country_list_assigner <- function(country){
  
  country.lower <- tolower(country)
  
  a <- as.list(subset(time_series_data, Country == country)$TimeSeriesID)
  names(a) <- rep(subset(time_series_data, Country == country)$Country[1], times = length(a))
  
  return(a)
}