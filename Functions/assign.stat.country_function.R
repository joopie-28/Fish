# Use this function to tag all species in a list of timeseries,
# it works on both basin and country data and takes a list of
# timeseries_ID's as input.

assign.stat.country <- function(country_list, country){
  
  output <- lapply(country_list[], 
                   function(TimeSeries_ID){
                     print((TimeSeries_ID))
                     
                     temp <- status_assignment_function(TimeSeries_ID, country)
                     
                     return(temp)
                   })
  
  for (i in 1:length(output)) {
    print(i)
    names(output)[i] <- country_list[[i]]
  }
  
  return(output)
  
}

assign.stat.country_nn <- function(country_list, country){
  
  output <- lapply(country_list[], 
                   function(TimeSeries_ID){
                     print((TimeSeries_ID))
                     
                     temp <- status_assignment_no_nn_function(TimeSeries_ID)
                     
                     return(temp)
                   })
  
  for (i in 1:length(output)) {
    print(i)
    names(output)[i] <- paste0(country_list[[i]],"-",names(country_list[i]))
  }
  
  return(output)
  
}

