# plotting doodles

# Nice function to make a graph, it's a bit convoluted but 
# it all works

a <- by(data = frame$novel, INDICES = frame$bins, FUN = mean)
errors <- by(data = frame$novel, INDICES = frame$bins, FUN = std.error)

  
b <- sapply(names(a), function(x){
  l <- a[[x]]
  return(l)
})
  
errors <- sapply(names(a), function(x){
  l <- errors[[x]]
  return(l)
})
  

b <- as.data.frame(b)
errors <- as.data.frame(errors)
a <- 2021-as.numeric(names(a))
 
 
b <- as.data.frame(cbind(b$b, errors$errors, a))
colnames(b) <- c("mean", "errors","year")
  
# Plot it!
ggplot() +
  geom_point(data = b, aes(y = mean, x = year), colour = 'gold', size = 2) +
  xlab("Year") + ylab('Probability of Novelty (%)') +
  #geom_ribbon(data = b, aes(year, mean, ymin = mean - errors, ymax = mean + errors),
              #fill = "gold", alpha = 0.3)
    geom_errorbar(
    data = b,
    aes(year, mean, ymin = mean - errors, ymax = mean + errors),
    colour = 'black', 
    width = 0.4
  ) 
  
  
  
  
 
  
  
  
  
  





