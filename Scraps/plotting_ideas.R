# plotting doodles

# Nice function to make a graph
novelty_plots <- function(frame){
  
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
  
  P <- ggplot() +
    geom_line(data = b[-36,] , aes(y = mean, x = year), colour = 'gold', size = 2) +
    xlab("Bin Age") + ylab('Probability of Novelty (%)') +
    geom_ribbon(data = b[-36,], aes(year, mean, ymin = mean - errors, ymax = mean + errors),
                fill = "gold", alpha = 0.3)
  
  scale_color_manual(name = "Data Type", values = c("Abundance" = "gold", "Presence/Absence" = "blue"))
  
  
  geom_errorbar(
    data = b[-30,],
    aes(year, mean, ymin = mean - errors, ymax = mean + errors),
    colour = 'black', 
    width = 0.4
  ) 
  
  
  
  
  return("Here's a nice overview of the data")
}

novelty_plots(com_frame)


