# Function for plotting custom forest plots for large glm's

forestPlotter <- function(model, labels){
  
  # Set theme
  set_theme(
    geom.outline.color = "antiquewhite4", 
    geom.outline.size = 1, 
    geom.label.size = 0.2,
    geom.label.alpha = 1,
    title.size = 1.5, 
    axis.angle.x = 0, 
    axis.textcolor = "black",
    axis.linecolor.x = 'black',
    axis.linecolor.y = 'black',
    base = theme_classic(),
    title.color = 'white'
    
  )
  
  # plot model
  plot_model(model, show.p = T,show.values = T, value.offset = 0.35, axis.labels =rev(labels),
             col = c('red', 'steelblue'), value.size = 4.2) +
    geom_hline(yintercept = 1, colour='black', lty = 2)
}