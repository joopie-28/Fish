# Fucntion that runs the effect of invader on novelty models

# Takes the feature data frame as input, and an additional arguement plot
# if you wish to plot the main result
inv.nov.glmm <- function(full.novel.mat, plot){

cat.var <- c("instant", "cumul", "novel")

# This function runs a GLM for each novelty category
invader.models <- lapply(cat.var, function(cat.1){
  
  # Had to do it this way, bit clunky but works fine
  full.novel.mat$response <-  full.novel.mat[,..cat.1]  
  
  print(paste0("Running ", cat.1, " invader model"))  
  
  # GLMM with interaction term and a nested random effect
  mod <- glmer(response ~ bin_lag + position + INC_increase + (1|country/basin),
               data = full.novel.mat, family = binomial)
  
  full.novel.mat$response <- NULL
  
  # Model Diagnostics with DHARMA
  disp.test <- testDispersion(simulateResiduals(mod))
  print(ifelse(disp.test$p.value <= 0.05, 
               paste0("Dispersal not OK ", round(disp.test$statistic, 2)),
               "Dispersal OK"))
  
  # Add coefficients to output to make inspection easier
  mod.1 <- list("model" = mod, "summary" = as.data.frame(summary(mod)$coef))
  
  return(mod.1)
  
})

names(invader.models) <- cat.var

if(plot){
  
  # if true, model plots for the invader effect are returned.
  model.plotter.6(invader.models, points = F)
}

return(invader.models)
}

# Plotting calls this function implicitly
model.plotter.6 <- function(invader.models, points){
  
  if (points){   
    plot_model(invader.models$model$novel, type = "pred",  terms=c("INC_increase [all]"), 
               title = "Probability of Novelty Emergence explained by Invader dynamics",
               axis.title = c("Net change in Invader Relative Abundance","Probability of Novel State (%)"),
               pred.type = "fe", colors = "green", show.data = T)
  }else{
    
    
    # Plot effects with ggplot
    
    df_novel <- ggpredict(invader.models$novel$model, type = "fe",  terms="INC_increase [all]")
    df_instant <- ggpredict(invader.models$instant$model, type = "fe" ,terms= "INC_increase [all]")
    df_cumul <- ggpredict(invader.models$cumul$model, type = "fe" ,terms= "INC_increase [all]")
    
    # Plot invader effect holding all else equal
    
    p1 <- ggplot(df_novel, aes(x, predicted)) +
      geom_line(color = "gold", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "yellow") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "True Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    p2 <- ggplot(df_instant, aes(x, predicted)) +
      geom_line(color = "red1", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "red3") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "Instantaneous Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    p3 <- ggplot(df_cumul, aes(x, predicted)) +
      geom_line(color = "steelblue4", lwd = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax = conf.high), alpha = 0.25, fill = "mediumturquoise") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.line = element_line("black"), panel.border = element_blank()) + 
      labs(y = "Cumulative Novel State (%)", x = "Invader Abundance Change") + 
      ylim(c(0:1))
    
    # 3-way plot
    grid.arrange(p1,p2,p3, nrow = 3)
  }
}