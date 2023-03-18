# Extract coefficients from glm/glmm
extract.coefs <- function(mod.list){
  mod.name <- names(mod.list)
  
  # unlist object
  mod <- mod.list[[1]]
  if(class(mod)[1] == 'glmerMod'){
    pred.df <- as.data.frame(summary(mod)$coefficients)
    n.random.effects <- length(VarCorr(mod))
    for(i in 1:(n.random.effects)){
      variable_name <- paste0("Variance", " ",as.data.frame(VarCorr(mod))$grp[i])
      pred.df[variable_name] <- round(as.data.frame(VarCorr(mod))$sdcor[i], digits = 2)
     
    }

  } else{
    pred.df <- as.data.frame(summary(mod)$coefficients)
  }
  
  # Return it as a csv file in the outputs folder.
  print(paste0("Storing the model results in", paste0("./outputs/", mod.name, ".csv")))
  write_csv(pred.df, file = paste0("./outputs/", paste0(mod.name, ".csv")))
  
  return(pred.df)
}
