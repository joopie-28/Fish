# This function checks for overdispersion in glmm's using the dharma package

CheckDispersion <- function(mod){
  print("Running post-model tests...")
  sim.res <- simulateResiduals(mod)
  disp.test <- testDispersion(sim.res)
  print(ifelse(disp.test$p.value <=0.05, 
               paste0("Dispersal not okay, please inspect model. ", round(disp.test$statistic, 2)),
               "Dispersal okay :)"))
  return(sim.res)
}