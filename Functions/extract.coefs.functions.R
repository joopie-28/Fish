# Extract coefficients from glm/glmm

extract.coefs.glmm <- function(mod){
  pred.df <- as.data.frame(summary(mod)$coefficients)
  pred.df$site.rand <- summary(emig.mod)$varcor$site[1,1]
  pred.df$quarter.rand <- summary(mod)$varcor$Quarter[1,1]
  return(pred.df)
}
extract.coefs.glm <- function(mod){
  pred.df <- as.data.frame(summary(mod)$coefficients)
  return(pred.df)
}