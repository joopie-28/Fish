# file to swiftly create csv's from model outputs

create_ouputCSV <- function(mod, name){
  pred.df <- round(as.data.frame(summary(mod)$coefficients),3)
  pred.df$Quarter.rand <- summary(mod)$varcor$Quarter[1,1]
  pred.df$Site.rand <- summary(mod)$varcor$site_ID[1,1]
  
  write.csv(pred.df,
            date.wrap(paste0("./outputs/Model_outputs/", 
                             name),
                      ".csv")) 
  print(paste0('csv stored in outputs'))
}



