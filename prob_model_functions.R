### Analysis of novelty detection as per Pandolfi et al 2020 ###

fixed.taxa.prob.models <- function(novel.freq.df, 
                                   test.model=FALSE){
  
  require(multcomp)
  
  taxa.prob.models <- lapply(1:6, function(n){
    
    success.var = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")[n]
    failure.var = c("non.all.instant", "non.all.cumul", "non.back", 
                    "non.instant", "non.cumul", "non.novel")[n]
    
    print(success.var)
    
    temp.df <- novel.freq.df  
    temp.df$success = temp.df[, success.var]
    temp.df$failure = temp.df[, failure.var]
    
    taxa.prob.m <- glm(cbind(success, failure) ~ -1 + taxa,
                       data=temp.df, family=binomial)
    
    # model tests
    if(test.model){
      print("Running post-model tests...")
      sim.res <- simulateResiduals(taxa.prob.m)
      disp.test <- testDispersion(sim.res)
      print(ifelse(disp.test$p.value <=0.05, 
                   paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),
                   "Dispersal okay :)"))
      
      gc()
    }
    
    # group comparisons
    pair.comp <- glht(taxa.prob.m, linfct=mcp(taxa = "Tukey"))
    pair.conf <- as.data.frame(summary(confint(pair.comp))$confint)
    pair.comp <- summary(pair.comp)
    
    pairs.df <- data.frame(test=names(pair.comp$test$coefficients),
                           coefs=pair.comp$test$coefficients,
                           sigma=pair.comp$test$sigma,
                           t=pair.comp$test$tstat,
                           p=pair.comp$test$pvalues,
                           lower = pair.conf$lwr,
                           upper = pair.conf$upr)
    
    write.csv(pairs.df,
              date.wrap(paste0(success.var,
                               " fixed taxa probability comparison"),
                        ".csv")) 
    
    write.csv(summary(taxa.prob.m)$coefficients,
              date.wrap(paste0(success.var,
                               " fixed taxa model summary"),
                        ".csv")) 
    
    # group predictions
    pred.df <- data.frame(taxa=factor(levels(novel.freq.df$taxa),
                                      levels = levels(novel.freq.df$taxa)),
                          var = success.var)
    
    pred.df <- cbind(pred.df,
                     as.data.frame(predict(taxa.prob.m, 
                                           newdata=pred.df, se.fit=TRUE)))
    
    return(list(model=taxa.prob.m,
                pred.df = pred.df,
                comp.df = pairs.df))
    
  })
  
  return(taxa.prob.models)
  
}

random.taxa.prob.models <- function(novel.freq.df, 
                                    test.model=FALSE){
  
  require(multcomp)
  require(lme4)
  require(merTools)
  
  taxa.prob.models <- lapply(1:6, function(n){
    
    success.var = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")[n]
    failure.var = c("non.all.instant", "non.all.cumul", "non.back", 
                    "non.instant", "non.cumul", "non.novel")[n]
    
    print(success.var)
    
    temp.df <- novel.freq.df
    temp.df$success = temp.df[, success.var]
    temp.df$failure = temp.df[, failure.var]
    
    taxa.prob.m <- glmer(cbind(success, failure) ~ 1 + (1|taxa),
                         data=temp.df, family=binomial)
    
    # model tests
    if(test.model){
      print("Running post-model tests...")
      sim.res <- simulateResiduals(taxa.prob.m)
      disp.test <- testDispersion(sim.res)
      print(ifelse(disp.test$p.value <=0.05, 
                   paste0("Dispersal not okay :( -- ", round(disp.test$statistic, 2)),
                   "Dispersal okay :)"))
      
      gc()
    }
    
    # group predictions
    pred.df <- as.data.frame(summary(taxa.prob.m)$coefficients)
    pred.df$taxa.rand <- summary(taxa.prob.m)$varcor$taxa[1,1]
    
    write.csv(pred.df,
              date.wrap(paste0(success.var,
                               " random taxa model summary"),
                        ".csv")) 
    
    return(list(model=taxa.prob.m,
                pred.df = pred.df))
    
  })
  
  return(taxa.prob.models)
  
}

# A little custom function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}
