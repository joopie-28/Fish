### Analyzing the novel communities using a slightly modified version of Pandolfi et Al's function. ###

# This function takes in a list of novelty detection matrices and returns a summary of novelty emergence. 

novel.probability <- function(novelty_list){
  
  # capture novel frequency for all taxa in list 
  novel.ts.freq <- do.call("rbind", 
                           lapply(1:5, function(n){
                             
                             x <- novelty_list[[n]]
                             
                             temp <- do.call("rbind", lapply(names(novelty_list[[n]][]),
                                                             function(ID){
                                                               
                                                               temp <- novelty_list[[n]][[which(names(novelty_list[[n]]) == ID)]]
                                                               
                                                               temp <- data.frame(TimeSeries_ID = ID,
                                                                                  t(sapply(c("back", "instant", "cumul", "novel"),
                                                                                           function(y){sum(temp$cat == y)})))
                                                               temp$taxa <- c("PAL", "NEA", "AFRO", "NEO", "AUS")[n]
                                                            
                                                               
                                                               return(temp)
                                                             }))
                           }))
 
  
  novel.ts.freq$taxa <- as.factor(novel.ts.freq$taxa)
  
  # We are running separate logistic regressions, rather than trying to run
  # multinomial models. This is because I think it's easier to talk about
  # raw probabilities of a particular classification, rather than
  # probabilities relative to a baseline category (which would have to be
  # 'background'). This is more intuitive to me.
  
  # modelling the probability of each classification occurring.
  novel.ts.freq$non.novel <- rowSums(novel.ts.freq[,c("back", "cumul", "instant")])
  novel.ts.freq$non.back <- rowSums(novel.ts.freq[,c("novel", "cumul", "instant")])
  novel.ts.freq$non.instant <- rowSums(novel.ts.freq[,c("novel", "cumul", "back")])
  novel.ts.freq$non.cumul <- rowSums(novel.ts.freq[,c("novel", "back", "instant")])
  
  ## modeling the probability of our two GAM tests.
  novel.ts.freq$all.instant <- rowSums(novel.ts.freq[,c("novel", "instant")])
  novel.ts.freq$all.cumul <- rowSums(novel.ts.freq[,c("novel", "cumul")])
  novel.ts.freq$non.all.instant <- rowSums(novel.ts.freq[,c("cumul", "back")])
  novel.ts.freq$non.all.cumul <- rowSums(novel.ts.freq[,c("instant", "back")])
  
  print("Fixed taxa models")
  fixed.prob.models <- fixed.taxa.prob.models(novel.freq.df = novel.ts.freq,
                                              test.model=TRUE)
  
  
  print("Random taxa models")
  random.prob.models <- random.taxa.prob.models(novel.freq.df = novel.ts.freq,
                                                test.model=TRUE)
  
  return(list(data = novel.ts.freq,
              fixed.prob.models = fixed.prob.models,
              random.prob.models = random.prob.models))
  
} 

