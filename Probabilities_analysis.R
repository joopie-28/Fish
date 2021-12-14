### Moving on with analysis, this is very unfinished ###

estimate.observed.expected <- function(prob.model.list,
                                       novel.list,
                                       dist.draws = 1e6){
  
  # extract observed classification probabilities from models
  print("Preparing observed transition data...")
  
  # condense community data into single data-frame.
  comm.dat <- do.call('rbind', lapply(1:5, function(n){
    x <- novel.list[[3]][[n]]
    temp <- do.call("rbind", x)
    temp$taxa <- c("PAL", "NEA", "AFRO", "NEO", "AUS")[n]
    return(temp)
  }))
  
  # set site to be a taxa-specific factor
  comm.dat$site <- as.factor(paste0(comm.dat$taxa,":", comm.dat$site))
  
  # observed transition data-frame
  obs.df <- as.data.frame(table(comm.dat$cat.bef,
                                comm.dat$cat,
                                comm.dat$site))
  colnames(obs.df) <- c("cat.bef", "cat.aft", "site", "obs")
  obs.df$taxa <- substr(obs.df$site, 1, regexpr(":", obs.df$site)-1)
  obs.df$taxa.site <- paste0(obs.df$taxa, ":", obs.df$site)
  
  # get all transition combinations
  trans.list <- expand.grid(levels(obs.df$cat.bef), 
                            levels(obs.df$cat.aft),
                            stringsAsFactors = FALSE)
  
  # calculate expected probabilities from shift probability models
  print("Calculating expected transition probabilities...")
  base.probs <- do.call("rbind", lapply(prob.model.list$random.prob.models,
                                        function(x){x$pred.df}))
  
  rownames(base.probs) = 1:nrow(base.probs)
  base.probs <- as.data.frame(base.probs)
  base.probs$cat = c("all.instant", "all.cumul", "back", "instant", "cumul", "novel")
  base.probs <- base.probs[-(1:2),]

}

novel.list <- Fish_Communities_B # Remember
prob.model.list <- novelty_analysis_output_A # Remember

obs.exp.df <- estimate.observed.expected(prob.model.list = novelty_analysis_output_A, novel.list = Fish_Communities_B)

# Comparing observed and expected probability of transitions
obs.exp.df <- estimate.observed.expected(prob.model.list = novelty_analysis_output_A,
                                         novel.list = longhurst.novel.cut)

obs.exp.df$obs.Estimate.back <- plogis(obs.exp.df$obs.Estimate)
obs.exp.df$obs.lower <- plogis(obs.exp.df$obs.Estimate - 1.96 * obs.exp.df$`obs.Std. Error`)
obs.exp.df$obs.upper <- plogis(obs.exp.df$obs.Estimate + 1.96 * obs.exp.df$`obs.Std. Error`)

write.csv(obs.exp.df, date.wrap(paste0("./outputs/",
                                       " obs vs exp transition"),
                                ".csv"))

figure2.plot(trans.df = obs.exp.df, 
             plot.name = "longhurst_complete",
             ylims=log(c(0.35,15)))