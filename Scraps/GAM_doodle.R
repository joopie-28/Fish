### Trying to make the GAM ###


# Create position variable for each bin

test <- lapply(names(Fish_Communities_A$BioRealm_Novelty_A_2$palearctic_novelty_A), 
               function(x){
                 
                 temp <- as.data.frame(Fish_Communities_A$BioRealm_Novelty_A_2$palearctic_novelty_A[[x]])
                 temp$position <- 0
                 
                 for (i in 1:nrow(temp)) {
                   temp$position[i] <- (i + 5)
                  
                    }
                 return(temp)
                 
                 })

test <- rbindlist(test)

# Convert Boolean to binary values
tester <- test[, c("instant", "cumul", "novel") ] + 0
test[, c("instant", "cumul", "novel")] <- tester[,c("instant", "cumul", "novel")]


test$bins <- as.numeric(test$bins)
test$position <- as.numeric(test$position)

gam_test <- gam(novel ~ 1 + bin.lag + s(bins, bs = "tp"),
                family = betar(), data = test,
                method = "REML")
summary(gam_test)

plot.gam(gam_test)
log(bin.lag) + log(position) + log(n) + s(bins, bs = "tp", k = -1)







