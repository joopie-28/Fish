### nMDS ordination of communities within the dataset ###

library(vegan)

### NOTE, all the variables here are named france because that's how I started it.
### they're just placeholder variables.

# First we create a commmunity by species matrix for the countries we're
# looking at.

nMDS_frame_function <- function(HydroBasin_ID){

  
  france_ID <- subset(time_series_data, HydroBasin == HydroBasin_ID)$TimeSeriesID

  france_matrices <- list()

  # Collect all the community matrices that belong to France

  for (i in 1:length(Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A)){
    name_list_2 <- list()
  
    if (names(Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A[i]) %in% france_ID){
    
      matrix <- Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A[i]
   
      if (typeof(matrix[[1]]) == "character"){
        matrix[[1]] <- NA
      }
    
      else {
        # Remove as we do not consider these. 
        matrix[[1]] <- matrix[[1]][-c(1:5), ]
        matrix[[1]]$category <- NA
      
        # Replace rownames with the novelty category so we can do a meaningful
        # visualization.
      
        identifier <- names(Fish_Communities_A$BioRealm_Matrices_A_2$palearctic_mat_A[i])
        n_m <- Fish_Communities_A$BioRealm_Novelty_A_2$palearctic_novelty_A[[identifier]]
      
        if (typeof(n_m) != "list"){
          next
        }
      
        else {  
          for (y in 1:nrow(matrix[[1]])) {
            for (k in 1:nrow(n_m)) {
          
              if (n_m$bins[k] == rownames(matrix[[1]][y,]))
              {
                name_list_2 <- append(name_list_2, 
                                      paste0(identifier,"_", n_m$cat[k],"_", n_m$bins[k]))
                
                matrix[[1]]$category[y] <- n_m$cat[k]
              }
            }
          }
        
          rownames(matrix[[1]]) <- name_list_2
        }
      }
      france_matrices <- c(france_matrices, matrix)
    }
  }


  # Clean up the data and combine into large matrix

  france_matrices <- france_matrices[!is.na(france_matrices)]


  france_frame <- rbindlist(lapply(france_matrices, 
                                   setDT, 
                                   keep.rownames = TRUE), 
                            fill = TRUE)

  # Ensure NA's are set to zero, as that is what they represent.
  france_frame[is.na(france_frame)] <- 0

  # Create a colour vector so we can colour code the plot
  france_frame$colours <- NA

  for (i in 1:nrow(france_frame)){
    if (france_frame$category[i] == "back"){
      france_frame$colours[i] <- "grey"
    }
    else if (france_frame$category[i] == "cumul"){
      france_frame$colours[i] <- "blue"
    }
    else if (france_frame$category[i] == "instant"){
      france_frame$colours[i] <- "red"
    }
    else {
      france_frame$colours[i] <- "gold"
    }
  }

  return(france_frame)
}

france_frame <- nMDS_frame_function(2080016510)

# Perform nMDS ordination.

NMDS=metaMDS(france_frame[,-c("rn", "category", "colours")], # Our community-by-species matrix
                     k=2, trymax = 1000)

# Plot it! I opted to use custom plotting here instead of
# the ordiplot function as it easier to colour code points
# this way.

plot(NMDS, type = "n")

points(x = NMDS$points[,"MDS1"], 
       y = NMDS$points[, "MDS2"], 
       col = france_frame$colours,
       pch = 16,
       cex = 1.1)

text(x = NMDS$species[, "MDS1"],
     y = NMDS$species[, "MDS2"], 
     labels= 1:nrow(NMDS$species), 
     cex=0.5)

legend("topleft", legend=c("Back", "Cumul", "Instant", "Novel", "Species"),
             col=c("grey", "blue", "red", "gold", "black"),pch = c(16,16,16,16,3), cex=0.8)

# Check stress if wanted.
stressplot(NMDS)









