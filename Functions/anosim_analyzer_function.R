# This function applies the ANOSIM algorithm to all novel comm.
# and returns a list of summary statistics which help us to
# objectively classify 'blips'. It takes the novelty framework
# output as an input. It is essentially an amalgamation of
# several different functions, which makes it somewhat complex.

# ANOSIM plots are created whilst the function is run, but these are 
# saved.

anosim.analyzer <- function(Fish_Communities_A){

rel.abu.1.com <- rbindlist(lapply(list(Fish_Communities_A$BioRealm_Novelty_A$palearctic_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$nearctic_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$afrotropics_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$neotropics_novelty_A,
                                       Fish_Communities_A$BioRealm_Novelty_A$australasia_novelty_A), FUN = "rbindlist"))

# List of sites in our community list
novel.list <- as.list(rel.abu.1.com[rel.abu.1.com$cat == "novel"]$site)


# Applies novelty.trajectory.plotter to a list of Time Series
anosim.plots.list <- novelty.trajectory.lists(novel.list)

return(anosim.plots.list)
}


