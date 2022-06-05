# Function that computes the lengths of novel states in our data,
# returns a useful list with metrics.

novel.length.calculator <- function(anosim.plots.list, cut.off){

# Filter out novel communities that show a lack of difference. It also 
# returns the data in a useful format. 

true.novelty.output <- cut.off.generator(anosim.plots.list)

# Calculate lengths of novel communities, if they exceed the threshold value. Returns length and R stat. 

full.novelty.lengths <- novel.length.checker(true.novelty.output, cut.off)

return(full.novelty.lengths)

}
