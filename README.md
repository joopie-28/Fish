# Novel Communities of Freshwater Fishes

This repository contains the R-code used to perform all main analyses for the "Novel Freshwater Fish Communities" project.

Work on this started 29-11-2021 and is still ongoing at the time of writing. In addition to main analyses pathways, the repository also includes relevant 
CSV output files and plots. It also contains a lot of 'scraps'; files which are no longer needed or include ideas that will not be further pursued in this
work.

Anyone who wishes to recreate the analyses done here should read the following:

There are two main scripts which are needed to execute the analyses:

### 1. complete_analysis_abundance.R
This script applies the novelty detection framework (Pandolfi et al.) to the 
rivfishtime database (note that the rivfishtime files are present in this repo).
The final output of running this script will be a large list containing all the 
relevant novelty output for each timeseries, at different bin sizes (1 or 2 years).
A copy of this script that does the same but first constructs presence/absence matrices
from the raw fish data is also available. However, we opted to use the relative abundance
pathway for further analyses.

### 2. novel_analyses.R
This script features all the analyses beyond simple novelty detection. It includes tagging
fish species as invaders or natives based on literature, determining lengths of novel
states using ANOSIM and also analyses how many novel communities are 'blips'. 

In addition to these scripts, there are a number of important files that are necessary to perform the
analyses, such as species status data, a host of smaller functions found in the 'functions' folder, to 
name a few. Any user should thus aim to run these scripts within an R-environment that includes all 
relevant data, or particular parts of the code will not execute.

On a final note, script layout is quite modular, which I feel improves readbility. However, this means that if you want to examine what a function is doing 'under the hood', you need to go the specific script where the function is declared (these scripts generally have the same name as the function so isn't too problematic).


