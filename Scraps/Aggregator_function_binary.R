### Building matrices based on varying bin widths as opposed to just taking time as is ###

# We are going to test different bin widths: 2 and 4 years.

# Our initial analysis will be done for binary, I think it makes sense to use some code that creates a complete binary
# dataframe, before any cut-offs are applied. 


matrix_list_B_bin2 <- list_matrix_B_bins_function(ID_list, 2)

# Calculate novelty and return output in a list

novelty_list_B_bin2 <- list_novelty_function(matrix_list_B_bin2)

novelty_analysis_output_B_bin2 <- novel.probability(novelty_list_B_bin2)



venn_plot_function(novelty_analysis_output_B_bin2)
 