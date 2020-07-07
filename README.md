# SealScent2020
Code used in analyses in "Chemical patterns of colony membership and mother-offspring similarity in Antarctic fur seals are reproducible".

## .RData Objects

'seal_raw_dfs.RData' - 
Contains a list with raw data frames of peak identification data for Antarctic fur seals sampled in 2016/2017 on Bird Island, South Georgia

'mom_pup_alignment_GCalignR.RData' - 
Contains an object for aligned chromatograms for scent data of Antarctic fur seal mother pup pairs from two colonies

'mom_pup_nmds_scaling.RData' - 
Contains multiple objects that were created in R in a non-metric multidimensional procedure on 'mom_pup_alignment_GCalignR.RData'

'pup_colonies_alignment_GCalignR.RData' - 
Contains an object for aligned chromatograms for scent data of Antarctic fur seal mother pups from six colonies

'pup_colonies_nmds_scaling.RData' - 
Contains multiple objects that were created in R in a non-metric multidimensional procedure on 'pup_colonies_alignment_GCalignR.RData'

'R2_initial_season_btrap.RData' - 
Holds a data frame for bootstrapped data of Antarctic fur seal mother pup-pair scent data sampled in 2010/2011 in two breeding colonies

'R2_replication_season_btrap' - 
Holds a data frame for bootstrapped data of Antarctic fur seal mother pup-pair scent data sampled in 2016/2017 in two breeding colonies

## .R scripts
'ReadFiles.R' - 
Script to read in a directory that contains individual reports about peak integration in OpenChrom Style. Reads in text file and cuts off 
everything except peak retention times

'align_furseal_peaks.R' - 
Aligns chromatograms. For that chromatograms are ideally saved in a list containing data frames with peak individual per individual.
List strucutre: chromatrograms.list <- individual.information <- peak_retentiontimes.df

'effect_size_estimate_plot.R' - 
Code to read in bootstrap data that ran on a server and subsequently plot meaningful values for effect size estimates.
