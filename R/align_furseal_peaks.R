require(stringr)
require(GCalignR)
require(vegan)
require(dplyr)
source("R/ReadFiles.R")

all_dfs <- suppressWarnings(ReadFiles(path = "data"))

for (i in 1:length(names(all_dfs))) {
  all_dfs[[i]] <- all_dfs[[i]] %>% filter(RT >= 15)
} 

# now check if any zero entries are in data tibbles
x <- NULL
for (i in 1:length(names(all_dfs))) {
  hold <- all_dfs[[i]]
  x[i] <- length(rownames(hold[,1]))
}
# delete zero entries
all_dfs <- all_dfs[-which(x == 0)]

seal_dfs.list <- all_dfs[c(20:69,90:199)]
save(seal_dfs.list, file = "RData/objects/seal_raw_dfs.Rdata")

# index correct MP & pup_colony subset
index_MP <- c(1:6, 20:69, 90:100, 111:150, 200:201)
index_pupcols <- c(1:6, 90:201)

# get blank names (same for index_MP and index_pupcols)
blank_samples <- names(all_dfs[index_MP])[stringr::str_detect(names(all_dfs[index_MP]), "DCM")] %>%
  c(., names(all_dfs[index_MP])[stringr::str_detect(names(all_dfs[index_MP]), "Syr")])

# Aligning mom-pup colonies
mom_pup_aligned <- align_chromatograms(data = all_dfs[index_MP], # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         rt_cutoff_low = 15, # remove peaks below 8 Minutes
                                         rt_cutoff_high = 54.7, # remove peaks exceeding 54.7 Minutes
                                         reference = "P13", # sample with highest overall peak number
                                         max_linear_shift = 0.05, # Premise 1: max. shift for linear corrections
                                         max_diff_peak2mean = 0.08, # Premise 2: max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.03, # Premise 3: min. expected distance between peaks 
                                         blanks = blank_samples, # negative control(s)
                                         delete_single_peak = T, # delete peaks that are present in just one sample 
                                         write_output = NULL) # add variable names to write aligned data to text files
save(mom_pup_aligned, file = "RData/objects/mom_pup_alignment_GCalignR.RData")

pup_colonies_aligned <- align_chromatograms(data = all_dfs[index_pupcols], # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         rt_cutoff_low = 15, # remove peaks below 8 Minutes
                                         rt_cutoff_high = 54.7, # remove peaks exceeding 54.7 Minutes
                                         reference = "P13",# sample with highest overall peak number 
                                         max_linear_shift = 0.05, # Premise 1: max. shift for linear corrections
                                         max_diff_peak2mean = 0.08, # Premise 2: max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.03, # Premise 3: min. expected distance between peaks 
                                         blanks = blank_samples, # negative control
                                         delete_single_peak = F, # delete peaks that are present in just one sample 
                                         write_output = NULL) # add variable names to write aligned data to text files
save(pup_colonies_aligned, file = "RData/objects/pup_colonies_alignment_GCalignR.RData")

## Preliminary NMDS Visualization
scent_factors_raw <- readr::read_delim("documents/metadata_seal_scent.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
scent_factors_raw <- as.data.frame(scent_factors_raw[-c(194:209),])

# set sample names as row names, ensure there are no duplicates
scent_factors <- scent_factors_raw[,-1]
rownames(scent_factors) <- scent_factors_raw[,1]

## check for empty samples, i.e. no peaks
x <- apply(peak_data_aligned$aligned$RT, 2, sum)
x <- which(x == 0)

## normalise area and return a data frame
scent <- norm_peaks(peak_data_aligned, conc_col_name = "Area",rt_col_name = "RT",out = "data.frame") 
## common transformation for abundance data to reduce the extent of mean-variance trends
scent <- log(scent + 1) 
# scent <- scent[-11,]

## subset scent_factors
scent_factors <- scent_factors[rownames(scent_factors) %in% rownames(scent),]
scent <- scent[rownames(scent) %in% rownames(scent_factors),]

## keep order of rows consistent
scent <- scent[match(rownames(scent_factors),rownames(scent)),] 

## get number of compounds for each individual sample after alignment
num_comp <- as.vector(apply(scent, 1, function(x) length(x[x>0]))) 

## bray-curtis similarity
scent_nmds.obj <- vegan::metaMDS(comm = scent, k = 2, try = 999, trymax = 3000, distance = "bray") #[166:184,]
# scent_nmds <- with(scent_factors, MDSrotate(scent_nmds, family))

## get x and y coordinates

scent_nmds <- as.data.frame(scent_nmds.obj[["points"]])  
## add the colony as a factor to each sample
scent_nmds <- cbind(scent_nmds,
                    age = scent_factors[["age"]],
                    tissue_tag = scent_factors[["tissue_tag"]],
                    colony = scent_factors[["colony"]],
                    family = as.factor(scent_factors[["family"]]),
                    clr = as.factor(scent_factors[["clr"]]),
                    shp = as.factor(scent_factors[["shp"]]),
                    gcms = as.factor(scent_factors[["gcms_run"]]),
                    peak_res = as.factor(scent_factors[["peak_res"]]),
                    sample_qlty = as.factor(scent_factors[["sample_qlty"]]),
                    vialdate = as.factor(scent_factors[["gcms_vialdate"]]),
                    captured = as.factor(scent_factors[["capture_date"]]),
                    sex = scent_factors[["sex"]],
                    num_comp = num_comp)
scent_nmds <- scent_nmds %>% mutate(BeachAge = str_c(colony, age, sep = "_"))

ggplot(data = scent_nmds) + 
  geom_point(size = 4.5, aes(MDS1, MDS2, color = BeachAge, shape = BeachAge)) + 
  scale_shape_manual(values = c(19, 1, 19, 1), 
                     labels = c("FWB mothers ", "FWB pups ", "SSB mothers ", "SSB pups ")) +
  scale_color_manual(values = c("#D55E00", "#D55E00", "#56B4E9", "#56B4E9"), 
                     labels = c("FWB mothers ", "FWB pups ", "SSB mothers ", "SSB pups ")) +
  theme_void() + 
  # ylim(-0.75,1.1) +
  # ggtitle("Colony Differences SSB & FWB '16/'17") +
  # annotate("text", x = 0.64, y = 1.1, label = "(A)", size = 5) +
  # annotate("text", x = 0.48, y = -0.75, label = "2D Stress: 0.23", size = 5) +
  theme(panel.background = element_rect(colour = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(size = 0.3, linetype = "solid", color = "black")) 
