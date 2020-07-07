load("RData/objects/mom_pup_alignment_GCalignR.RData")

scent_factors_raw <- read_delim("documents/metadata_seal_scent.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
scent_factors_raw <- as.data.frame(scent_factors_raw[-c(194:209),])

# set sample names as row names, ensure there are no duplicates
scent_factors <- scent_factors_raw[,-1]
rownames(scent_factors) <- scent_factors_raw[,1]

## check for empty samples, i.e. no peaks
x <- apply(mom_pup_aligned$aligned$RT, 2, sum)
x <- which(x == 0)

## normalise area and return a data frame
scent <- norm_peaks(mom_pup_aligned, conc_col_name = "Area",rt_col_name = "RT",
                    out = "data.frame") 
## common transformation for abundance data to reduce the extent of mean-variance trends
scent <- log(scent + 1) 

## subset scent_factors
scent_factors <- scent_factors[rownames(scent_factors) %in% rownames(scent),]
scent <- scent[rownames(scent) %in% rownames(scent_factors),]

## keep order of rows consistent
scent <- scent[match(rownames(scent_factors),rownames(scent)),] 

## get number of compounds for each individual sample after alignment
num_comp <- as.vector(apply(scent, 1, function(x) length(x[x>0]))) 

## bray-curtis similarity
scent_nmds.obj <- vegan::metaMDS(comm = scent, k = 2, try = 999, 
                                 trymax = 9999, distance = "bray") 

scent_nmds <- as.data.frame(scent_nmds.obj[["points"]])  

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
# creates & adds new variable BeachAge and simplifies plotting