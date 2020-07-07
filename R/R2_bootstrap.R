# Bootstrap to track R2 values for randomized subsets. In addition,
# bootstrap cannot only be used to randomize the chemical data frame, 
# but also to evaluate R2 change for different subsets based on different
# premises. 1) Frequent peaks 2) Strong concentrations 3) Peaks identified by SIMPER

require(vegan)

# path: file path to scent_nmds-mompup2017_ssbfwb.RData", 
#objects: scent_nmds, scent_nmds.obj, scent_factors, scent
# df.permutations: number of times the scent.df from loaded data will be permuted
# nmds.permutations: number of permutation in nMDS using Bray-Curtis
# btrap.iterations: number of procedure repeats

scent_btrap_r2_swarm_data <- function(path, df.permutations = 15, 
                                      nmds.permutations = 999, 
                                      btrap.iterations = 5000){
  # Create a data frame by permuting the data for scent
  # compounds data and also ensure that each population*age occur
  # same amounts of time in the permutation data frame. 
  #---------------------------------------
  
  # load data frame with data of aligned fur seal chromatograms
  load(path)
  scent_factors <- peak_factors
  # transfer BeachAge Column from scent_nmds to meta data.frame scent_factors
  scent_factors <- cbind(scent_factors, 
                         BeachAge = scent_nmds$BeachAge)
  
  # create index column for meta data frame
  scent_factors <- cbind(scent_factors,
                         SampleIndex = 1:length(rownames(scent_factors)))
  
  
  # create data.frame to track PERMANOVA results over repeated tests
  nonsubset_results_paov <- data.frame(R2_age = double(), p_colfam = double(),
                                       R2_residual = double(), 
                                       F_Het = double(), p_Het = double())
  promcomp_results_paov <- data.frame(R2_age = double(), p_colfam = double(),
                                      R2_residual = double(), 
                                      F_Het = double(), p_Het = double())
  highcomp_results_paov <- data.frame(R2_age = double(), p_colfam = double(),
                                      R2_residual = double(), 
                                      F_Het = double(), p_Het = double())
  simper_results_paov <- data.frame(R2_age = double(), p_colfam = double(),
                                    R2_residual = double(), 
                                    F_Het = double(), p_Het = double())
  
  # create list to store created objects in an iteration
  iter_object_container <- list()
  
  for (i in 1:btrap.iterations) {
    
    
    # create data.frame subsets (colony subset) by indexing the meta data.frame 
    scent.f.ssb.m <- scent_factors[scent_factors$BeachAge == "SSB_1",]
    scent.f.fwb.m <- scent_factors[scent_factors$BeachAge == "FWB_1",]
    scent.f.ssb.p <- scent_factors[scent_factors$BeachAge == "SSB_2",]
    scent.f.fwb.p <- scent_factors[scent_factors$BeachAge == "FWB_2",]
    
    # int vector of row index number of permuted scent.ssb data.frame
    # row numbers will be used to create a permuted data.frame of 
    # evenly distributed draws of individuals
    permute_rows_ssb_m <- sample(scent.f.ssb.m$SampleIndex, df.permutations, replace = T)
    permute_rows_fwb_m <- sample(scent.f.fwb.m$SampleIndex, df.permutations, replace = T)
    permute_rows_ssb_p <- sample(scent.f.ssb.p$SampleIndex, df.permutations, replace = T)
    permute_rows_fwb_p <- sample(scent.f.fwb.p$SampleIndex, df.permutations, replace = T)
    
    # create overall index number that can be used to 
    #index data.frame(scent): index corresponds to correct individual 
    perm_index_all <- c(permute_rows_ssb_m, 
                        permute_rows_fwb_m,
                        permute_rows_ssb_p,
                        permute_rows_fwb_p)
    
    # create new data.frame with indeces found in permutation results vector perm_index_all
    scent.permute <- scent[perm_index_all,]
    scent_factors.permute <- scent_factors[perm_index_all,]
    # rownames(scent.permute) == rownames(scent_factors.permute) # TRUE
    #---------------------------------------
    
    # Perform analysis to find 3 subsets based on different premises with the permuted data frame. 
    # Track 15 best performing compounds of an analysis
    #---------------------------------------
    
    ## NDMS scale results
    ## count number of peaks that are not 0 per column
    peak_count <- as.vector(apply(scent.permute, 2, function(x) length(x[x>0]))) 
    
    ## add peaks in a column that are not 0 to estimate highest concentration peak sum
    peak_add <- as.vector(apply(scent.permute, 2, function(x) sum(x))) #[x>0]
    
    ## create dataframe with same name properties as scent.RData
    compound_subset <- data.frame(name = colnames(scent.permute) , peak_count, peak_add)
    
    ## sort data frame for most prominent compounds over all samples
    most_abundant <- compound_subset %>% arrange(desc(peak_count))
    
    ## shorten scent matrix to only the 15 most abundant compounds
    scent.promcomp <- scent.permute[colnames(scent.permute) %in% most_abundant$name[1:15]]
    
    ## sort data frame for most highly concentrated compounds over all samples
    most_concentration <- compound_subset %>% arrange(desc(peak_add))
    
    ## shorten scent matrix to only the 15 most abundant compounds
    scent.highcomp <- scent.permute[colnames(scent.permute) %in% most_concentration$name[1:15]]
    
    ## simper
    # simper analysis and results array
    sim <- with(scent_factors.permute, 
                simper(scent.permute, colony))
    best.compounds.simper.btrap <- summary(sim)[[1]]
    #filter 15 compounds that contribute most towards dissimilarity of individuals
    simper_comps <- as.numeric(rownames(best.compounds.simper.btrap))
    best_comps <- simper_comps[1:15] 
    # subset peak data matrix {scent}
    scent.simper.btrap <- scent.permute[,which(colnames(scent.permute) %in% as.character(best_comps))]
    #---------------------------------------
    
    # Take 15 identified compounds and limit nMDS of the permuted data frame (scent.permute) to only those compounds
    #---------------------------------------
    
    # bray-curtis similarity
    scent_nmds_regular.obj <- vegan::metaMDS(comm = scent.permute, k = 2, 
                                             try = df.permutations, distance = "bray")
    scent_nmds_count.obj <- vegan::metaMDS(comm = scent.promcomp, k = 2, 
                                           try = df.permutations, distance = "bray") 
    scent_nmds_add.obj <- vegan::metaMDS(comm = scent.highcomp, k = 2, 
                                         try = df.permutations, distance = "bray") 
    scent_nmds_simper.obj <- vegan::metaMDS(comm = scent.simper.btrap, k = 2, 
                                            try = df.permutations, distance = "bray") #[166:184,]
    
    ## get x and y coordinates
    scent_nmds_regular <- as.data.frame(scent_nmds_regular.obj[["points"]])
    scent_nmds_count <- as.data.frame(scent_nmds_count.obj[["points"]]) 
    scent_nmds_add <- as.data.frame(scent_nmds_add.obj[["points"]])
    scent_nmds_simper <- as.data.frame(scent_nmds_simper.obj[["points"]])
    
    ## add the colony as a factor to each sample
    scent_nmds <- data.frame(MDS1r = scent_nmds_regular[["MDS1"]],
                             MDS2r = scent_nmds_regular[["MDS2"]],
                             MDS1c = scent_nmds_count[["MDS1"]],
                             MDS2c = scent_nmds_count[["MDS2"]],
                             MDS1a = scent_nmds_add[["MDS1"]],
                             MDS2a = scent_nmds_add[["MDS2"]],
                             MDS1s = scent_nmds_simper[["MDS1"]],
                             MDS2s = scent_nmds_simper[["MDS2"]], 
                             age = scent_factors.permute[["age"]],
                             colony = scent_factors.permute[["colony"]],
                             family = scent_factors.permute[["family"]],
                             BeachAge = scent_factors.permute[["BeachAge"]]
    )
    #---------------------------------------
    
    # Perform PERMANOVA on distance matrix based limited scent compounds data
    #---------------------------------------
    
    # not subsetted
    nonsubset.df_permanova <- adonis(scent.permute ~ age + colony + colony:family,
                                     data = scent_factors.permute,
                                     permutations = 9999)
    nonsubset.df_hetgeneity <- anova(betadisper(vegdist(scent.permute), 
                                                scent_factors.permute$colony))
    
    # track important values of statistical analysis in this run
    nonsubset_iter_res_paov <- cbind(R2_age = nonsubset.df_permanova$aov.tab$R2[1],
                                     R2_colony = nonsubset.df_permanova$aov.tab$R2[2],
                                     R2_famcol = nonsubset.df_permanova$aov.tab$R2[3],
                                     R2_residual = nonsubset.df_permanova$aov.tab$R2[4],
                                     F_Het = nonsubset.df_hetgeneity$`F value`[1],
                                     p_Het = nonsubset.df_hetgeneity$`Pr(>F)`[1])
    
    # bind run values to track changes over iterations in the for-loop
    nonsubset_results_paov <- rbind(nonsubset_results_paov,
                                    nonsubset_iter_res_paov)
  
    #prom comps
    promcomp.df_permanova <- adonis(scent.promcomp ~ age + colony + colony:family, 
                                    data = scent_factors.permute, 
                                    permutations = 9999) 
    promcomp.df_hetgeneity <- anova(betadisper(vegdist(scent.promcomp), 
                                               scent_factors.permute$colony)) 
    
    promcomp_iter_res_paov <- cbind(R2_age = promcomp.df_permanova$aov.tab$R2[1],
                                    R2_colony = promcomp.df_permanova$aov.tab$R2[2],
                                    R2_famcol = promcomp.df_permanova$aov.tab$R2[3],
                                    R2_residual = promcomp.df_permanova$aov.tab$R2[4],
                                    F_Het = promcomp.df_hetgeneity$`F value`[1],
                                    p_Het = promcomp.df_hetgeneity$`Pr(>F)`[1])
    
    promcomp_results_paov <- rbind(promcomp_results_paov,
                                   promcomp_iter_res_paov)
    
    # high comps
    highcomp.df_permanova <- adonis(scent.highcomp ~ age + colony + colony:family, 
                                    data = scent_factors.permute, 
                                    permutations = 9999) 
    highcomp.df_hetgeneity <- anova(betadisper(vegdist(scent.highcomp), scent_factors.permute$colony)) 
    
    highcomp_iter_res_paov <- cbind(R2_age = highcomp.df_permanova$aov.tab$R2[1],
                                    R2_colony = highcomp.df_permanova$aov.tab$R2[2],
                                    R2_famcol = highcomp.df_permanova$aov.tab$R2[3],
                                    R2_residual = highcomp.df_permanova$aov.tab$R2[4],
                                    F_Het = highcomp.df_hetgeneity$`F value`[1],
                                    p_Het = highcomp.df_hetgeneity$`Pr(>F)`[1])
    
    highcomp_results_paov <- rbind(highcomp_results_paov,
                                   highcomp_iter_res_paov)
    # SIMPER
    simper.df_permanova <- adonis(scent.simper.btrap ~ age + colony + colony:family, 
                                  data = scent_factors.permute, 
                                  permutations = 9999) 
    simper.df_hetgeneity <- anova(betadisper(vegdist(scent.simper.btrap), scent_factors.permute$colony)) 
    
    simper_iter_res_paov <- cbind(R2_age = simper.df_permanova$aov.tab$R2[1],
                                  R2_colony = simper.df_permanova$aov.tab$R2[2],
                                  R2_famcol = simper.df_permanova$aov.tab$R2[3],
                                  R2_residual = simper.df_permanova$aov.tab$R2[4],
                                  F_Het = simper.df_hetgeneity$`F value`[1],
                                  p_Het = simper.df_hetgeneity$`Pr(>F)`[1])
    
    simper_results_paov <- rbind(simper_results_paov,
                                 simper_iter_res_paov)
    #---------------------------------------
    
    
    # # pack all this in a list to be later on stored in a list that can be saved again
    # create name giving the iteration step
    iteration_count <- paste0("iter_", i)
    
    # create list that stores relevant workspace elements for an iteration step
    iter_objects <- list(scent.permute = scent.permute,
                         scent_factors.permute = scent_factors.permute,
                         scent.promcomp = scent.promcomp,
                         scent.highcomp = scent.highcomp,
                         sim = sim,
                         scent.simper.btrap = scent.simper.btrap,
                         scent_nmds_regular.obj = scent_nmds_regular.obj,
                         scent_nmds_count.obj = scent_nmds_count.obj,
                         scent_nmds_add.obj = scent_nmds_add.obj,
                         scent_nmds_simper.obj = scent_nmds_simper.obj,
                         scent_nmds_regular = scent_nmds_regular,
                         scent_nmds_count = scent_nmds_count,
                         scent_nmds_add = scent_nmds_add,
                         scent_nmds_simper = scent_nmds_simper,
                         promcomp.df_permanova = promcomp.df_permanova,
                         promcomp.df_hetgeneity = promcomp.df_hetgeneity,
                         highcomp.df_permanova = highcomp.df_permanova,
                         highcomp.df_hetgeneity = highcomp.df_permanova,
                         simper.df_permanova = simper.df_permanova,
                         simper.df_hetgeneity = simper.df_hetgeneity)
    
    # save everything as a list in a container list, that stores 
    # information/elements of all iteration steps
    iter_object_container[[i]] <- iter_objects
    names(iter_object_container)[i] <- iteration_count
    
  } # end i
  
  paov_r2_results <- list(regular = nonsubset_results_paov,
                          promcomp = promcomp_results_paov,
                          highcomp = highcomp_results_paov,
                          simper_res = simper_results_paov)
  return(list(paov_r2_results = paov_r2_results, 
              iter_object_container = iter_object_container))
} # end function "scent_btrap_r2_swarm_data.R"
