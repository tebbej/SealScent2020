ReadFiles <- function(path = NULL) {

library(readr)
library(dplyr)
library(purrr)
# test <- read_delim(file = "Desktop/DCM7_26092018_3.txt", delim = "/t")

# Get all the filenames from your working directory and place it in all_files as string/character entries. In this case read in all names of
# an OpenChrom batch report file directory.
all_files <- list.files(path = path)
sample_retries <- all_files[grepl("*retry*", all_files)]
all_files <- all_files[-(which(all_files %in% sample_retries))]

names_samples <- lapply(all_files, function(x) {
  if (stringr::str_detect(x, "Env")) {
    out <- stringr::str_split(x, "_")[[1]][1:2] %>% 
      paste0(., collapse = "_")
  # } else if (stringr::str_detect(x, "neu")) {
  #   out <- stringr::str_split(x, "_")[[1]][1:2] %>% 
  #     paste0(., collapse = "_")
  # } else if (stringr::str_detect(x, "resolved")) {
  #   out <- stringr::str_split(x, "_")[[1]][1:2] %>% 
  #     paste0(., collapse = "_")
  }  else {
    out <- stringr::str_split(x, "_")[[1]][1] %>% 
      paste0(., collapse = "_")
  }
}) %>% unlist()  


# Read in a file patch to perform read_chroma_files.
read_chroma_files <- function(file_path) {
    # read in the file as one object
  file_full <- read_lines(file = file_path)
    # peak quant line is now able to give the exact position when the Peak Quantitation Summary Information is reached in the file. 
    # (Varies for different lengths of chromatograms!)
  peak_quant_line <- which(file_full == "PEAK QUANTITATION SUMMARY") + 1
  
    # number_of_records estimates the length of the records in Peak Quantitation Summary
  number_of_records <- length(rownames(file_full)) - peak_quant_line
  
  options(digits=15)
    # Using readr we can read in a tab delimited file that skips all file entries to the start of the Peak Quantitation Summary
  sample_quant <- read_delim(file_path, 
                           "\t", escape_double = FALSE, trim_ws = TRUE, comment = "~",
                           skip = peak_quant_line, col_types = "cc") %>%
     # All OpenChrom Batch Report Files end with 2 extra lines that must be cut to get a correct tibble.
    #.[-nrow(.), ] %>%
    rename(RT = "RT (Minutes)") %>%
      # For some reason, read_delim reads data
      # as characters and not as numeric variables in this case.
    mutate_all(as.numeric)
  
  
} #end read_chroma_files

  # with library(purrr) we can perform our function read_chroma_files on all files that have been loaded in with all_files 
  # (Chromatogram files of interest in your working directory)

all_files <- file.path(path, all_files)
all_dfs <- map(all_files, read_chroma_files)

  # rename tibbles in the list to the corresponding sample names
names(all_dfs) <- names_samples
return(all_dfs)
}