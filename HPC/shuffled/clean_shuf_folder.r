# ----------clean_shuf_folder.r---------
# Used to clean files from a shuffling HPC folder:
# arguments: 1. -all or -clean
#            2. <directory_path_to_clean>
#
# How to ectivate:
# Rscript clean_shuf_folder.r [-clean\-all] ./shuffled_folder_name
# ----------------

# ------ include ----


# ------ consts ----
is_nuke_state <- FALSE
shuff_folder_name <- ""

# files to remove on CLEANING state
list_of_endings_to_remove <- c("\\.R",
                               "\\.sh",
                               "\\.rds")
list_of_contains_to_remove <- c("experiments\\.csv",
                                "\\.e",
                                "\\.o",
                                "Infomap",
                                "infomap_multilayer\\.txt")
list_of_starting_to_remove <- c("core_ASV_","shuff_farm_")

# files to keep in NUKING state
list_of_endings_to_keep <- c("_edge_list\\.csv",
                             "_log\\.txt") # in files
list_of_starting_to_keep <- c("shuff_farm_", 
                              "run_summary\\.csv",
                              "experiments\\.csv") # in main folder

# ------ argument handling -----
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop('No arguments were found, cleaning stopped.')
} else {
  is_nuke_state <- args[1] == "-all" # remove ALL but edge list and log.
  shuff_folder_name <- args[2] # use this folder for running
  if (!dir.exists(shuff_folder_name)) {
    stop("Directory doesn't exist, cleaning stopped.")
  }
}

# ---- run the cleaning -----
folder_names <- list.dirs(shuff_folder_name)[-1]

if  (is_nuke_state) { # Remove everything but the networks themselves
  # build regex to get all files to keep - nuke mode
  regex_ends <- paste(list_of_endings_to_keep, "$",sep='', collapse = '|')
  regex_begin <- paste("^",list_of_starting_to_keep, sep='', collapse = '|')
  
  # remove files from each directory
  for (inner_folder_name in folder_names) {
    keep <- list.files(path = inner_folder_name,
                       pattern = regex_ends, full.names = T)
    all <- list.files(path = inner_folder_name, full.names = T)
    remove <- all[!(all %in% keep)]
    unlink(remove)
    
  }
  
  # remove all but abundance files etc from the main folder
  keep_main <- list.files(path = shuff_folder_name,
                     pattern = regex_begin, full.names = T)
  all_main <- list.files(path = shuff_folder_name, full.names = T)
  remove_main <- all_main[!(all_main %in% keep_main)]
  unlink(remove_main)
  
} else{ # Only needed to clean - analysis results are to be kept 
  
  # build regex to get all files to remove - cleaning mode
  regex_begin <- paste("^",list_of_starting_to_remove, sep='', collapse = '|')
  regex_contains <- paste(list_of_contains_to_remove, sep='', collapse = '|')
  regex_ends <- paste(list_of_endings_to_remove, "$",sep='', collapse = '|')
  all_regex <- paste(regex_begin, regex_contains, regex_ends, sep='|')
  
  for (inner_folder_name in folder_names) {
    # get list to remove
    files_to_remove <- list.files(path = inner_folder_name,
                                  pattern = all_regex, full.names = T)
    
    unlink(files_to_remove)
  }
}
