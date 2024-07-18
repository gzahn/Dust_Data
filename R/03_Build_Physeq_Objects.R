# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(dada2)
library(decontam)

# functions
source("./R/functions.R")

# metadata
metadata <- readRDS("./data/full_clean_metadata.RDS")

# list of sequencing runs
all_runs <- metadata$seq_run %>% as.character() %>% unique

# # EXAMPLE ####
# # Run 1
# build_asv_table(metadata=metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
#                 run.id.colname = "seq_run", # name of column in metadata indicating which sequencing run a sample comes from
#                 run.id = "1", # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
#                 amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
#                 amplicon = "ITS", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
#                 sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
#                 fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
#                 rev.fp.colname = "rev_filepath", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
#                 fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
#                 rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
#                 maxEE = c(2,2), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
#                 truncQ = 2, # special value denoting "end of good quality sequence"
#                 rm.phix = TRUE, # remove phiX sequences?
#                 compress = TRUE, # gzip compression of output?
#                 multithread = (parallel::detectCores() -1), # how many cores to use? Set to FALSE on windows
#                 single.end = TRUE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
#                 filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
#                 asv.table.dir = "testing", # path to directory where final ASV table will be saved
#                 random.seed = 666
#                 )

# RUN ALL ITS DATA ####
for(seqrun in all_runs){
  build_asv_table(metadata=metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                  run.id.colname = "seq_run", # name of column in metadata indicating which sequencing run a sample comes from
                  run.id = i, # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
                  amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                  amplicon = "ITS", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                  sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
                  fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
                  rev.fp.colname = "rev_filepath", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
                  fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
                  rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
                  maxEE = c(2,2), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
                  truncQ = 2, # special value denoting "end of good quality sequence"
                  rm.phix = TRUE, # remove phiX sequences?
                  compress = TRUE, # gzip compression of output?
                  multithread = (parallel::detectCores() -1), # how many cores to use? Set to FALSE on windows
                  single.end = TRUE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
                  filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
                  asv.table.dir = "testing", # path to directory where final ASV table will be saved
                  random.seed = 666
  )
}

# RUN ON ALL SSU DATA ####
for(seqrun in all_runs){
  build_asv_table(metadata=metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                  run.id.colname = "seq_run", # name of column in metadata indicating which sequencing run a sample comes from
                  run.id = i, # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
                  amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                  amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                  sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
                  fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
                  rev.fp.colname = "rev_filepath", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
                  fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
                  rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
                  maxEE = c(2,2), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
                  truncQ = 2, # special value denoting "end of good quality sequence"
                  rm.phix = TRUE, # remove phiX sequences?
                  compress = TRUE, # gzip compression of output?
                  multithread = (parallel::detectCores() -1), # how many cores to use? Set to FALSE on windows
                  single.end = FALSE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
                  filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
                  asv.table.dir = "testing", # path to directory where final ASV table will be saved
                  random.seed = 666
  )
}
