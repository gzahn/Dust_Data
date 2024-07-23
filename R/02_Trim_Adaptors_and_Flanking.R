# SETUP ####

## Packages ####
library(tidyverse); packageVersion('tidyverse')
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel); packageVersion("parallel")

## Functions ####
source("./R/functions.R")

## Data ####
metadata <- readRDS("./data/full_clean_metadata.RDS") %>% 
  dplyr::filter(!is.na(fwd_filepath))

# if running on local machine  
metadata$fwd_filepath <- paste0("./data/raw/",basename(metadata$fwd_filepath))
metadata$rev_filepath <- paste0("./data/raw/",basename(metadata$rev_filepath))


## primer sequences ####

# SSU
WANDAf <- "CAGCCGCGGTAATTCCAGCT"  
AML2r <- "GAACCCAAACACTTTGGTTTCC"  
# ITS1
ITS1f <- "CTTGGTCATTTAGAGGAAGTAA" 
ITS2r <- "GCTGCGTTCTTCATCGATGC" 


# RUN CUTADAPT ####

## On SSU Samples ####
remove_primers(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
               amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
               amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
               sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
               fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to raw data
               rev.fp.colname = "rev_filepath",
               fwd_pattern="_R1_",
               rev_pattern="_R2_",
               fwd_primer=WANDAf,
               rev_primer=AML2r,
               multithread=parallel::detectCores()-1)

## On ITS Samples ####
remove_primers(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
               amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
               amplicon = "ITS", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
               sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
               fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to raw data
               rev.fp.colname = "rev_filepath",
               fwd_pattern="_R1_",
               rev_pattern="_R2_",
               fwd_primer=ITS1f,
               rev_primer=ITS2r,
               multithread=parallel::detectCores()-1)


# Run itsxpress
  #on just ITS samples (and just FWD ITS1 samples, actually)




# testing file output on local computer
cutadapt_ftest_fwd <- paste0("./data/raw/cutadapt/",metadata$index,"_ITS_cutadapt_fwd.fastq.gz")
cutadapt_ftest_rev <- paste0("./data/raw/cutadapt/",metadata$index,"_ITS_cutadapt_rev.fastq.gz")

missing_fwds <- cutadapt_ftest_fwd[!file_test("-f",cutadapt_ftest_fwd)]
missing_revs <- cutadapt_ftest_rev[!file_test("-f",cutadapt_ftest_fwd)]

metadata$local_fwd_cutadapt_paths <- cutadapt_ftest_fwd
metadata$local_rev_cutadapt_paths <- cutadapt_ftest_rev

saveRDS(metadata,"./data/cutadapt_metadata.RDS")
