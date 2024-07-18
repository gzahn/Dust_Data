# SETUP ####

## Packages ####
library(tidyverse)
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(parallel); packageVersion("parallel")

## Functions ####
source("./R/functions.R")

## Data ####
metadata <- readRDS("./data/full_clean_metadata.RDS")
  
#parse filenames

## primer sequences ####

# SSU
WANDAf <- "CAGCCGCGGTAATTCCAGCT"  
AML2r <- "GAACCCAAACACTTTGGTTTCC"  
# ITS1
ITS1f <- "CTTGGTCATTTAGAGGAAGTAA" 
ITS2r <- "GCTGCGTTCTTCATCGATGC" 

# RUN CUTADAPT ####
  
## On SSU Samples ####
remove_primers(metadata2, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
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
remove_primers(metadata2, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
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
  #on just ITS samples

