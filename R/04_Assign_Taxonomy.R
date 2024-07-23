# SETUP ####

## packages ####
library(tidyverse)
library(dada2)

## functions ####
source("./R/functions.R")

## file paths ####
  # asv tables

  # databases

  # output files



# ASSIGN TAXONOMY AND EXPORT ####

# assign taxonomy
assign_taxonomy_to_asv_table(asv.table, # asv table object name
                             tax.database, # path to taxonomic database (fasta)
                             multithread=(parallel::detectCores()-1), # set to FALSE on Windows
                             random.seed=666,
                             try.rc = TRUE, # attempt revComplement assignments as well? (doubles time)
                             min.boot=50)
# export as RDS


