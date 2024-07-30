# SETUP ####

## packages ####
library(tidyverse)
library(dada2)

## functions ####
source("./R/functions.R")


# get taxonomy database files
ITS_DB <- "./taxonomy/Eukaryome_General_ITS_v1.8_reformatted_maarjam.fasta.gz"
SSU_DB <- "./taxonomy/Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta.gz"

# get possible asv table files
# turn this into list.files() for final version
run_1_its <- "./data/ASV_Tables/Run_1_ITS_ASV_Table.RDS"
run_2_its <- "./data/ASV_Tables/Run_2_ITS_ASV_Table.RDS"
run_3_its <- "./data/ASV_Tables/Run_3_ITS_ASV_Table.RDS"
run_4_its <- "./data/ASV_Tables/Run_4_ITS_ASV_Table.RDS"
run_5_its <- "./data/ASV_Tables/Run_5_ITS_ASV_Table.RDS"
run_6_its <- "./data/ASV_Tables/Run_6_ITS_ASV_Table.RDS"
run_7_its <- "./data/ASV_Tables/Run_7_ITS_ASV_Table.RDS"

run_1_ssu <- "./data/ASV_Tables/Run_1_SSU_ASV_Table.RDS"
run_2_ssu <- "./data/ASV_Tables/Run_2_SSU_ASV_Table.RDS"
run_3_ssu <- "./data/ASV_Tables/Run_3_SSU_ASV_Table.RDS"
run_4_ssu <- "./data/ASV_Tables/Run_4_SSU_ASV_Table.RDS"
run_5_ssu <- "./data/ASV_Tables/Run_5_SSU_ASV_Table.RDS"
run_6_ssu <- "./data/ASV_Tables/Run_6_SSU_ASV_Table.RDS"
run_7_ssu <- "./data/ASV_Tables/Run_7_SSU_ASV_Table.RDS"



for(run in ls(pattern="run_._")){
  x <- get(run)
  if(file.exists(x)){
    asv <- readRDS(x)
    outfile <- str_replace(x,"_ASV","_Taxonomy")
    
    tax <- assign_taxonomy_to_asv_table(asv.table=asv,
                                        tax.database=ifelse(grepl("_SSU_ASV_Table",x),SSU_DB,ITS_DB),
                                        multithread=parallel::detectCores(),
                                        random.seed=666,
                                        try.rc = TRUE,
                                        min.boot=50)
    # export as RDS
    saveRDS(tax,outfile)
    
  } else {next}
}

