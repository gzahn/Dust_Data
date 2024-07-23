# SETUP ####

## packages ####
library(tidyverse)
library(dada2)

## functions ####
source("./R/functions.R")

## file paths ####
  # asv tables
run_1_its <- readRDS("./data/ASV_Tables/Run_1_ITS_ASV_Table.RDS")
run_2_its <- readRDS("./data/ASV_Tables/Run_2_ITS_ASV_Table.RDS")
run_3_its <- readRDS("./data/ASV_Tables/Run_3_ITS_ASV_Table.RDS")
run_4_its <- readRDS("./data/ASV_Tables/Run_4_ITS_ASV_Table.RDS")
run_5_its <- readRDS("./data/ASV_Tables/Run_5_ITS_ASV_Table.RDS")
run_6_its <- readRDS("./data/ASV_Tables/Run_6_ITS_ASV_Table.RDS")
run_1_ssu <- readRDS("./data/ASV_Tables/Run_1_SSU_ASV_Table.RDS")
run_2_ssu <- readRDS("./data/ASV_Tables/Run_2_SSU_ASV_Table.RDS")
run_3_ssu <- readRDS("./data/ASV_Tables/Run_3_SSU_ASV_Table.RDS")
run_4_ssu <- readRDS("./data/ASV_Tables/Run_4_SSU_ASV_Table.RDS")
run_5_ssu <- readRDS("./data/ASV_Tables/Run_5_SSU_ASV_Table.RDS")
run_6_ssu <- readRDS("./data/ASV_Tables/Run_6_SSU_ASV_Table.RDS")

# databases
ITS_DB <- "./taxonomy/Eukaryome_General_ITS_v1.8_reformatted.fasta.gz"
SSU_DB <- "./taxonomy/Eukaryome_General_SSU_v1.8_reformatted.fasta.gz"

# output files
run_1_its_out <- "./data/ASV_Tables/Run_1_ITS_Taxonomy_Table.RDS"
run_2_its_out <- "./data/ASV_Tables/Run_2_ITS_Taxonomy_Table.RDS"
run_3_its_out <- "./data/ASV_Tables/Run_3_ITS_Taxonomy_Table.RDS"
run_4_its_out <- "./data/ASV_Tables/Run_4_ITS_Taxonomy_Table.RDS"
run_5_its_out <- "./data/ASV_Tables/Run_5_ITS_Taxonomy_Table.RDS"
run_6_its_out <- "./data/ASV_Tables/Run_6_ITS_Taxonomy_Table.RDS"
run_1_ssu_out <- "./data/ASV_Tables/Run_1_SSU_Taxonomy_Table.RDS"
run_2_ssu_out <- "./data/ASV_Tables/Run_2_SSU_Taxonomy_Table.RDS"
run_3_ssu_out <- "./data/ASV_Tables/Run_3_SSU_Taxonomy_Table.RDS"
run_4_ssu_out <- "./data/ASV_Tables/Run_4_SSU_Taxonomy_Table.RDS"
run_5_ssu_out <- "./data/ASV_Tables/Run_5_SSU_Taxonomy_Table.RDS"
run_6_ssu_out <- "./data/ASV_Tables/Run_6_SSU_Taxonomy_Table.RDS"



# ASSIGN TAXONOMY AND EXPORT ####

## ITS Tables ####

# assign taxonomy
x <- assign_taxonomy_to_asv_table(asv.table=run_1_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_1_its_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_2_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_2_its_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_3_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_3_its_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_4_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_4_its_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_5_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_5_its_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_6_its,
                                  tax.database=ITS_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_6_its_out)


## SSU Tables ####

# assign taxonomy
x <- assign_taxonomy_to_asv_table(asv.table=run_1_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_1_ssu_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_2_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_2_ssu_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_3_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_3_ssu_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_4_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_4_ssu_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_5_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_5_ssu_out)

x <- assign_taxonomy_to_asv_table(asv.table=run_6_ssu,
                                  tax.database=SSU_DB,
                                  multithread=parallel::detectCores(),
                                  random.seed=666,
                                  try.rc = TRUE,
                                  min.boot=50)
# export as RDS
saveRDS(x,run_6_ssu_out)
