# SETUP ####

## packages ####
library(tidyverse)
library(dada2)
library(phyloseq)

## functions ####
source("./R/functions.R")

## paths ####
asv_table_dir <- "./data/ASV_Tables"

# find asv table files
asv_tables <- list.files(asv_table_dir,full.names = TRUE, pattern = "_ASV_Table.RDS")

# find taxonomy table files
tax_tables <- list.files(asv_table_dir,full.names = TRUE, pattern = "_Taxonomy_Table.RDS")

# out path
out_paths <- paste0(asv_tables %>% str_remove("_ASV_Table.RDS"),"_physeq_object.RDS")

# DATA ####

# get metadata
meta <- readRDS("./data/full_clean_metadata.RDS")



# BUILD PS OBJECTS ####
for(i in seq_along(asv_tables)){
  
  # load asv table
  asv <- readRDS(asv_tables[i])
  
  # load taxonomy table
  taxa <- readRDS(tax_tables[i])
  
  # subset metadata to match
  meta_subset <- meta[meta$library_id %in% row.names(asv),]
  
  # set up physeq components
  otu <- otu_table(asv,taxa_are_rows = FALSE)
  met <- sample_data(meta_subset)
  sample_names(met) <- row.names(asv)
  tax <- tax_table(taxa)
  
  # build physeq object
  physeq <- phyloseq(otu,
                     met,
                     tax
                     )

  # save physeq object
  saveRDS(physeq,out_paths[i])

}


# JOIN PHYSEQ OBJECTS ####

## SSU ####
ssu_ps_files <- list.files(asv_table_dir,full.names = TRUE,pattern="SSU_physeq_object.RDS")  
ssu_ps_list <- map(ssu_ps_files,readRDS)
ssu_ps <- purrr::reduce(ssu_ps_list,merge_phyloseq)

## ITS ####
its_ps_files <- list.files(asv_table_dir,full.names = TRUE,pattern="ITS_physeq_object.RDS")  
its_ps_list <- map(its_ps_files,readRDS)
its_ps <- purrr::reduce(its_ps_list,merge_phyloseq)

# EXPORT ####
dir.create("./data/physeq_objects")
saveRDS(ssu_ps,"./data/physeq_objects/full_ssu_ps_raw.RDS")
saveRDS(its_ps,"./data/physeq_objects/full_its_ps_raw.RDS")


