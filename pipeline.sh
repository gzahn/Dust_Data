#!/bin/bash


# Tested on Ubuntu 22.04 and RedHat 8.1.8 Linux systems

Rscript R/00_Download_Sample_Metadata.R
Rscript R/00_Download_Raw_Seq_Data.R
Rscript R/01_Clean_Sample_Metadata.R
Rscript R/02_Trim_Adaptors_and_Flanking.R
Rscript R/03_Build_ASV_Tables.R
Rscript R/04_Assign_Taxonomy.R

