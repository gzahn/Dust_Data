#!/bin/bash

# this will download and build the custom databases
# The Eukaryome database, supplemented with maarjaam VXT assignments for glomeromycota taxa

# Run this from the directory where you will keep your databases for your project...


# get Eukaryome databases
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_SSU_v1.8.zip
unzip General_EUK_SSU_v1.8.zip && rm General_EUK_SSU_v1.8.zip
wget https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v1.8.zip
unzip General_EUK_ITS_v1.8.zip && rm General_EUK_ITS_v1.8.zip


# reformat Eukaryome to play nicely with dada2 assignTaxonomy() function
zcat General_EUK_SSU_v1.8.zip | sed 's/>\([^;]*;\)/>/' | sed 's/>\([^|]*|\)/>/' | sed 's/|/;/g' > Eukaryome_General_SSU_v1.8_reformatted.fasta
zcat General_EUK_ITS_v1.8.zip | sed 's/>\([^;]*;\)/>/' | sed 's/>\([^|]*|\)/>/' | sed 's/|/;/g' > Eukaryome_General_ITS_v1.8_reformatted.fasta


# get maarjAM databases (fasta)
wget https://maarjam.ut.ee/resources/maarjam_database_SSU_TYPE.fasta.2021.zip -O maarjam_database_SSU_TYPE_2021.fasta.gz
wget https://maarjam.ut.ee/resources/maarjam_database_SSU.fasta.2021.zip -O maarjam_database_SSU_2021.fasta.gz
wget https://maarjam.ut.ee/resources/maarjam_database_full_ITS.fasta.2021.zip -O maarjam_database_full_ITS_2021.fasta.gz

# decompress
gunzip maarjam_database_*gz

# Run custom R script to add taxonomy to maarjam databases
Rscript ./format_maarjAM.R

# combine the databases to new files
cat maarjam_database_SSU*newheaders.fasta Eukaryome_General_SSU_v1.8_reformatted.fasta > Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta
cat maarjam_database_full_ITS_2021_newheaders.fasta Eukaryome_General_ITS_v1.8_reformatted.fasta > Eukaryome_General_ITS_v1.8_reformatted_maarjam.fasta

# compress the final files
gzip *.fasta

