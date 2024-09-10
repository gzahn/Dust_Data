# Dust_Data

Code to process and analyze Chaudhary Lab dust and soil data from NEON sites.


Full pipeline can be run via:

```./pipeline.sh```


Tested on Ubuntu 22.04 and Red Hat 8.10 Linux systems using R v 4.4.0


# Steps:

  1. Build the SSU & ITS taxonomy databases

  This will download the Eukaryome databases, format them for DADA2 classification, and append the Marjaam VTX taxa to the SSU database.
  
     ```
     cd ./taxonomy/
     ./build_taxonomy_database.sh
     cd -
     ```
  2.  R/00_Download_Raw_Seq_Data.R
  
  This script downloads the raw sequence data from the Sequence Read Archive (SRA - PRJNAXXXXXXX)

  3. R/01_Clean_Sample_Metadata.R

  This script prepares and cleans the sample metadata

  4. R/02_Trim_Adaptors_and_Flanking.R

  This script trims Illumina adaptors from all reads and removes flanking 18S regions from ITS reads

  5. R/03_Build_ASV_Tables.R

  This script builds ASV tables using DADA2 from trimmed reads. Possible contaminant sequences are detected from negative controls and removed from samples.
  Separate ASV tables are computed for each sequencing run in the project; SSU and ITS amplicons are handled separately.

  6. R/04_Assign_Taxonomy.R

  This script assigns taxonomy to all ASV tables using our constructed SSU and ITS databases, respectively.

  7. R/05_Build_Physeq_Objects.R

  This script combines cleaned sample metadata, ASV tables, & taxonomic assignments into phyloseq objects. Separate phyloseq objects are constructed for each ASV table, and can be combined downstream if desired.
  


___



List of sequence data file paths:

./file_list.txt

This file was generated with

```find /dartfs/rc/lab/C/ChaudharyB/SequencingRuns -name "*.fastq.gz" -type f > file_list.txt```
