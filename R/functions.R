
# ra() ####
# relative abundance transformation
ra <- function(x){x/sum(x)}

# plot_bar2() ####
# phyloseq bar plot without lines for each ASV
plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
          title = NULL, facet_grid = NULL, width = 0.9) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", width = width)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# remove_primers() ####
# Function to remove primers from raw amplicon files

remove_primers <- function(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                           amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                           amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                           sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
                           fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to raw data
                           rev.fp.colname = "rev_filepath",
                           fwd_pattern="_R1_",
                           rev_pattern="_R2_",
                           fwd_primer="GTGCCAGCMGCCGCGGTAA",
                           rev_primer="GGACTACHVGGGTWTCTAAT",
                           multithread=parallel::detectCores()-1){
  
  library(tidyverse); packageVersion("tidyverse")
  library(dada2); packageVersion("dada2")
  library(purrr); packageVersion("purrr")
  library(Biostrings); packageVersion("Biostrings")
  library(ShortRead); packageVersion("ShortRead")
  library(parallel); packageVersion("parallel")
  
  # tests
  if(all(is.na(metadata[[fwd.fp.colname]]))){
    stop("Fwd filepath column has no filenames in it.")
  }
  
  # File parsing
  
  # subset metadata to just that amplicon and get fwd and rev file paths
  x <- metadata[metadata[[amplicon.colname]] == amplicon & !is.na(metadata[[fwd.fp.colname]]),]
  # x <- metadata2[metadata2[[amplicon.colname]] == amplicon & !is.na(metadata2[[fwd.fp.colname]]),]
  fnFs <- x[[fwd.fp.colname]]
  fnRs <- x[[rev.fp.colname]]

  FWD <- fwd_primer # Sequence of FWD primer
  REV <- rev_primer  # Sequence of REV primer
  
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  # Prefilter to remove reads with ambiguous (N) bases ####
  
  # build new directory names
  filtN_paths <- dirname(fnFs)
  filtN_paths <- filtN_paths[!is.na(filtN_paths)]
  filtN_paths <- file.path(filtN_paths,"filtN")
  
  # get sample_names
  sample_names <- x[[sampleid.colname]]
  # make output file paths
  fnFs.filtN <- file.path(filtN_paths, paste0(sample_names,"_filtN_fwd.fastq.gz")) # Put N-filterd files in filtN/ subdirectory
  fnRs.filtN <- file.path(filtN_paths, paste0(sample_names,"_filtN_rev.fastq.gz"))
  
  # create directories, if needed
  # create new directories as needed
  for(i in unique(filtN_paths)){
    if(!file_test("-d", i)){dir.create(i)}
  }
  
  # Remove reads with any Ns
  filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, 
                maxN = 0, 
                multithread = ifelse(multithread>1,TRUE,FALSE)) # on Windows, set multithread = FALSE
  
  # build cutadapt file structure
  path.cut <- file.path(dirname(fnFs),"cutadapt")
  # create directories as needed
  for(i in unique(path.cut)){
    if(!dir.exists(i)) dir.create(i)
  }
  
  # build file names for cutadapt output
  fnFs.cut <- file.path(path.cut, paste0(sample_names,"_cutadapt_fwd.fastq.gz"))
  fnRs.cut <- file.path(path.cut, paste0(sample_names,"_cutadapt_rev.fastq.gz"))
  
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run Cutadapt
  for(i in seq_along(fnFs)) {
    system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", "--cores 0", # -n 2 required to remove FWD and REV from reads
                                 "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                                 fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
  
  


}


# run_itsxpress() ####
# Isolate ITS region 

# Valid taxa_group arguments: 
# Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All

run_itsxpress <- function(directory, # where cutadapted reads live
                          itsregion="ITS1", # must be "ITS1" or "ITS2"
                          taxa_group="All",
                          nthreads=(parallel::detectCores()-1),
                          fwd_pattern="_R1_"){
  
  # find the "cutadapted" files
  fwds <- list.files(directory,pattern = fwd_pattern,full.names = TRUE)
  # build names for outfiles
  outs <- paste0(tools::file_path_sans_ext(fwds) %>% 
                   tools::file_path_sans_ext(), 
                 "_ITS.fastq.gz")
  its_dir <- directory %>% str_replace('cutadapt','ITSx')

  if(!dir.exists(its_dir)){dir.create(its_dir)}
  
  outs <- file.path(its_dir,basename(outs))
  
  # build the ITSxpress command and run it on each file in turn

  for(i in 1:length(fwds)){
      itsxpress <- paste0("itsxpress --fastq ",fwds[i],
                        " --outfile ",outs[i],
                        " --region ",itsregion,
                        " --taxa ",taxa_group,
                        " --threads ",nthreads,
                        " --log ",outs[i],".log",
                        " --single_end")
    
    system(command = itsxpress)
  }
  
}



# build_asv_table() ####
# function to run dada2 pipeline on trimmed amplicon reads
# error profiling and correction should be done by sequencing run
# this function will run dada2 on files from a single sequencing run, given a metadata sheet that has samples from multiple runs

build_asv_table <- function(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                            run.id.colname = "run_id", # name of column in metadata indicating which sequencing run a sample comes from
                            run.id, # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
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
                            asv.table.dir = "./ASV_Tables", # path to directory where final ASV table will be saved
                            random.seed = 666){
  
  # tests
  stopifnot("data.frame" %in% class(metadata))
  if(is.null(metadata[[run.id.colname]])){stop("run.id.colname is empty or not found in metadata.")}
  if(class(run.id) != 'character'){stop("run.id must be a character vector of length 1.")}
  if(class(fwd.fp.colname) != 'character'){stop("fwd.fp.colname must be a character vector of length 1, matching the column name in metadata.")}
  if(single.end & class(rev.fp.colname) != 'character'){stop("rev.fp.colname must be a character vector of length 1, matching the column name in metadata.")}
  
  # parse options
  if(single.end){
    paired <- FALSE
    maxEE <- maxEE[1]
  } else {
    paired <- TRUE
    maxEE <- maxEE
  }
  set.seed(random.seed)
  
  # subset to sequencing run and only rows with file names
  metadata <- 
  metadata[as.character(metadata[[run.id.colname]]) == as.character(run.id) & 
             !is.na(metadata[[fwd.fp.colname]]) & metadata[[amplicon.colname]] == amplicon,]
  
  
  # parse filepaths and sample names
  fns <- metadata[[fwd.fp.colname]]
  if(paired){rns <- metadata[[rev.fp.colname]]}
  
  sample_names <- metadata[[sampleid.colname]]
  names(fns) <- sample_names
  if(paired){names(rns) <- sample_names}
  
  filtered_paths <- file.path(dirname(fns),filtered.dir)
  
  filts_f <- file.path(filtered_paths,paste0(sample_names,"_filtered_fwd.fastq.gz"))
  names(filts_f) <- sample_names
  if(paired){
    filts_r <- file.path(filtered_paths,paste0(sample_names,"_filtered_rev.fastq.gz"))
    names(filts_r) <- sample_names
  }
  
  
  # create new directories as needed
  for(i in unique(filtered_paths)){
    if(!file_test("-d", i)){dir.create(i)}
  }
  
  
  
  
  # filter and trim
  if(paired){
    out <- filterAndTrim(fns, filts_f, rns, filts_r, 
                       maxN=0, 
                       maxEE=maxEE, 
                       truncQ=truncQ,
                       rm.phix=rm.phix, 
                       compress=compress,
                       multithread=multithread)
  } else {
    out <- filterAndTrim(fns, filts_f,
                         maxN=0,
                         maxEE=maxEE[1],
                         truncQ=truncQ,
                         rm.phix=rm.phix,
                         compress=compress,
                         multithread=multithread)
    }
    
  # In case some samples may have had zero reads pass QC, reassign filts
  # for each unique filtered path, go in, find all files, stick them together in one big vector,
  # then sort them to match the metadata sheet by finding any missing and removing those from the filts_f filts_r and samplenames

  if(paired){
    new_filts_f <- c()
    new_filts_r <- c()
    for(i in unique(filtered_paths)){
      x <- list.files(i,full.names = TRUE,pattern = "_filtered_fwd.fastq.gz")  
      y <- list.files(i,full.names = TRUE,pattern = "_filtered_rev.fastq.gz")  
      new_filts_f <- c(new_filts_f,x)
      new_filts_r <- c(new_filts_r,y)
    }
  } else {
    new_filts_f <- c()
    # new_filts_r <- c()
    for(i in unique(filtered_paths)){
      x <- list.files(i,full.names = TRUE,pattern = "_filtered_fwd.fastq.gz")  
      # y <- list.files(i,full.names = TRUE,pattern = "_filtered_rev.fastq.gz")  
      new_filts_f <- c(new_filts_f,x)
      # new_filts_r <- c(new_filts_r,x)
    }
  }
  
  
  # find any missing and pull them out
  filts_f <- filts_f[which((filts_f %in% new_filts_f))]
  if(paired){filts_r <- filts_r[which((filts_r %in% new_filts_r))]}
  
  
  
  # learn errors
  errF <- learnErrors(filts_f, multithread=ifelse(multithread>1,TRUE,FALSE), 
                      MAX_CONSIST = 20,verbose = 1,
                      randomize = TRUE) # set multithread = FALSE on Windows
  errF_out <- paste0("Run_",as.character(run.id),"_",amplicon,"_err_Fwd.RDS")
  saveRDS(errF,file.path(asv.table.dir,errF_out))
  
  
  if(paired){
  errR <- learnErrors(filts_r, multithread=ifelse(multithread>1,TRUE,FALSE), 
                      MAX_CONSIST = 20,verbose = 1,
                      randomize = TRUE) # set multithread = FALSE on Windows
  errR_out <- paste0("Run_",as.character(run.id),"_",amplicon,"_err_Rev.RDS")
  saveRDS(errR,file.path(asv.table.dir,errR_out))
  }
    
  
  # dereplication
  derepF <- derepFastq(filts_f, verbose=FALSE)
  if(paired){derepR <- derepFastq(filts_r, verbose=FALSE)}
  
  
  # SAMPLE INFERRENCE ####
  dadaFs <- dada(derepF, err=errF, multithread=ifelse(multithread>1,TRUE,FALSE), 
                 selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
  if(paired){
    dadaRs <- dada(derepR, err=errR, multithread=ifelse(multithread>1,TRUE,FALSE), 
                   selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows  
  }
  
  
  # MERGE FWD and REV READS ####
  if(paired){
    mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=FALSE)
  } 
    
  # MAKE SEQUENCE TABLE ####
  if(paired){
    seqtab <- makeSequenceTable(mergers)
  } else {
    seqtab <- makeSequenceTable(dadaFs)
    }
  
  
  # REMOVE CHIMERAS ####
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=ifelse(multithread>1,TRUE,FALSE), verbose=TRUE)
  
  # REMOVE CONTAMINANTS ####
  
  # make metadata match files that made it through all previous steps!
  
  # find negative control samples, if any
  metadata[["control"]] <- metadata[["sample_type"]] == "neg_control"
  
  
  # only run if there are negative control(s) that have at least some reads
  if(any(metadata[["control"]]) & sum(seqtab.nochim[metadata[["control"]],]) > 0){
    # Find and remove contaminants
    contams = decontam::isContaminant(seqtab.nochim, neg = metadata$control, normalize = TRUE)

    # remove contaminant sequences and control samples from both tables, respectively ####
    seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
    seqtab.nochim = seqtab.nochim[!metadata$control,]
    print(paste0(sum(contams$contaminant)," likely contaminants were removed from the ASV table of Run ",as.character(run.id),"."))
  }
  
  # make output name for ASV table
  asv_out <- paste0(asv.table.dir,"Run_",as.character(run.id),"_",amplicon,"_ASV_Table.RDS")
  
  saveRDS(seqtab.nochim,file.path(asv.table.dir,asv_out))

}



# assign_taxonomy_to_asv_table()
# function to assign taxonomy to an asv table from the dada2 pipeline
# inputs: asv table | path-to-database | 

assign_taxonomy_to_asv_table <- function(asv.table, # asv table object name
                                         tax.database, # path to taxonomic database (fasta)
                                         multithread=(parallel::detectCores()-1), # set to FALSE on Windows
                                         random.seed=666,
                                         try.rc = TRUE, # attempt revComplement assignments as well? (doubles time)
                                         min.boot=50 # bootstrap of 50% recommended for seqs shorter than 250nt
                                         ){
  x <- assignTaxonomy(seqs = asv.table,
                      refFasta = tax.database,
                      minBoot = min.boot,
                      tryRC = try.rc,
                      outputBootstraps = FALSE,
                      multithread = multithread,
                      verbose = FALSE)
  
  return(x)
  
}
