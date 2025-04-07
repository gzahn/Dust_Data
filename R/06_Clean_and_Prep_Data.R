# SETUP ####

## packages ####
library(tidyverse)
library(phyloseq)
library(fungaltraits)

## functions ####
source("./R/functions.R")

## paths ####
asv_table_dir <- "./data/ASV_Tables"
ssu_ps_path <- "./data/physeq_objects/full_ssu_ps_raw.RDS"
its_ps_path <- "./data/physeq_objects/full_its_ps_raw.RDS"

## data ####
ssu <- readRDS(ssu_ps_path)
its <- readRDS(its_ps_path)

# CLEAN ####

# remove list-cols from metadata
its@sam_data$total_precip_jun_nov_mm <- NULL
ssu@sam_data$total_precip_jun_nov_mm <- NULL

## subset to sample types of interest (dust, soil?) ####
its <- 
  its %>% 
  subset_samples(sample_type %in% c("dust","soil"))

ssu <- 
  ssu %>% 
  subset_samples(sample_type %in% c("dust","soil"))

## subset ITS to only fungi #### 
its <- clean_ps_taxonomy(its)
its <- its %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum))

## find guild assignments ####

# get funguild database
funguild <- FUNGuildR::get_funguild_db()

# get fungaltraits database
fungaltraits_db <- fungal_traits()

# get search queries
query_taxonomy <- 
  data.frame(
    asv = taxa_names(its),
    Class = its@tax_table[,3],
    Order = its@tax_table[,4],
    Family = its@tax_table[,5],
    Genus = its@tax_table[,6] %>% str_remove("(Fungi)"),
    Species = its@tax_table[,7]
  ) %>% 
  mutate(Genus = case_when(Genus == "unclassified" ~ NA,
                           TRUE ~ Genus),
         Species = case_when(Species == "unclassified" ~ "sp.",
                             TRUE ~ Species)) %>% 
  mutate(taxon = paste0(Genus," ",Species) %>% 
           str_remove(" NA") %>% 
           str_remove(" VTX.*") %>% 
           str_remove(" unclassified") %>% 
           str_remove("_gen.*") %>% 
           str_remove(".var.*"))

# build funguild assignment list
guilds <- list()
for(query in seq_along(query_taxonomy$taxon)){
  rows <- funguild[funguild$taxon == query_taxonomy$taxon[query],]
  if(nrow(rows) < 1){
    rows <- funguild[funguild$taxon == query_taxonomy$taxon[query] %>% str_split(" ") %>% map_chr(1),] # try just genus
  }
  # if still no assignment found in funguild, move on
  if(nrow(rows) < 1){guilds[query_taxonomy$asv[query]] <- NA;next}
  if(nrow(rows) > 0){
    guilds[query_taxonomy$asv[query]] <- unique(rows[["guild"]])
  }
}


funguild_df <- 
  data.frame(asv = names(guilds),
             funguild_guild = unlist(guilds))

# get unique list of mycorrizal guild names
mycorrhizal_guild_names <- 
  grep(pattern = "mycorrhiz",x = unique(funguild_df$funguild_guild),ignore.case = TRUE,value = TRUE)


# compare funguild w/ fungaltraits
fungaltraits_db$guild_fg
fungaltraits_db$Genus
  # convert "unclassified" to NA
query_taxonomy$Genus[query_taxonomy$Genus == "unclassified"] <- NA


guild_assignments <- 
  query_taxonomy %>% 
  full_join(funguild_df) %>% 
  mutate(Guild = funguild_guild) %>% 
  simplify_fungal_guilds()


ft_guild_lookup <- 
  fungaltraits_db %>% 
  dplyr::select(Genus,guild_fg) %>% 
  unique.data.frame() %>% 
  dplyr::filter(!is.na(guild_fg))

left_join(guild_assignments, ft_guild_lookup) %>% 
  dplyr::filter(grepl(pattern="ycorrhiz",major_guild))
# FunGuild assignments are more complete, stick with them!


# if they all match up in the same order as physeq tax_table...
if(all(taxa_names(its) == guild_assignments$asv)){
  #                                      ...add guild to phyloseq object
  its@tax_table[,1] <- guild_assignments$major_guild
} else (print("Out of order...double-check!"))


its@tax_table[,1] %>% unique %>% unname
# remove spurious capitalization
its@tax_table[,1] <- its@tax_table[,1] %>% str_to_sentence()

## subset to only mycorrhizal (ECM, Ericoid)
# rename kingdom to "guild" first
colnames(tax_table(its))[1] <- "Guild"
## Examine guild assignments for all fungal taxa ####

## save copies of physeqs containing all fungal taxa ####
its_alltaxa <- its
ssu_alltaxa <- clean_ps_taxonomy(ssu) %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum))
colnames(tax_table(ssu_alltaxa))[1] <- "Guild"
ssu_alltaxa@tax_table[which(ssu_alltaxa@tax_table[,2] == "Glomeromycota"),1] <- "Arbuscular mycorrhizal"
ssu_alltaxa@tax_table[which(ssu_alltaxa@tax_table[,2] != "Glomeromycota"),1] <- "Other"
its_alltaxa@tax_table[grep(pattern="mycorrhiz",x=its_alltaxa@tax_table[,1],invert = TRUE),1] <- "Other"
ssu_alltaxa@tax_table[,1] %>% table
its_alltaxa@tax_table[,1] %>% table

saveRDS(its_alltaxa,"./data/physeq_objects/full_its_ps_alltaxa_clean.RDS")
saveRDS(ssu_alltaxa,"./data/physeq_objects/full_ssu_ps_alltaxa_clean.RDS")

# I would like to model the probability of recovering ANY mycorrhizal taxa
  # merge two objects
# melting is too slow...
# its_alltaxa_melt <- psmelt(its_alltaxa)
# ssu_alltaxa_melt <- psmelt(ssu_alltaxa) %>% 
#   dplyr::filter(Guild == "Fungi")
# build data frame by hand...
its_alltaxa@tax_table[,1] %>% unique %>% unname
ssu_alltaxa@tax_table[,1] %>% table
am <- ssu_alltaxa@otu_table[,ssu_alltaxa@tax_table[,1] == "Arbuscular mycorrhizal"] %>% rowSums()
ecm <- its_alltaxa@otu_table[,its_alltaxa@tax_table[,1] == "Ectomycorrhizal"] %>% rowSums()
erm <- its_alltaxa@otu_table[,its_alltaxa@tax_table[,1] == "Ericoid mycorrhizal"] %>% rowSums()
om <- its_alltaxa@otu_table[,its_alltaxa@tax_table[,1] == "Orchid mycorrhizal"] %>% rowSums()
am[am>1] <- 1
ecm[ecm>1] <- 1
erm[erm>1] <- 1
om[om>1] <- 1

e_mycorrhizal_presence_lgl_df <- 
  data.frame(
  index = names(ecm),
  ecm = as.logical(ecm),
  erm = as.logical(erm),
  om = as.logical(om)
)
all(sample_names(its) == e_mycorrhizal_presence_lgl_df$index)
its@sam_data$ecm <- e_mycorrhizal_presence_lgl_df$ecm
its@sam_data$erm <- e_mycorrhizal_presence_lgl_df$erm
its@sam_data$om <- e_mycorrhizal_presence_lgl_df$om



a_mycorrhizal_presence_lgl_df <- 
  data.frame(
    index = names(am),
    am = as.logical(am)
  )
all(sample_names(ssu) == a_mycorrhizal_presence_lgl_df$index)
ssu@sam_data$am <- a_mycorrhizal_presence_lgl_df$am
    
# save data frames for random forest
saveRDS(its@sam_data %>% as("data.frame"),"./output/mycorrhizal_type_logical_df_em.RDS")
saveRDS(ssu@sam_data %>% as("data.frame"),"./output/mycorrhizal_type_logical_df_am.RDS")



# barplot of guild assignments
its_alltaxa %>% 
  merge_samples("sample_type") %>%  
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Guild")
ggsave("./output/figs/ITS_guild_assignments_relabund.png")

data.frame(
  genus = its_alltaxa@tax_table[,6],
  guild = its_alltaxa@tax_table[,1]) %>% 
  unique.data.frame() %>% 
  dplyr::filter(is.na(Guild) & !is.na(Genus)) %>% 
  write_csv("./output/ITS_genera_w_unknown_guilds.csv")


# subset to mycorrhizal only
its@tax_table[,1] %>% unique %>% unname
its <- 
  its %>% 
  subset_taxa(grepl(pattern="mycorrhizal",its@tax_table[,1]))

## subset SSU to only AM fungi ####
ssu <- clean_ps_taxonomy(ssu)
ssu <- ssu %>% 
  subset_taxa(Phylum == "Glomeromycota")
  
### sanity checks
# check that only mycorrhizal fungi remain in both physeq objects
if(ssu@tax_table[,1] %>% unique %>% as.character() != "Fungi" | 
   !any(ssu@tax_table[,4] %>% unique %>% as.character() != "Glomerales") |
   !any(grepl(pattern = "mycorrhizal",its@tax_table[,1] %>% unique %>% as.character()))){
  stop("Recheck subsetting. Non-fungi are still present in taxonomy table(s).")
} else {
  cat("Subsetting to mycorrhizal fungi seems to have worked.")
}



## remove empty samples/taxa ####
# must have at least 3 reads (use this cutoff for now)
# This is permissive for any community analyses, but conservative for presence/absence...
its %>% 
  subset_samples(sample_sums(its) > 2)
ssu %>% 
  subset_samples(sample_sums(ssu) > 2)
# ... this REALLY cuts down on the number of samples available

## make presence/absence versions of physeqs ####

### add "amf" guild to ssu physeq
colnames(ssu@tax_table[,1]) <- "Guild"
ssu@tax_table[,1] <- "Arbuscular mycorrhizal"

its_pa <- 
  its %>% 
  subset_samples(sample_sums(its) > 0)
ssu_pa <- 
  ssu %>% 
  subset_samples(sample_sums(ssu) > 0)

# clean up to useful metadata columns
keeper_cols <- c("index","sample_id","sample_type","amplicon","library_id","site","run_id",
                 "height","height_cm","lat_dd","long_dd","am_em_dom","site_avg_elevation_m","year",
                 "mean10m_wind_annual_m_s","mean10m_minimum_wind_annual_m_s","mean10m_maximum_wind_annual_m_s",
                 "mean10m_wind_jun_nov_m_s","mean10m_minimum_wind_jun_nov_m_s","mean10m_maximum_wind_jun_nov_m_s",
                 "mean_canopy_wind_annual_m_s","mean_canopy_minimum_wind_annual_m_s","mean_canopy_maximum_wind_annual_m_s",
                 "mean_canopy_wind_jun_nov_m_s","mean_canopy_minimum_wind_jun_nov_m_s","mean_canopy_maximum_wind_jun_nov_m_s",
                 "mean_temp_annual_c","max_temp_annual_c","min_temp_annual_c","mean_temp_jun_nov_c","max_temp_jun_nov_c",
                 "min_temp_jun_nov_c","total_precip_mm")
its_pa@sam_data <- 
  its_pa@sam_data %>% 
  as('data.frame') %>% 
  dplyr::select(all_of(keeper_cols)) %>% 
  sample_data()
ssu_pa@sam_data <- 
  ssu_pa@sam_data %>% 
  as('data.frame') %>% 
  dplyr::select(all_of(keeper_cols)) %>% 
  sample_data()
# make taxonomic rank names match
colnames(ssu_pa@tax_table) <- c("Guild","Phylum","Class","Order","Family","Genus","Species")

# merged version (presence/absence for all mycorrhizal guild taxa, merged SSU and ITS)
full_pa <- 
  merge_phyloseq(its_pa,ssu_pa)

full_pa <- 
  full_pa %>% 
  transform_sample_counts(function(x){ifelse(x>0,1,0)})



# save as RDS
saveRDS(full_pa,"./data/physeq_objects/merged_ps_mycorrhizal_taxa_only_presenceabsence.RDS")

# melt phyloseq (relabund, merged)
full_melt <- psmelt(full_pa)
names(full_melt)
# save melted data frame
saveRDS(full_melt,"./data/physeq_objects/merged_ps_mycorrhizal_taxa_only_presenceabsence_melted_df.RDS")

