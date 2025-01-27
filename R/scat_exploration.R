# SETUP ####

# packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(ecodist)
library(geosphere)
library(gdm)
library(patchwork)
library(fungaltraits)
library(FUNGuildR)
library(ggmap)
library(ggimage)
library(dada2)

# functions and themes
source("./R/functions.R")

set.seed(666)

theme_set(theme_bw() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face='bold',size=16),
                  axis.title = element_text(face='bold',size=16),
                  axis.text = element_text(face='bold',size=12),
                  legend.title = element_text(face='bold',size=16),
                  legend.text = element_text(face='bold',size=12)))
# PREP DATA ####

# metadata
meta <- readRDS("./data/full_clean_metadata.RDS")

# ASV tables (from Sarah)
its <- read.csv("./data/ASV_Tables/ASV_Table_ITS_SCAT.csv")
ssu <- read.csv("./data/ASV_Tables/ASV_Table_SSU_SCAT.csv")

# ASV sequences
its_seqs <- colnames(its)[-1]
ssu_seqs <- colnames(ssu)[-1]

# sample names
its_samplenames <- its[,1]
ssu_samplenames <- ssu[,1]
identical(its_samplenames,ssu_samplenames)
its_ids <- paste0(its_samplenames,"_ITS")
ssu_ids <- paste0(ssu_samplenames,"_SSU")


# strip off first columns
its <- its[,-1]
ssu <- ssu[,-1]

# subset metadata
its_meta <- meta[meta$library_id %in% its_ids,]
ssu_meta <- meta[meta$library_id %in% ssu_ids,]
# reorder to match
its_meta <- its_meta[its_meta$library_id %>% order(its_ids),]
ssu_meta <- ssu_meta[ssu_meta$library_id %>% order(ssu_ids),]


## assign taxonomy ####
its_tax <- assignTaxonomy(its_seqs,
                          refFasta = "./taxonomy/Eukaryome_General_ITS_v1.8_reformatted_maarjam.fasta.gz",
                          multithread = TRUE,verbose = TRUE)
saveRDS(its_tax,"./data/ASV_Tables/Scat_ITS_Taxonomy_Table.RDS")
ssu_tax <- assignTaxonomy(ssu_seqs,
                          refFasta = "./taxonomy/Eukaryome_General_SSU_v1.8_reformatted_VTX.fasta.gz",
                          multithread = TRUE,verbose = TRUE)
saveRDS(ssu_tax,"./data/ASV_Tables/Scat_SSU_Taxonomy_Table.RDS")



## build phyloseq ####
# metadata
its_met <- sample_data(its_meta)
sample_names(its_met) <- its_met$library_id
ssu_met <- sample_data(ssu_meta)
sample_names(ssu_met) <- ssu_met$library_id
# otu tables
its_otu <- otu_table(its,taxa_are_rows = FALSE)
sample_names(its_otu) <- sample_names(its_met)
ssu_otu <- otu_table(ssu,taxa_are_rows = FALSE)
sample_names(ssu_otu) <- sample_names(ssu_met)
# taxonomy tables
its_tax <- readRDS("./data/ASV_Tables/Scat_ITS_Taxonomy_Table.RDS")
ssu_tax <- readRDS("./data/ASV_Tables/Scat_SSU_Taxonomy_Table.RDS")
its_tax <- tax_table(its_tax)
ssu_tax <- tax_table(ssu_tax)
# build
its_ps <- phyloseq(its_otu,its_met,its_tax) %>% clean_ps_taxonomy()
ssu_ps <- phyloseq(ssu_otu,ssu_met,ssu_tax) %>% clean_ps_taxonomy()
# remove nonfungi & nonamf, respectively
its_ps <- its_ps %>% 
  subset_taxa(Kingdom == "Fungi")
# check for any glomeromycota in ITS and remove if found
its_ps@tax_table %>% 
  as.data.frame() %>% 
  map(unique) %>% grep(pattern="glomero",ignore.case = TRUE)
ssu_ps <- ssu_ps %>% 
  subset_taxa(Phylum == "Glomeromycota")
its_ps <- its_ps %>% 
  subset_taxa(taxa_sums(its_ps) > 0)
ssu_ps <- ssu_ps %>% 
  subset_taxa(taxa_sums(ssu_ps) > 0)




# export separate ps objects
saveRDS(its_ps,"./data/physeq_objects/scat_ITS_physeq.RDS")
saveRDS(ssu_ps,"./data/physeq_objects/scat_SSU_physeq.RDS")

# merge
# rename samples before merging
sample_names(its_ps) <- its_samplenames
sample_names(ssu_ps) <- ssu_samplenames
# remove amplicon info
its_ps@sam_data$amplicon <- NULL
ssu_ps@sam_data$amplicon <- NULL
# merge
ps <- merge_phyloseq(its_ps,ssu_ps)

# combine taxa at species level (not ASV level)
ps_species <- tax_glom(ps,taxrank = "Species",NArm = FALSE,bad_empty=c(NA, "", " ", "\t","unclassified"))

## export ####
saveRDS(ps,"./data/physeq_objects/scat_FULL_physeq.RDS")
saveRDS(ps_species,"./data/physeq_objects/scat_FULL_physeq_species-level.RDS")
ps
ps_species


# (reload point here) ####

# load 
ps <- readRDS("./data/physeq_objects/scat_FULL_physeq.RDS") # combined data
ps_species <- readRDS("./data/physeq_objects/scat_FULL_physeq_species-level.RDS") # combined, species-level
its_ps <- readRDS("./data/physeq_objects/scat_ITS_physeq.RDS") # just its fungi
ssu_ps <- readRDS("./data/physeq_objects/scat_SSU_physeq.RDS") # just ssu amf

# light cleaning
ps@sam_data$year <- ps@sam_data$year %>% factor(levels = c("2019","2021","2022"))
ps_species@sam_data$year <- ps_species@sam_data$year %>% factor(levels = c("2019","2021","2022"))
its_ps@sam_data$year <- its_ps@sam_data$year %>% factor(levels = c("2019","2021","2022"))
ssu_ps@sam_data$year <- ssu_ps@sam_data$year %>% factor(levels = c("2019","2021","2022"))

ps_species@tax_table[,1] %>% unique %>% unname

ssu_ps %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Family") +
  facet_wrap(~year,scales='free')
amf_list <- ssu_ps@tax_table[,c(6,7)] %>% unname %>% 
  as.data.frame %>% 
  mutate(taxon = paste0(V1," ",V2))
amf_list$asv = taxa_names(ssu_ps)
ssu_melt <- psmelt(ssu_ps)
its_melt <- psmelt(its_ps)
# export files for Sarah
write_csv(ssu_melt,"./output/scat_amf_melt.csv")
write_csv(its_melt,"./output/scat_its_melt.csv")
write_csv(amf_list,"./output/scat_amf_taxa_list.csv")


## Guilds ####

# get funguild database
funguild <- FUNGuildR::get_funguild_db()

# get fungaltraits database
fungaltraits_db <- fungal_traits()

ps_species@sam_data$sample_type %>% unique()
# get search queries
query_taxonomy <- 
  data.frame(
    asv = taxa_names(ps_species),
    Genus = ps_species@tax_table[,6],
    Species = ps_species@tax_table[,7]
  ) %>% 
  mutate(taxon = paste0(Genus," ",Species) %>% 
           str_remove(" NA") %>% 
           str_remove(" VTX.*") %>% 
           str_remove(" unclassified") %>% 
           str_remove("_gen.*"))


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

guild_assignments <- 
  query_taxonomy %>% 
  full_join(funguild_df) 

# if they all match up in the same order as physeq tax_table...
if(all(taxa_names(ps_species) == guild_assignments$asv)){
  # add guild to phyloseq object
  ps_species@tax_table[,1] <- guild_assignments$funguild_guild
} else (print("Out of order...double-check!"))
colnames(ps_species@tax_table)[1] <- "Guild"

# Simplify fungal guilds to major groupings on relabund transformed data
ps_species@sam_data %>% 
  as('data.frame') %>% map(class)


ps_species@sam_data$total_precip_jun_nov_mm <- NULL # get rid of this 'list' col for now. just need to hurry


# make merge variable
ps_species@sam_data$mergevar <- paste(ps_species@sam_data$site,
                                      ps_species@sam_data$am_em_dom,
                                      sep="_")
# merge
ps_guilds_ra_melt <- 
  ps_species %>% 
  merge_samples('mergevar')
# repair metadata
ps_guilds_ra_melt@sam_data$site <- sample_names(ps_guilds_ra_melt) %>% str_split("_") %>% map_chr(1)
ps_guilds_ra_melt@sam_data$am_em_dom <- sample_names(ps_guilds_ra_melt) %>% str_split("_") %>% map_chr(2)

ps_guilds_ra_melt <- 
ps_guilds_ra_melt %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  simplify_fungal_guilds()
# reorder sites by am_em_dom
siteorder <- 
ps_guilds_ra_melt %>% 
  dplyr::select(site,am_em_dom) %>% 
  unique.data.frame() %>% 
  arrange(am_em_dom) %>% 
  pluck('site')
ps_guilds_ra_melt$site <- ps_guilds_ra_melt$site %>% factor(levels=siteorder)

# plot summarized info
ps_guilds_ra_melt %>% 
  mutate(major_guild = case_when(is.na(major_guild) ~ "Unassigned",
                                 TRUE ~ major_guild)) %>% 
  ggplot(aes(x=site,y=Abundance,fill=major_guild)) +
  geom_col() +
  geom_point(aes(y=1.05,color=am_em_dom),
             shape=15,size=10) +
  theme(axis.text.x = element_text(angle=270,hjust=0,vjust=.5)) +
  labs(x="Site",y="Relative abundance",fill="Guild",color="Dominant\nplant type") +
  scale_fill_viridis_d(option = 'turbo') +
  scale_color_viridis_d(option = 'mako',begin=.2,end=.9) +
  guides(fill = guide_legend(override.aes = list(color = NA)))



## Ordinations ####
ord <- 
  ps_species %>% 
  subset_samples(sample_sums(ps_species) > 0) %>% 
  transform_sample_counts(ra) %>% 
  ordinate(method = "PCoA")
plot_ordination(ps_species,ord,
                color = "year")
ps_species@sam_data %>% names


ssu_ord <- 
  ssu_ps %>% 
  subset_samples(sample_sums(ssu_ps) > 0) %>% 
  transform_sample_counts(ra) %>% 
  ordinate(method = "NMDS")
plot_ordination(ssu_ps,ssu_ord,
                color = "year")
