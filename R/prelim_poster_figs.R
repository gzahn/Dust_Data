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

pal.discrete <- viridis::turbo(21,begin = .1,end=.8)

# data
ssu <- readRDS("./data/physeq_objects/full_ssu_ps_raw.RDS") %>% 
  subset_samples(sample_type == "dust") # just dust samples for now
its <- readRDS("./data/physeq_objects/full_its_ps_raw.RDS") %>% 
  subset_samples(sample_type == "dust")


# LIGHT CLEANING ####

# agglomerate ASVs at species level
# save ASV resolution for later if needed
ssu_asv <- ssu
its_asv <- its
rank_names(ssu);rank_names(its)
ssu <- 
  ssu %>% tax_glom("Species")
its <- 
  its %>% tax_glom("Species")

# remove non-fungi
ssu <- 
  ssu %>% 
  subset_taxa(Kingdom == "k__Fungi")
its <- 
  its %>% 
  subset_taxa(Kingdom == "k__Fungi")

# remove empty samples and taxa
ssu <- ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0)
its <- its %>% 
  subset_taxa(taxa_sums(its) > 0)
ssu <- ssu %>% 
  subset_samples(sample_sums(ssu) > 0)
its <- its %>% 
  subset_samples(sample_sums(its) > 0)

# build full physeq object (ASVs grouped by species assignment)
full <- merge_phyloseq(ssu_asv,its_asv) %>% 
  tax_glom("Species",bad_empty = c("k__unclassified",
                                   "p__unclassified",
                                   "c__unclassified",
                                   "o__unclassified",
                                   "f__unclassified",
                                   "g__unclassified",
                                   "s__unclassified",
                                   NA, "", " ", "\t")) %>% 
  subset_taxa(Kingdom == "k__Fungi")


# ALPHA DIVERSITY ####

# estimate Shannon div and richness
ssu_alpha <- estimate_richness(ssu,measures = c("Shannon","Observed"))
its_alpha <- estimate_richness(its,measures = c("Shannon","Observed"))

# add to data frame for plotting
ssu_alpha_df <- 
  microbiome::meta(ssu) %>% 
  mutate(shannon = ssu_alpha$Shannon,
         richness = ssu_alpha$Observed)
# compare previous and new richness measures
plot(ssu_alpha_df$ssu_asv_richness,ssu_alpha_df$richness)

its_alpha_df <- 
  microbiome::meta(its) %>% 
  mutate(shannon = its_alpha$Shannon,
         richness = its_alpha$Observed)
# compare previous and new richness measures
plot(its_alpha_df$its_asv_richness,its_alpha_df$richness)

alpha_df <- full_join(ssu_alpha_df,its_alpha_df)

richness_v_weather_plots <- 
alpha_df %>% 
  pivot_longer(contains("_jun_nov_"),
               names_to = "jun_nov_variable",
               values_to = "weather_value") %>% 
  ggplot(aes(x=weather_value,y=richness)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~amplicon*jun_nov_variable,
             scales='free')
shannon_v_weather_plots <- 
alpha_df %>% 
  pivot_longer(contains("_jun_nov_"),
               names_to = "jun_nov_variable",
               values_to = "weather_value") %>% 
  ggplot(aes(x=weather_value,y=shannon)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~amplicon*jun_nov_variable,
             scales='free')

## taxonomy plots ####

# merge by site, then transform counts to relabund
full_site <- 
full %>% 
  merge_samples("site",fun=sum) %>%
  transform_sample_counts(ra)

# arrange sites by latitude

full_site@sam_data$site <- factor(sample_names(full_site),
                                  levels = data.frame(site = sample_names(full_site),lat = full_site@sam_data$lat_dd) %>% 
                                    arrange(lat) %>% 
                                    pluck('site')
                                  )

full_site@tax_table[,2] <- full_site@tax_table[,2] %>% str_remove("p__")
full_site %>% 
  plot_bar2(x="site",fill = "Phylum") +
  labs(x="Site",y="Relative abundance") +
  scale_fill_viridis_d(begin=0,end=1,option='turbo')


# GUILDS ####

# vector of species names to match fungaltraits database
species_names <- 
paste0(
  full@tax_table[,6] %>% str_remove("g__"),
  "_",
  full@tax_table[,7] %>% str_remove("s__")
)
species_names[duplicated(species_names)]

# get fungaltraits database
fungaltraits_db <- fungal_traits()

guilds <- list()
for(spp in species_names){
  rows <- fungaltraits_db[fungaltraits_db$species == spp,]
  if(nrow(rows) < 1){guilds[spp] <- NA;next}
  if(nrow(rows) > 0){
    guilds[spp] <- unique(rows[["guild_fg"]])
  }
}

guild_df <- data.frame(genus=species_names %>% str_split("_") %>% map_chr(1),
                       species = species_names,
                       fungaltraits_guild = unlist(guilds))

# get funguild database
funguild <- FUNGuildR::get_funguild_db()

# try again to compare assignments at species taxonomy level
guilds <- list()
for(spp in species_names){
  rows <- funguild[funguild$taxon == str_replace(spp,"_"," "),]
  if(nrow(rows) < 1){guilds[spp] <- NA;next}
  if(nrow(rows) > 0){
    guilds[spp] <- unique(rows[["guild"]])
  }
}
guild_df$funguild_guild <- unlist(guilds)

# try again to compare assignments at genus taxonomy level

# find rows where there is no assignment yet
part1 <- is.na(guild_df$fungaltraits_guild)
part2 <- is.na(guild_df$funguild_guild)
both <- part1 & part2
guild_df$need_guild <- both



guild_df <- 
guild_df %>% 
  mutate(funguild_genus_guild = FUNGuildR::funguild_assign(guild_df,db=funguild,tax_col = "genus") %>% 
           pluck('guild')) %>% 
  mutate(guild_assignment = case_when(need_guild == TRUE ~ funguild_genus_guild,
                                      need_guild != TRUE ~ funguild_guild)) %>% 
  dplyr::select(species,guild_assignment)

guild_df$mycorrhizal <- grepl("Mycorrhizal",ignore.case = TRUE,guild_df$guild_assignment)

#build full taxonomy with guild assignments
taxonomy <- 
tax_table(full) %>% 
  as.data.frame %>% 
  cbind(guild_df) %>% 
  mutate(across(everything(),function(x){str_remove(x,".__")}))
# add guild to phyloseq object
full@tax_table[,1] <- taxonomy$guild_assignment
colnames(full@tax_table)[1] <- "Guild"

## Can subset to mycorrhizal-only ####
mycorrhizal <- 
  full %>% 
  subset_taxa(grepl("Mycorrhizal",ignore.case = TRUE,full@tax_table[,1]))
  
mycorrhizal_ra <- 
  full %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(grepl("Mycorrhizal",ignore.case = TRUE,full@tax_table[,1]))

mycorrhizal@tax_table[,1] <- ifelse(grepl("Arbusc",mycorrhizal@tax_table[,1],ignore.case = TRUE),
                                    "Arbuscular","Ecto")


sam.order <- mycorrhizal %>%
  tax_glom("Guild") %>% 
  merge_samples('site') %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  dplyr::select(sample=Sample,precip=total_precip_jun_nov_mm) %>% 
  arrange(precip) %>% 
  unique.data.frame() %>% 
  pluck('sample')

mycorrhizal %>%
  tax_glom("Guild") %>% 
  merge_samples('site') %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  mutate(Sample = factor(Sample,levels=sam.order)) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Guild)) +
  geom_col() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        axis.text.y = element_text(face='bold',size=12),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5)) +
  scale_fill_manual(values = pal.discrete[c(1,20)])
  
ssu_melted_raw <- 
  ssu %>% 
  psmelt()
its_melt_raw <- 
  its %>% 
  psmelt()




# DISTANCE MEASURES ####
ssu_comm_dist <- 
ssu %>% 
  transform_sample_counts(ra) %>% 
  otu_table() %>% 
  vegdist()
its_comm_dist <- 
  its %>% 
  transform_sample_counts(ra) %>% 
  otu_table() %>% 
  vegdist()

# using distHaversine (meters separation)

# build geographic data frames and subset to complete cases
ssu_latlon <- data.frame(lon=ssu_alpha_df$long_dd,
                         lat=ssu_alpha_df$lat_dd)
ssu_comm_dist <- as.matrix(ssu_comm_dist)[complete.cases(ssu_latlon),complete.cases(ssu_latlon)]
ssu_latlon <- ssu_latlon[complete.cases(ssu_latlon),]

its_latlon <- data.frame(lon=its_alpha_df$long_dd,
                         lat=its_alpha_df$lat_dd)
its_comm_dist <- as.matrix(its_comm_dist)[complete.cases(its_latlon),complete.cases(its_latlon)]
its_latlon <- its_latlon[complete.cases(its_latlon),]

ssu_geog_dist <- ssu_latlon %>% 
  geosphere::distm()
its_geog_dist <- its_latlon %>% 
  geosphere::distm()

dist_df <- 
data.frame(geog_dist_m=c(ssu_geog_dist),
           comm_dist=c(ssu_comm_dist),
           region="SSU") %>% 
  bind_rows(
    data.frame(geog_dist_m=c(its_geog_dist),
               comm_dist=c(its_comm_dist),
               region="ITS")
  )
## plot geog vs community distances
dist_df %>% 
  ggplot(aes(x=geog_dist_m/1000,y=comm_dist,color=region)) +
  geom_point(color='black',alpha=.05) +
  geom_smooth(se=FALSE,linewidth=3) +
  labs(x="Geographic distance (km)",y="Community distance (Bray)",color="Amplicon") +
  scale_color_viridis_d(option = 'turbo',end=.8,begin = .2)


# # multiple regresison on matrices
ssu_mrm <- ecodist::MRM(as.dist(ssu_comm_dist) ~ as.dist(ssu_geog_dist))
its_mrm <- ecodist::MRM(as.dist(its_comm_dist) ~ as.dist(its_geog_dist))
ssu_mrm 
its_mrm 



# ORDINATIONS ####
ssu_ord <- ssu %>% 
  ordinate(method="PCoA")
its_ord_rda <- its %>% 
  ordinate(method="RDA")
its_ord_dca <- its %>% 
  ordinate(method = "DCA")
its_ord_pcoa <- its %>% 
  ordinate(method = "PCoA")


# exploratory plots
alpha_df %>% names
# canister height
ssu %>% 
  plot_ordination(ssu_ord,color="canister_height")
its %>% 
  plot_ordination(its_ord_pcoa,color="canister_height")
# jun-nov total precip
ssu %>% 
  plot_ordination(ssu_ord,color="total_precip_jun_nov_mm")
its %>% 
  plot_ordination(its_ord_pcoa,color="total_precip_jun_nov_mm")

# jun-nov mean wind speed
ssu %>% 
  plot_ordination(ssu_ord,color="mean10m_wind_jun_nov_m_s")
its %>% 
  plot_ordination(its_ord_pcoa,color="mean10m_wind_jun_nov_m_s")

# site and height
ssu %>% 
  plot_ordination(ssu_ord,color="site") +
  facet_wrap(~canister_height)
its %>% 
  plot_ordination(its_ord_pcoa,color="site") +
  facet_wrap(~canister_height)

its %>% 
  plot_ordination(its_ord_pcoa,color="site")


# why do some of these have so much missing metadata?
its@sam_data[which(is.na(its@sam_data$site)),]
# ^^^ These don't even have the NEON site listed. Go back and check metadata import/cleaning
its %>% 
  plot_ordination(its_ord_pcoa,color="site_avg_elevation_m")



# PERMANOVAS ####

# this function performs a permanova using a vector of predictor names
# it automatically handles missing values and subsets accordingly
auto_permanova(physeq = ssu,
                data = ssu_alpha_df,
                pred.cols = c("mean10m_wind_jun_nov_m_s",
                              "total_precip_jun_nov_mm",
                              "site_avg_elevation_m"),
                strata.col = "site",
                mod.type = "additive")

auto_permanova(physeq = ssu,
                data = ssu_alpha_df,
                pred.cols = c("mean10m_wind_jun_nov_m_s",
                              "total_precip_jun_nov_mm",
                              "site_avg_elevation_m"),
                strata.col = "site",
                mod.type = "interactive")

auto_permanova(physeq = its,
                data = its_alpha_df,
                pred.cols = c("mean10m_wind_jun_nov_m_s",
                              "total_precip_jun_nov_mm",
                              "site_avg_elevation_m"),
                strata.col = "site",
                mod.type = "additive")

# GDM ####

# extract species by site info

# need tax table to do this
ssu_melt <- 
  ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0) %>% 
  transform_sample_counts(ra) %>% 
  psmelt()
its_melt <- 
  its %>% 
  subset_taxa(taxa_sums(its) > 0) %>% 
  transform_sample_counts(ra) %>% 
  psmelt()


sppTab_ssu <-  ssu_melt %>% dplyr::select(OTU,library_id,long_dd,lat_dd,Abundance)
sppTab_its <- its_melt %>% dplyr::select(OTU,library_id,long_dd,lat_dd,Abundance)

predictors <- c("height_cm","total_precip_jun_nov_mm",
                "mean10m_wind_jun_nov_m_s","collection_year",
                "min_temp_jun_nov_c","max_temp_jun_nov_c")
# get columns with site ID, env. data, and xy-coordinates
envTab_ssu <- ssu_melt %>% dplyr::select(library_id, long_dd, lat_dd, all_of(predictors))
envTab_its <- its_melt %>% dplyr::select(library_id, long_dd, lat_dd, all_of(predictors))

# remove rows with missing values
good_ssu <- complete.cases(envTab_ssu)
good_its <- complete.cases(envTab_its)
envTab_ssu <- envTab_ssu[good_ssu,]
envTab_its <- envTab_its[good_its,]
sppTab_ssu <- sppTab_ssu[good_ssu,]
sppTab_its <- sppTab_its[good_its,]

# format for gdm
gdmTab_ssu <- formatsitepair(bioData=sppTab_ssu, 
                             bioFormat=2, #x-y spp list
                             XColumn="long_dd", 
                             YColumn="lat_dd",
                             sppColumn="OTU", 
                             siteColumn="library_id", 
                             predData=envTab_ssu,
                             abundance = TRUE,
                             abundColumn = "Abundance")

gdmTab_its <- formatsitepair(bioData=sppTab_its, 
                             bioFormat=2, #x-y spp list
                             XColumn="long_dd", 
                             YColumn="lat_dd",
                             sppColumn="OTU", 
                             siteColumn="library_id", 
                             predData=envTab_its,
                             abundance = TRUE,
                             abundColumn = "Abundance")

# fit GDM models
gdm_ssu <- gdm(data = gdmTab_ssu,geo = TRUE)
gdm_its <- gdm(data = gdmTab_its,geo = TRUE)
# quick look at model fits
summary(gdm_ssu)
summary(gdm_its)

# predictions from model (using same distances)
gdm_ssu_pred <- predict(object=gdm_ssu, data=gdmTab_ssu)
gdm_its_pred <- predict(object=gdm_its, data=gdmTab_its)

ssu_preds <- data.frame(observed = gdmTab_ssu$distance,
                        predicted = gdm_ssu_pred,
                        dist_from_edge = gdm_ssu$ecological,
                        sample_type = "SSU")
its_preds <- data.frame(observed = gdmTab_its$distance,
                         predicted = gdm_its_pred,
                         dist_from_edge = gdm_its$ecological,
                         sample_type = "ITS")

full_join(ssu_preds,its_preds) %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~sample_type,scales = 'free')

ssu_splines <- gdm::isplineExtract(gdm_ssu) %>% as.data.frame() %>% mutate(amplicon='SSU')
its_splines <- gdm::isplineExtract(gdm_its) %>% as.data.frame() %>% mutate(amplicon='ITS')

isplines <- full_join(ssu_splines,its_splines)
names(isplines) <- 
  names(isplines) %>% 
  str_replace("x\\.","actual_") %>% 
  str_replace("y\\.","partial_")


p1 <- isplines %>% 
  ggplot(aes(x=actual_Geographic,y=partial_Geographic,color=amplicon)) +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none',
        axis.title = element_text(face='bold',size=16),
        axis.text.y = element_text(face='bold',size=12),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5)) +
  scale_color_viridis_d(begin = .1,end=.8,option='turbo') +
  labs(y="f(Geog. distance)",x="Geo. distance",color="Amplicon")
p1

p2 <- isplines %>% 
  ggplot(aes(x=actual_total_precip_jun_nov_mm,y=partial_total_precip_jun_nov_mm,color=amplicon)) +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        axis.text.y = element_text(face='bold',size=12),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5)) +
  scale_color_viridis_d(begin = .1,end=.8,option='turbo') +
  labs(y="f(Precip. (Jun-Nov))",x="Precip. (Jun-Nov)",color="Amplicon")
p2

p3 <- isplines %>% 
  ggplot(aes(x=actual_mean10m_wind_jun_nov_m_s,y=partial_mean10m_wind_jun_nov_m_s,color=amplicon)) +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none',
        axis.title = element_text(face='bold',size=16),
        axis.text.y = element_text(face='bold',size=12),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5)) +
  scale_color_viridis_d(begin = .1,end=.8,option='turbo') +
  labs(y="f(Wind - 10m (Jun-Nov))",x="Wind - 10m (Jun-Nov)",color="Amplicon")
p3

p4 <- isplines %>% 
  ggplot(aes(x=actual_height_cm,y=partial_height_cm,color=amplicon)) +
  geom_smooth(se=FALSE,linewidth=2) +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none',
        axis.title = element_text(face='bold',size=16),
        axis.text.y = element_text(face='bold',size=12),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5)) +
  scale_color_viridis_d(begin = .1,end=.8,option='turbo') +
  labs(y="f(Canister height (cm))",x="Canister height (cm)",color="Amplicon")
p4


(p1 | p2) / (p3 | p4) & patchwork::plot_layout(guides = 'collect')
