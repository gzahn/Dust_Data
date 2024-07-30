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

# functions and themes
source("./R/functions.R")

# envir. variables
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

set.seed(666)

theme_set(theme_bw() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face='bold',size=16),
                  axis.title = element_text(face='bold',size=16),
                  axis.text = element_text(face='bold',size=12),
                  legend.title = element_text(face='bold',size=16),
                  legend.text = element_text(face='bold',size=12)))
readRDS("./data/ASV_Tables/Run_1_SSU_ASV_Table.RDS") %>%
  sum
readRDS("./data/ASV_Tables/Run_2_SSU_ASV_Table.RDS") %>%
  sum
readRDS("./data/ASV_Tables/Run_5_SSU_ASV_Table.RDS") %>%
  sum

# data
ssu <- readRDS("./data/physeq_objects/full_ssu_ps_raw.RDS")
its <- readRDS("./data/physeq_objects/full_its_ps_raw.RDS")
sample_sums(ssu) %>% sum

site_info <- read_csv("./data/siteinfo_updated.csv") %>% 
  arrange(AM_EM_dom)

# CLEANING ####

# clean up taxonomy names
ssu <- clean_ps_taxonomy(ssu)
its <- clean_ps_taxonomy(its)

# dump everything except glomeromycota from SSU
ssu <- 
  ssu %>% 
  subset_taxa(Phylum == "Glomeromycota")
sample_sums(ssu) %>% sum

# get rid of non-fungi from ITS, along with any Glomeromycota
its <- 
  its %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum) & Phylum != "Glomeromycota")

# get rid of control samples
its <- 
its %>% 
  subset_samples(sample_names(its) %in% grep("^CON|^POS|neg_",sample_names(its),invert = TRUE,value = TRUE))
ssu <- 
  ssu %>% 
  subset_samples(sample_names(ssu) %in% grep("^CON|^POS|neg_",sample_names(ssu),invert = TRUE,value = TRUE))



# dump empty samples/taxa
its <- 
its %>% 
  subset_samples(sample_sums(its) > 1)
ssu <- 
ssu %>% 
  subset_samples(sample_sums(ssu) > 1)
its <- 
  its %>% 
  subset_taxa(taxa_sums(its) > 0)
ssu <- 
  ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0)


# JOIN PHYSEQS ####
full_asv <- merge_phyloseq(ssu,its)


# MERGE ASVs BY TAXONOMY ####
full <- 
  full_asv %>% 
  tax_glom("Species",bad_empty = c("unclassified",NA, "", " ", "\t"),NArm = FALSE)

full@sam_data$run_id

# get rid of samples that (for now) aren't in the "for_analysis" sheet with metadata
full_asv <- 
  full_asv %>% 
  subset_samples(!is.na(site))
full <- 
  full %>% 
  subset_samples(!is.na(site))


# two versions of the physeq object:
# full_asv (asv resolution); full (species resolution)
full_asv
full


# remove empty samples and taxa
full <- full %>% 
  subset_taxa(taxa_sums(full) > 0)
full_asv <- full_asv %>% 
  subset_taxa(taxa_sums(full_asv) > 0)
full <- full %>% 
  subset_samples(sample_sums(full) > 0)
full_asv <- full_asv %>% 
  subset_samples(sample_sums(full_asv) > 0)

# bit of tidying
full@sam_data$total_precip_jun_nov_mm <- full@sam_data$total_precip_jun_nov_mm %>% unlist %>% as.numeric()
full_asv@sam_data$total_precip_jun_nov_mm <- full_asv@sam_data$total_precip_jun_nov_mm %>% unlist %>% as.numeric()


# ALPHA DIVERSITY ####

# estimate Shannon div and richness
full_alpha <- estimate_richness(full,measures = c("Shannon","Observed"))
full_asv_alpha <- estimate_richness(full_asv,measures = c("Shannon","Observed"))

# add to data frame for plotting
full_alpha_df <- 
  microbiome::meta(full) %>% 
  mutate(shannon = full_alpha$Shannon,
         richness = full_alpha$Observed)

full_asv_alpha_df <- 
  microbiome::meta(full_asv) %>% 
  mutate(shannon = full_asv_alpha$Shannon,
         richness = full_asv_alpha$Observed)

alpha_df <- full_join(full_alpha_df,full_asv_alpha_df)

## richness vs weather ####
richness_v_weather_plots <- 
alpha_df %>% 
  pivot_longer(contains("_jun_nov_"),
               names_to = "jun_nov_variable",
               values_to = "weather_value") %>% 
  ggplot(aes(x=weather_value,y=richness)) +
  geom_point(alpha=.25) +
  geom_smooth(aes(color=factor(year)),method='lm') +
  facet_wrap(~amplicon*jun_nov_variable,
             scales='free') +
  theme(strip.text = element_text(face='bold',size=10)) +
  scale_color_viridis_d(begin=.1,end = .9,option='turbo')
saveRDS(richness_v_weather_plots,"./output/figs/richness_v_weather_plots.RDS")


## taxonomy plots ####

# merge by site/year, then transform counts to relabund

full@sam_data$mergevar <- paste(full@sam_data$site,full@sam_data$year,sep="_")

# pull metadata for easy access
meta <- as(full@sam_data,"data.frame")
# find yearly site means for precip
precip_summary <- 
meta %>% 
  group_by(site,year) %>% 
  summarize(precip = mean(total_precip_jun_nov_mm,na.rm = TRUE))
# precip_summary %>% View
meta$total_precip_jun_nov_mm
full_site <- 
full %>% 
  merge_samples("mergevar",fun=sum) %>%
  transform_sample_counts(ra)
# repair important metadata
full_site@sam_data$site <- sample_names(full_site) %>% str_split("_") %>% map_chr(1)
full_site@sam_data$year <- sample_names(full_site) %>% str_split("_") %>% map_chr(2)
full_site@sam_data$mean_precip_jun_nov_mm <- precip_summary$precip


# arrange sites by precip
full_site@sam_data$site <- factor(sample_names(full_site),
                                  levels = data.frame(site = sample_names(full_site),precip = full_site@sam_data$mean_precip_jun_nov_mm) %>% 
                                    arrange(precip) %>% 
                                    pluck('site')
                                  )

phylum_by_site_precip <- 
full_site %>% 
  plot_bar2(x="site",fill = "Phylum") +
  labs(x="Site",y="Relative abundance") +
  scale_fill_viridis_d(begin=0,end=1,option='turbo') +
  facet_wrap(~year,scales='free') +
  labs(caption = "Within each year, sites arranged by increasing precipitation.")
saveRDS(phylum_by_site_precip,"./output/figs/phylum_by_site_precip_barplot.RDS")


# AMF species list ####

full %>% 
  subset_taxa(Phylum == "Glomeromycota") %>% 
  tax_table() %>% 
  as("matrix") %>% 
  as.data.frame() %>% 
  mutate(count=full %>% 
           subset_taxa(Phylum == "Glomeromycota") %>% 
           taxa_sums(),
         ASV = row.names(.)) %>% 
  write_csv("./output/glomeromycota_species_list.csv")



# GUILDS ####

# get funguild database
funguild <- FUNGuildR::get_funguild_db()

# get fungaltraits database
fungaltraits_db <- fungal_traits()


# get search queries
query_taxonomy <- 
data.frame(
  asv = taxa_names(full),
  Genus = full@tax_table[,6],
  Species = full@tax_table[,7]
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
if(all(taxa_names(full) == guild_assignments$asv)){
  # add guild to phyloseq object
  full@tax_table[,1] <- guild_assignments$funguild_guild
} else (print("Out of order...double-check!"))
colnames(full@tax_table)[1] <- "Guild"

#

# transform to relabund
full_ra <- 
  full %>% 
  transform_sample_counts(ra)

# melt
full_ra_melt <- 
  psmelt(full_ra)
full_melt <- 
  psmelt(full)


full_ra_melt$Guild %>% unique %>% 
  str_remove_all("\\|") %>% 
  str_split("-") %>% 
  unlist %>% 
  unique
# simplify guild assignments
full_ra_melt <- 
  simplify_fungal_guilds(full_ra_melt)
full_melt <- 
  simplify_fungal_guilds(full_melt)

# add simplified guild to physeq



x <- 
data.frame(taxon=full_melt$OTU,
           guild=full_melt$major_guild) %>% 
  unique.data.frame() %>% 
  full_join(data.frame(taxon=taxa_names(full)))
row.names(x) <- x$taxon

if(all(taxa_names(full) == x[taxa_names(full),'taxon'])){
  full@tax_table[,1] <- x[taxa_names(full),'guild']
} else { print("Not in same order. recheck!")}

full_guild <- 
full %>% 
  tax_glom("Guild")



ssu <- 
full %>% 
  subset_samples(amplicon == "SSU")

elev_levels <- 
ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0) %>% 
  sample_data() %>% 
  as("data.frame") %>% 
  select(site_avg_elevation_m,site) %>% 
  unique.data.frame() %>% 
  arrange(site_avg_elevation_m) %>% 
  pluck("site")

guild_by_site <- 
  full_guild %>%
  merge_samples("site")
elev_levels2 <- 
  guild_by_site %>% 
  subset_taxa(taxa_sums(ssu) > 0) %>% 
  sample_data() %>% 
  as("data.frame") %>% 
  select(site_avg_elevation_m) %>%
  mutate(site = row.names(.)) %>% 
  unique.data.frame() %>% 
  arrange(site_avg_elevation_m) %>% 
  pluck("site")


guild_by_site %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  mutate(Sample = factor(Sample,levels=elev_levels2)) %>% 
  dplyr::filter(Abundance > 0) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Guild)) +
  geom_col() + 
  scale_fill_manual(values = pal$pal.okabe) +
  labs(y="Relative abundance",x="NEON Site") +
  theme(axis.text.x = element_text(angle = 270,face='bold',hjust=0,vjust=.5))
ggsave("./output/figs/major_guilds_by_site_barplot.png",height = 6, width = 8,dpi=400)
saveRDS(guild_by_site,"./output/figs/major_guilds_by_site_barplot.RDS")

# MYCORRHIZAL GROUPS ####
full@tax_table[,1] %>% unique %>% unname
ECM <- 
full %>% 
  subset_taxa(Guild %in% c("Ectomycorrhizal"))
AM <- 
  full %>% 
  subset_taxa(Guild %in% c("Arbuscular mycorrhizal"))
ErM <- 
  full %>% 
  subset_taxa(Guild %in% c("Ericoid mycorrhizal"))


# measure richness and add to metadata
full@sam_data$ecm_richness <- 
  ECM %>% 
  estimate_richness(measures = "Observed") %>% 
  pluck("Observed")
full@sam_data$amf_richness <- 
  AM %>% 
  estimate_richness(measures = "Observed") %>% 
  pluck("Observed")
full@sam_data$erm_richness <- 
  ErM %>% 
  estimate_richness(measures = "Observed") %>% 
  pluck("Observed")


# add to melted data frame as well
full_melt$ecm_richness <- full %>% psmelt %>% pluck("ecm_richness")
full_melt$amf_richness <- full %>% psmelt %>% pluck("amf_richness")
full_melt$erm_richness <- full %>% psmelt %>% pluck("erm_richness")


full_melt %>% names
  group_by(site) %>% 
  summarize(ecm_richness = sum(ecm_richness))

  ## envir. plot ####  
  
# plot ecm richness against envir variables
mycor_df <- 
  full_melt %>% 
  pivot_longer(c(contains("jun_nov"),site_avg_elevation_m),
               names_to = "measure",
               values_to = "value") %>% 
  dplyr::select(site,measure,value,ecm_richness,amf_richness,erm_richness,collection_year) %>% 
  unique.data.frame() %>% 
  mutate(measure = measure %>% 
           str_remove("_jun_nov.*") %>% 
           str_replace_all("_"," ") %>% 
           str_replace("10m"," 10 m ") %>% 
           str_to_sentence()) %>% 
  left_join(site_info) # add Bala's site info (canopy, veg type)

mycor_df %>% 
  pivot_longer(contains("richness"),
               names_to="mycorrhizal_group",
               values_to = "richness") %>% 
  mutate(mycorrhizal_group = mycorrhizal_group %>% 
           str_remove("_richness") %>% 
           str_replace("amf","AM") %>% 
           str_replace("ecm","EcM") %>% 
           str_replace("erm","ErM")) %>% 
  ggplot(aes(x=value,y=richness,color=mycorrhizal_group)) +
  geom_point(alpha=.25,size=2) +
  geom_smooth(se=FALSE,method='lm',linewidth=2) +
  facet_wrap(~measure,scales='free') +
  labs(y="Ecto- and ericoid-\nmycorrhizal richness",
       x="Envir. variable value",
       color="Mycorrhizal\ngroup") +
  theme(strip.text = element_text(face='bold',size=12)) +
  scale_color_manual(values = pal$pal.okabe[c(1,2,6)])
ggsave("./output/figs/mycorrhizal_groups_vs_environ_variables.png",width = 10, height = 7)  


## barplot ####

# cluster sites by veg type
site_order_veg <- site_info$site

# build df for all types of mycorrhizal taxa
mycor_merged <- 
  merge_phyloseq(ECM,AM,ErM) %>%
  merge_samples('site')
mycor_merged@sam_data$site <- sample_names(mycor_merged)

mycor_merged <- 
  mycor_merged %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  dplyr::select(-canopy) %>% 
  left_join(site_info) 

# find first and last sites in each group...
am_sites <-   
site_info[site_info$AM_EM_dom == "AM",] %>% 
  pluck('site')
am_sites <- am_sites[c(1,length(am_sites))]

ecm_sites <-   
  site_info[site_info$AM_EM_dom == "EcM",] %>% 
  pluck('site')
ecm_sites <- ecm_sites[c(1,length(ecm_sites))]


mycor_merged %>% 
  mutate(site=factor(site,levels=site_order_veg)) %>% 
  ggplot(aes(x=site,y=Abundance,fill=Guild)) +
  geom_col() +
  geom_point(size=8,shape = 15,aes(y=1.05,color=AM_EM_dom)) +
  scale_fill_manual(values = pal$pal.okabe[c(1,2,6)]) +
  scale_color_manual(values = pal$pal.earthtones) +
  labs(x="Site",y="Relative abundance",fill="Mycorrhizal\ngroup",color="Vegetation\ndominance") +
  theme(axis.text.x = element_text(face='bold',size=14,angle=270,hjust=0,vjust=0.5)) +
  guides(fill = guide_legend(override.aes = list(color = NA)))


amf
full_melt[grep("VTX00113",full_melt$Species),"run_id"]
  
amf_by_site <- 
ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0) %>% 
  merge_samples("site",fun = "sum") %>% 
  transform_sample_counts(ra)
amf_by_site@sam_data$site <- factor(sample_names(amf_by_site),levels = elev_levels)

amf_by_site <- 
amf_by_site %>% 
  plot_bar2(fill="Family",x='site') +
  scale_fill_manual(values = pal$pal.okabe) +
  labs(y="Relative abundance",x="NEON Site")
amf_by_site
saveRDS(amf_by_site,"./output/figs/amf_family_relabund_by_site_barplot.RDS")
ggsave("./output/figs/amf_family_relabund_by_site_barplot.png",height = 6,width = 7,dpi=400) 


# SITE MAP ####
# map of sites, showing amf present

# re-set ssu
ssu <- 
  full %>% 
  subset_samples(amplicon == "SSU")
ssu <- 
  ssu %>% 
  subset_samples(sample_names(ssu) %in% grep("^CON|^POS|neg_",sample_names(ssu),invert = TRUE,value = TRUE))
ssu <- 
  ssu %>% 
  subset_samples(sample_sums(ssu) > 1)
ssu <- 
  ssu %>% 
  subset_taxa(taxa_sums(ssu) > 0)

x <- data.frame(lon = full@sam_data$long_dd,
           lat = full@sam_data$lat_dd,
           site = full@sam_data$site) %>% 
  unique.data.frame() 

x <- 
x %>% 
mutate(amf_present = lat %in% unique(ssu@sam_data$lat_dd) &
           lon %in% unique(ssu@sam_data$long_dd)) %>% 
  group_by(site) %>% 
  reframe(amf_present = unique(amf_present),
            lon=mean(lon),lat=mean(lat))
present <- 
  x %>% 
  group_by(site) %>% 
  summarize(TF = sum(amf_present)) %>% 
  dplyr::filter(TF > 0) %>% pluck('site')

df_for_amf_map <- 
x %>% 
  select(site,lon,lat) %>% 
  unique.data.frame() %>% 
  mutate(amf_present = case_when(site %in% present ~ TRUE,
                                 TRUE ~ FALSE))

# make individual pie charts for each site
unique_sites <- full@sam_data$site %>% unique
SITE <- unique_sites[1]

for(SITE in unique_sites){
  x <- full %>% 
    subset_taxa(Phylum != "Glomeromycota") %>% 
    subset_samples(amplicon == "ITS") %>% 
    subset_samples(site == SITE) %>% 
    merge_samples('site')
  x <- x %>% 
    transform_sample_counts(ra) %>%
    psmelt()
  x %>% 
    ggplot(aes(x=" ",y=Abundance,fill=Phylum)) +
    geom_bar(stat = 'identity') +
    coord_polar("y",start=0) +
    scale_fill_manual(values = pal$pal.okabe[-1]) +
    theme_void() +
    theme(legend.position = 'none')
  ggsave(paste0("./output/figs/grobs/",SITE,".png"),width = .75,height = .75)  
}

# add image paths to df

grobs <- list.files("./output/figs/grobs",full.names = TRUE)
df_for_amf_map$piechart <- grobs


# save the basic plot as an image
# use ggimage to add them to the map

# build map styling
mapstyle <- rjson::fromJSON(file = "./R/mapstyle1.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# Get googlemap
area <- 
  ggmap::get_googlemap(center = c(lon = df_for_amf_map$lon %>% mean(na.rm=TRUE), 
                                  lat = df_for_amf_map$lat %>% mean(na.rm=TRUE)),
                       zoom = 3,
                       scale = 2,
                       style=mapstyle)

ggmap::ggmap(area) +
  geom_image(data=df_for_amf_map,aes(image=piechart)) +
  geom_point(data=dplyr::filter(df_for_amf_map,amf_present),
             aes(x=lon,y=lat),color='black',size=1.5,alpha=1,shape=19) +
  geom_point(data=dplyr::filter(df_for_amf_map,amf_present),
             aes(x=lon,y=lat),color='black',size=2,alpha=1,shape=8)

  geom_text_s(data=df_for_amf_map,color.target = 'all',
              aes(x=long_dd,y=lat_dd,label=site),
              nudge_x = 1,
              nudge_y = sites$nudge_y,
              color='white',size=2,
              point.padding = 0.4)
ggsave("./output/figs/map_with_piecharts_phylum.png",height = 10,width = 10,dpi=400)



amf_ra <- 
full_ra_melt %>% 
  dplyr::filter(sample_type == "dust" & amplicon == "SSU" & Abundance > 0)
amf_ra %>% names

amf_ra %>% 
  ggplot(aes())


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
