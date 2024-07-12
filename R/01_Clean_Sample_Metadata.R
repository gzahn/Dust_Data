#
# note: redo this, checking for variable significance,
#       once all data have been updated
#

# SETUP ####

# load packages
library(tidyverse)
library(janitor)
library(skimr)
library(ranger)
library(vip)
library(broom)

# handy shortcut functions
'%ni%' <- Negate('%in%')

# load metadata that was downloaded from Google Drive
for_analysis <- readRDS("./data/dust_data_for_analysis.RDS")
collection_data <- readRDS("./data/collection_data_responses.RDS")
lib_data <- readRDS("./data/library_data.RDS")

# placeholder for library prep data frame

# EXPLORE AND CLEAN ####

## "for analysis" data frame ####
glimpse(for_analysis)
skimr::skim(for_analysis)

# list of column names that are character that should be numeric
num_cols <- c("height","collection_year","nano_ng_u_l",
  "a260_a280","a260_a230","height_cm","lat_dd",
  "mwac_condition","gelband_pa",
  "its_reads","its_asv_richness","ssu_reads","ssu_asv_richness",
  "site_avg_elevation_m","mean_temperature_o_c","mean_precipitation_mm",
  "x30yr_percip","x30yr_temp")

# list of "notes" columns
notes_cols <- c("water_pa","plantmatter_pa","insect_pa","insect","anomoly","lab_error")

# do class conversions
for_analysis <- 
for_analysis %>% 
  mutate(long_dd = long_dd %>% # convert longitude from list to numeric
           map_chr(1) %>% # don't know why this is a 'list' but it was
           as.numeric()) %>% 
  mutate(across(all_of(num_cols), as.numeric)) %>% 
  mutate(date_nanodrop = date_nanodrop %>% 
           as.POSIXct(format = '%m/%d/%Y, %I:%M:%S %p'))

# find unique notes
notes_values <- for_analysis %>% 
  dplyr::select(all_of(notes_cols)) %>% 
  lapply(unique)

# unique insect names from column
insect_options <- notes_values$insect %>% 
  str_split(",") %>% 
  map(str_squish) %>% 
  unlist %>% 
  unique()
# find logicals for each insect name
# first, remove the "no" insect from the options
insect_options <- insect_options[insect_options != "no"]  

insect_lgl_columns <- list()
for(i in insect_options){
  x <- grepl(pattern = i,for_analysis$insect)
  insect_lgl_columns[[i]] <- x
}

insect_columns <- insect_lgl_columns %>% as.data.frame() %>% janitor::clean_names()
# add to main data frame
for_analysis <- bind_cols(for_analysis,insect_columns)
# pull names of insect column for use in RF predictions
insect_columns <- names(insect_columns)

# change insect_pa to logical
for_analysis$insect_pa <- if_else(for_analysis$insect_pa == "yes",TRUE,FALSE)

# quick random forest to see which insect columns are predictive of richness, if any
# make new clean data frames for ssu and its (independently)
ssu <- for_analysis %>% 
  dplyr::select(ssu_asv_richness,all_of(insect_columns))
ssu <- ssu[complete.cases(ssu),]
its <- for_analysis %>% 
  dplyr::select(its_asv_richness,all_of(insect_columns))
its <- its[complete.cases(its),]
# run ranger model to see if any insect columns are predictive of richness
ssu.mod <- ranger::ranger(formula = ssu_asv_richness ~ .,
               data=ssu,importance = "permutation")
vip::vip(ssu.mod)
ssu.mod$r.squared
its.mod <- ranger::ranger(formula = its_asv_richness ~ .,
               data=its,importance = "permutation")
vip::vip(its.mod)
its.mod$r.squared
# those are some tiny r squared values...these insects aren't explaining much for richness

# try a full glm
mod.insect.its <- 
its %>% 
  pivot_longer(-its_asv_richness, names_to = "insect") %>% 
  glm(its_asv_richness ~ insect * value, data = .)

mod.insect.ssu <- 
  ssu %>% 
  pivot_longer(-ssu_asv_richness, names_to = "insect") %>% 
  glm(ssu_asv_richness ~ insect * value, data = .)

# show only significant terms (looking at whether ANY insect is predictive of its richness)
broom::tidy(mod.insect.its) %>% 
  dplyr::filter(p.value < 0.06)
broom::tidy(mod.insect.ssu) %>% 
  dplyr::filter(p.value < 0.06)

# Note: both RF and GLM show no significant impact of insect presence in spore trap on 
# either ITS or SSU fungal richness.
# Can remove those columns...

# remove unimportant columns from "for_analysis" data frame
# If we decide to leave these in, skip this
for_analysis <- 
  for_analysis %>% 
  dplyr::select(-all_of(c(insect_columns,"insect")))

# fix misspelled column name - domian -> domain
names(for_analysis)[names(for_analysis) == "domian"] <- "domain"
glimpse(for_analysis)

# check with Paul what these mean and if they're good to go...
for_analysis$am_em_dom %>% table

# convert more lgl columns to T/F
lgl_cols <- c("water_pa","plantmatter_pa","lab_error")
for_analysis[,lgl_cols] %>% lapply(unique) # all consistent no/yes values
for_analysis <- 
  for_analysis %>% 
  mutate(across(all_of(lgl_cols),function(x){if_else(x == "no",FALSE,TRUE)}))

# minor adjustment to "anomoly" column
for_analysis <- 
  for_analysis %>% 
  mutate(anomoly = case_when(anomoly == "no" ~ "none",
                             anomoly == "slowToDrain" ~ "slow_to_drain",
                             anomoly == "lotsOfSoil" ~ "lots_of_soil"),
         anomoly = factor(anomoly,levels = c('none','lots_of_soil','slow_to_drain')))

# is anomoly predictive of richness? can we remove it?
for_analysis %>% 
  glm(its_asv_richness ~ anomoly,data=.) %>%
  summary()
for_analysis %>% 
  glm(ssu_asv_richness ~ anomoly,data=.) %>%
  summary()
# can safely remove it...no predictive power for its or ssu richness
for_analysis <- 
  for_analysis %>% 
  dplyr::select(-anomoly)
for_analysis$collection_year %>% table
collection_data$year %>% table

# test predictive ability of remaining lgl cols
its.plantmatter <- 
for_analysis %>% 
  glm(its_asv_richness ~ plantmatter_pa, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)
its.water <- 
for_analysis %>% 
  glm(its_asv_richness ~ water_pa, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)
its.laberror <- 
  for_analysis %>% 
  glm(its_asv_richness ~ lab_error, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)
ssu.plantmatter <- 
  for_analysis %>% 
  glm(ssu_asv_richness ~ plantmatter_pa, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)
ssu.water <- 
  for_analysis %>% 
  glm(ssu_asv_richness ~ water_pa, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)
ssu.laberror <- 
  for_analysis %>% 
  glm(ssu_asv_richness ~ lab_error, data = .) %>% 
  broom::tidy() %>% 
  dplyr::filter(p.value < 0.05)

lgl_predictive_capacity <- 
(its.plantmatter %>% 
  full_join(its.water) %>% 
  full_join(its.laberror) %>% 
  mutate(region = "ITS") %>% 
  dplyr::filter(term != "(Intercept)")) %>% 
  bind_rows(
    (ssu.plantmatter %>% 
       full_join(ssu.water) %>% 
       full_join(ssu.laberror) %>% 
       mutate(region = "SSU") %>% 
       dplyr::filter(term != "(Intercept)"))
  )
lgl_predictive_capacity
# seems like plant matter and lab error correlate with ITS (but not SSU) richness
# plot them...
for_analysis$plantmatter_pa %>% table


for_analysis %>%
  dplyr::select(its_asv_richness,plantmatter_pa,lab_error) %>% 
  pivot_longer(all_of(c("plantmatter_pa","lab_error")),
               names_to = "condition") %>% 
  ggplot(aes(x=value,y=its_asv_richness)) +
  geom_boxplot() +
  facet_wrap(~condition) +
  theme_bw()
# leave in those columns, and remove "water_pa" and "insect_pa" (since it was already dealt with)
for_analysis$water_pa <- NULL
for_analysis$insect_pa <- NULL

## "collection data" data frame ####
collection_data %>% glimpse

# build "index" column to match up with "for_analysis" data frame
# 001-1-21 : dust_collector_label-canister_height-year(2-digit)

# first, clean up each column
collection_data$dust_collector_label %>% unique # some entries left off leading '0'
# collection_data$dust_collector_label <- # undo this 
# collection_data$dust_collector_label %>% 
#   str_replace_all("^01$","001") %>% 
#   str_replace_all("^02$","002") 
for_analysis$index
collection_data$canister_height %>% unique # imported as a list, some entries have label appended
collection_data$canister_height <- # just pull last digit
collection_data$canister_height %>% 
  unlist %>% 
  str_sub(start = -1)

collection_data$year %>% unique # looks good

# combine to build "index"
collection_data <- 
  collection_data %>% 
  mutate(index = paste0(dust_collector_label,'-',canister_height,'-',str_sub(year,start=-2)))

# now, we can merge with for_analysis based on that column
# but what is the overlap between the two?
print(paste(collection_data$index %in% for_analysis$index %>% sum,
            "samples from collection_data are also present in for_analysis."))
# pretty close, but where are the other collection_data samples coming from and why are
# they not in for_analysis?  (Ask Paul)

# for now, go ahead and subset to just those that are in "for_analysis"
collection_data_subset <- 
  collection_data %>% 
  dplyr::filter(collection_data$index %in% for_analysis$index)

# what are the "notes_" columns?
collection_data_subset %>% 
  dplyr::select(starts_with("notes_")) %>% 
  lapply(unique)
# okay, they're just a bunch of notes ;)
# at least we can combine them together into a single column
notes_cols <- 
collection_data_subset %>% 
  dplyr::select(starts_with("notes_")) %>% names
collection_data_subset <- 
  collection_data_subset %>% 
  mutate(notes = paste(notes_9,notes_11,notes_13,notes_15,notes_17,notes_19,notes_21,notes_23,notes_25,sep=';')) %>% 
  dplyr::select(-all_of(notes_cols))
collection_data_subset$notes

# further cleanup
collection_data_subset <- 
  collection_data_subset %>% 
  mutate(notes = notes %>% 
           str_remove_all("NA;") %>% 
           str_remove_all(";NA"),
         notes = if_else(notes == "NA",NA,notes)) 

collection_data_subset$date %>% 
  unlist %>% 
  unique # ...oh for fuck's sake!
# not in excel date format either :(

# gotta do this manually
collection_data_subset$date %>% unique

collection_data_subset$date <- 
collection_data_subset$date %>% 
  map(function(x){x %>% as.character %>% str_remove(" UTC")}) %>% 
  map(as.character) %>% 
  unlist() %>% 
  str_replace("11/10/21","2021-11-10") %>% # one weird date version doesn't match others .. recheck w new data
  as.POSIXct(format = '%Y-%m-%d',tz = "UTC")

# clean up "position"
collection_data_subset$position_e_g_ne_corner_s_side %>% unique
# nope...these are all over the place. Just leave it for now.
collection_data_subset$position_e_g_ne_corner_s_side %>% table %>% as.data.frame()
## "combined" data sets ####

# full join with main data frame
full <- full_join(for_analysis,collection_data_subset,by="index") %>% 
  dplyr::select(-all_of(c("x26","x27")))
  
# clean up a couple column names...
names(full)[names(full) == "x30yr_percip"] <- "precip_30yr"
names(full)[names(full) == "x30yr_temp"] <- "temp_30yr"
# well, those columns are empty right now anyway, so delete them
full$precip_30yr <- NULL
full$temp_30yr <- NULL

# figure out why row numbers differ
full$index %ni% for_analysis$index

duplicated_samples <- 
full %>% 
  dplyr::filter(index %in% full$index[which(duplicated(full$index))])
duplicated_samples



# See if any of the logical note columns from the second spreadsheet have any impact on richness
new_lgl_cols <- 
c("is_the_collector_bent_over_or_not_positioned_straight",
"is_the_vegetation_preventing_the_mast_from_spinning_freely",
"otherwise_can_the_mast_spin_freely",
"does_the_collector_look_damaged_in_any_other_way",
"are_the_jars_tight_within_the_hose_clamps_are_the_hose_clamps_intact",
"are_the_antler_lids_all_attached",
"are_the_antler_tubes_postitioned_parallel_to_the_ground",
"is_the_longer_tube_facing_away_from_the_fin",
"do_the_antler_tubes_look_clear_of_insects_debris")

its_rf <- full[,c("its_asv_richness",new_lgl_cols)]
its_rf <- its_rf[complete.cases(its_rf),]
mod.rf.its <- ranger::ranger(formula = its_asv_richness ~ .,data=its_rf,importance = "permutation")
mod.rf.its$r.squared
vi(mod.rf.its)[1,] %>% pluck("Variable")
glm(data = its_rf,
    formula = its_asv_richness ~ does_the_collector_look_damaged_in_any_other_way) %>% 
  summary

ssu_rf <- full[,c("ssu_asv_richness",new_lgl_cols)]
ssu_rf <- ssu_rf[complete.cases(ssu_rf),]
mod.rf.ssu <- ranger::ranger(formula = ssu_asv_richness ~ .,data=ssu_rf,importance = "permutation")
mod.rf.ssu$r.squared
vi(mod.rf.ssu)[1,] %>% pluck("Variable")
glm(data = ssu_rf,
    formula = ssu_asv_richness ~ otherwise_can_the_mast_spin_freely) %>% 
  summary

full$otherwise_can_the_mast_spin_freely %>% table


# potentially important notes columns (for richness in its and ssu)
# "does_the_collector_look_damaged_in_any_other_way",
# "otherwise_can_the_mast_spin_freely"

# remove other vars
full <- 
full %>%
  dplyr::select(-all_of(new_lgl_cols[new_lgl_cols %ni% c("does_the_collector_look_damaged_in_any_other_way",
                                                         "otherwise_can_the_mast_spin_freely")]))

# clean up random list column
# distance_from_ground_to_center_of_top_canister_cm
full$distance_from_ground_to_center_of_top_canister_cm %>% 
  map(as.numeric) %>% unlist

# convert NULL values to NA
full$distance_from_ground_to_center_of_top_canister_cm[full$distance_from_ground_to_center_of_top_canister_cm %>% 
                                                         map(is.null) %>% 
                                                          unlist] <- NA
# now convert to numeric
full$distance_from_ground_to_center_of_top_canister_cm <- 
  full$distance_from_ground_to_center_of_top_canister_cm %>% unlist %>% as.numeric


## get rid of SCBI2 (alleged forest scary bear) samples ####
full <- 
  full %>% dplyr::filter(site != "SCBI2")


## add library info ####

# create "index" to match data sheet
# pull run #
# pull dust vs soil sample type
# pull gel bands 
# check all worksheets in those files (each row in each tab is separate sample)




## add filepath info ####
# Need to be able to link this sample metadata to the ASV table


## add NEON site envir. data ####
# local temp
# soil moisture
# wind direction/speed
# ...



# EXPORT ####

saveRDS(full,"./data/full_clean_metadata.RDS")

glimpse(full)
names(full)
