# SETUP ####

# load packages
library(tidyverse)
library(googlesheets4)
library(googledrive)
library(janitor)


# authenticate Google drive
# use browser option ... sign in to Dartmouth account in browser, giving full permissions to API
googledrive::drive_auth() # select option 1, use browser to complete

# google drive URLs
dust_data_for_analysis <- "https://docs.google.com/spreadsheets/d/18weKsvAKS9TMZeKMTzM85gti8NYsPXcxd33pgnI_xkw/edit?gid=1423804161#gid=1423804161"
collection_data_responses <- "https://docs.google.com/spreadsheets/d/1XDKWhWwQIUKnYm699nhuVxhNXHJW7qVrcqQ_HqIpwho/edit?gid=207128096#gid=207128096"
library_data <- "https://docs.google.com/spreadsheets/d/15y2TKOAPEYAAUtAy0jAHdi0q45WhEjTnWWGwvnlGE_U/edit?gid=0#gid=0"


# LOAD DATA FROM GOOGLE DRIVE ####

# dust data for analysis (main sheet)
for_analysis <- googledrive::drive_get(dust_data_for_analysis) %>% 
  googlesheets4::read_sheet(sheet = "ForAnalysis") %>% 
  janitor::clean_names() # automatically clean the names of columns

# dust data for analysis (NEON Weather station data by site)
wind <- googledrive::drive_get(dust_data_for_analysis) %>% 
  googlesheets4::read_sheet(sheet = "tristanWind") %>% 
  janitor::clean_names() # automatically clean the names of columns

temp <- googledrive::drive_get(dust_data_for_analysis) %>% 
  googlesheets4::read_sheet(sheet = "tristanTemp") %>% 
  janitor::clean_names() # automatically clean the names of columns

precip <- googledrive::drive_get(dust_data_for_analysis) %>% 
  googlesheets4::read_sheet(sheet = "tristanPrecip") %>% 
  janitor::clean_names() # automatically clean the names of columns



# collection data responses - all years (main sheet)
collection_data <- googledrive::drive_get(collection_data_responses) %>% 
  googlesheets4::read_sheet(sheet = "Edited and most complete copy of collection") %>% 
  janitor::clean_names()

# library prep sheets
lib_data <- googledrive::drive_get(library_data) %>% 
  googlesheets4::read_sheet(sheet = "all_data") %>% 
  janitor::clean_names()


# SAVE METADATA FILES ####
saveRDS(for_analysis, "./data/dust_data_for_analysis.RDS")
saveRDS(collection_data, "./data/collection_data_responses.RDS")
saveRDS(library_data, "./data/library_data.RDS")
saveRDS(wind,"./data/wind_data.RDS")
saveRDS(temp,"./data/temp_data.RDS")
saveRDS(precip,"./data/precip_data.RDS")

# export library prep data as well, but only as a single merged data frame


