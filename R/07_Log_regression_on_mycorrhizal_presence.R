library(tidyverse)
library(ranger)
library(vip)

# load data
am <- readRDS("./output/mycorrhizal_type_logical_df_am.RDS")
em <- readRDS("./output/mycorrhizal_type_logical_df_em.RDS")

# pick predictor columns
pred_vars <- 
  c("am","height_cm","lat_dd","long_dd","site","am_em_dom","site_avg_elevation_m","canopy","year","total_precip_mm","sample_type")
class(am$year)
# subset
am <- 
  am %>% 
  select(all_of(pred_vars)) %>% 
  mutate(year = case_when(year == 4023 ~ 2023,
                          year == 4021 ~ 2021,
                          year == 4022 ~ 2022,
                          TRUE ~ year))
pred_vars <- 
  c("ecm","height_cm","lat_dd","long_dd","site","am_em_dom","site_avg_elevation_m","year","total_precip_mm","sample_type")
ecm <- 
  em %>% 
  select(all_of(pred_vars)) %>% 
  mutate(year = case_when(year == 4023 ~ 2023,
                          year == 4021 ~ 2021,
                          year == 4022 ~ 2022,
                          TRUE ~ year))
pred_vars <- 
  c("erm","height_cm","lat_dd","long_dd","site","am_em_dom","site_avg_elevation_m","year","total_precip_mm","sample_type")
erm <- 
  em %>% 
  select(all_of(pred_vars)) %>% 
  mutate(year = case_when(year == 4023 ~ 2023,
                          year == 4021 ~ 2021,
                          year == 4022 ~ 2022,
                          TRUE ~ year))
pred_vars <- 
  c("om","height_cm","lat_dd","long_dd","site","am_em_dom","site_avg_elevation_m","year","total_precip_mm","sample_type")
om <- 
  em %>% 
  select(all_of(pred_vars)) %>% 
  mutate(year = case_when(year == 4023 ~ 2023,
                          year == 4021 ~ 2021,
                          year == 4022 ~ 2022,
                          TRUE ~ year))

# Random Forest models
# subset to complete cols
mod <- ranger(data=am[complete.cases(am),],formula = am ~ .,importance = "permutation")
vip(mod)
mod <- ranger(data=ecm[complete.cases(ecm),],formula = ecm ~ .,importance = "permutation")
vip(mod)
mod <- ranger(data=erm[complete.cases(erm),],formula = erm ~ .,importance = "permutation")
vip(mod)
mod <- ranger(data=om[complete.cases(om),],formula = om ~ .,importance = "permutation")
vip(mod)


# Logistic regressions
mod_log <- 
  glm(data=am,
      formula = am ~ lat_dd + long_dd + am_em_dom + year + total_precip_mm + sample_type,
      family=binomial)
am$preds <- predict(mod_log,am, type='response')
broom::tidy(mod_log) %>% 
  write.csv("./output/probability_of_AMF_presence_model_table.csv")
mod_log_ecm <- 
  glm(data=ecm,
      formula = ecm ~ lat_dd + long_dd + am_em_dom + year + total_precip_mm ,
      family=binomial)
broom::tidy(mod_log_ecm) %>% 
  write.csv("./output/probability_of_ECM_presence_model_table.csv")
mod_log_erm <- 
  glm(data=erm,
      formula = erm ~ lat_dd + long_dd + am_em_dom + year + total_precip_mm,
      family=binomial)
broom::tidy(mod_log_erm) %>% 
  write.csv("./output/probability_of_ERM_presence_model_table.csv")
mod_log_om <- 
  glm(data=om,
      formula = om ~ lat_dd + long_dd + am_em_dom + year + total_precip_mm,
      family=binomial)
broom::tidy(mod_log_om) %>% 
  write.csv("./output/probability_of_OM_presence_model_table.csv")
