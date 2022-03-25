## ---------------------------
##
## Script name: Wrangle steroids administration
##
## Purpose of script:
##
## Author: Katherine Hoffman
##
## Date Created: 2021-04-05
##
## Copyright (c) Katherine Hoffman, 2021
## Email: kah2797@med.cornell.edu
##
##
## Notes:
##   
##

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation


# load packages -----------------------------------------------------------

library(lubridate)
library(wesanderson)
library(tictoc)
library(dtplyr)
library(tidyverse)
library(tidylog)
library(zoo)
library(gtsummary)
library(gt)
library(ggridges)
library(labelled)
library(snakecase)
library(janitor)


# load data -----------------------------------------------------------

file_date <- "2020-08-17"
file_path <- paste0("data/raw/",file_date)

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds"))


# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
df_days <- read_rds(here::here("data/derived/df_days.rds")) %>%
  filter(empi %in% cohort$empi)

# vital signs, with only the relevant respiratory + general vs kept
vitals <- read_rds(here::here(file_path, "wcm_vitals_filtered_rel.rds"))

vital_params <- c(
  "vs_hr_hr",
  "xp_resp_spo2",
  "xp_resp_rate_pt",
  "vs_bp_noninvasive (s)",
  "vs_bp_noninvasive (d)",
  "VS_Temperature (C) calc"
)


vitals_low_long <-
  vitals %>%
  filter(observation_name %in% vital_params) %>%
  select(empi, recordeddtm, observation_name, valuetext) %>%
  right_join(df_days) %>%
  filter(recordeddtm %within% interval(day_start_dt, day_end_dt)) %>%
  mutate(valuetext = parse_number(valuetext),
         observation_name = case_when(observation_name == "vs_hr_hr" ~ "heart_rate",
                                 observation_name == "xp_resp_spo2" ~ "spo2",
                                 observation_name == "xp_resp_rate_pt" ~ "resp_rate",
                                 observation_name == "vs_bp_noninvasive (s)" ~ "blood_pressure_s",
                                 observation_name == "vs_bp_noninvasive (d)" ~ "blood_pressure_d",
                                 observation_name == "vs_bp_noninvasive (m)" ~ "blood_pressure_m",
                                 observation_name == "VS_Temperature (C) calc" ~ "temperature"))  %>%
  group_by(empi, observation_name, day) %>%
  arrange(desc(valuetext)) %>%
  filter(row_number() == 1)

vitals_high_long <-
  vitals %>%
  filter(observation_name %in% vital_params) %>%
  select(empi, recordeddtm, observation_name, valuetext) %>%
  right_join(df_days) %>%
  filter(recordeddtm %within% interval(day_start_dt, day_end_dt)) %>%
  mutate(valuetext = parse_number(valuetext),
         observation_name = case_when(observation_name == "vs_hr_hr" ~ "heart_rate",
                                      observation_name == "xp_resp_spo2" ~ "spo2",
                                      observation_name == "xp_resp_rate_pt" ~ "resp_rate",
                                      observation_name == "vs_bp_noninvasive (s)" ~ "blood_pressure_s",
                                      observation_name == "vs_bp_noninvasive (d)" ~ "blood_pressure_d",
                                      observation_name == "vs_bp_noninvasive (m)" ~ "blood_pressure_m",
                                      observation_name == "VS_Temperature (C) calc" ~ "temperature"))  %>%
  group_by(empi, observation_name, day) %>%
  arrange(desc(valuetext)) %>%
  filter(row_number() == n())

vitals_high_wide <-
  vitals_high_long %>%
  pivot_wider(id_cols = c("empi","day","max_day","day_start_dt","day_end_dt"),
              names_from = observation_name,
              values_from = valuetext,
              names_prefix = "daily_highest_"
              ) %>%
  clean_names() %>%
  arrange(empi, day)

vitals_low_wide <-
  vitals_low_long %>%
  pivot_wider(id_cols = c("empi","day","max_day","day_start_dt","day_end_dt"),
              names_from = observation_name,
              values_from = valuetext,
              names_prefix = "daily_lowest_"
  ) %>%
  clean_names() %>%
  arrange(empi, day)

vitals_high_long_clean <- vitals_high_wide %>% 
  pivot_longer(cols = starts_with("daily_")) %>% 
  drop_na(value)

vitals_low_long_clean <- vitals_low_wide %>% 
  pivot_longer(cols = starts_with("daily_")) %>% 
  drop_na(value)

vitals_long <- full_join(vitals_high_long_clean, vitals_low_long_clean) %>%
  ungroup() %>%
  mutate(vs = str_remove(name, "daily_"))

vitals_long %>%
  left_join(cohort %>% select(empi, red_cap_source)) %>%
  count(day, vs, red_cap_source) %>%
  ggplot(aes(day, n, fill=red_cap_source)) +
  geom_bar(stat="identity") +
  facet_wrap(~vs) +
  theme_bw()

vitals_wide <- full_join(vitals_high_wide, vitals_low_wide)

write_rds(vitals_long, "data/derived/vitals_long.rds")
write_rds(vitals_wide, "data/derived/vitals_wide.rds")
