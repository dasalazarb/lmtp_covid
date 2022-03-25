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
cohort <- read_rds(here::here("data/derived/cohort_w_mar.rds"))

df_days <- read_rds(here::here("data/derived/df_days.rds"))

labs <-
  read_rds(here::here(paste0(file_path, "/wcm_labs_filtered.rds"))) %>%
  mutate(empi = as.character(empi))

cohort <- cohort %>%
  filter(empi %in% labs$empi)

low_lab_params <- c(
  "Platelet",
  "Absolute Neutrophil Count",
  "Neutrophil Automated",
  "Glucose Whole Blood  Meter POC",
  "Glucose Whole Blood Meter POC",
  "pO2 (Arterial) - EPOC",
  "pO2 Arterial",
  "pCO2 (Arterial) - EPOC",
  "pCO2 Arterial",
  "P CO2 Arterial",
  "P O2 Arterial"
)

high_lab_params <- c(
  "Platelet",
  "Lymphocyte Automated",
  "Absolute Lymphocyte Count",
  "Fibrinogen",
  "D-Dimer",
  "C-Reactive Protein",
  "C-Reactive Protein High Sensitivity",
  "Ferritin",
  "Absolute Neutrophil Count",
  "Neutrophil Automated",
  "Prothrombin Time",
  "Activated Partial Thromboplastin Time",
  "Troponin-I",
  "Glucose Whole Blood  Meter POC",
  "Glucose Whole Blood Meter POC",
  "Bilirubin Total",
  "BUN/Creatinine Ratio",
  "Creatine Kinase",
  "Creatinine"
)


# labs %>%
#     filter(str_detect(tolower(result_name), "reactive")) %>%
#   count(result_name)

labs_low_long <-
  labs %>%
  filter(result_name %in% low_lab_params) %>%
  select(empi, result_time, result_name, ord_value) %>%
  right_join(df_days) %>%
  filter(result_time %within% interval(day_start_dt, day_end_dt)) %>%
  mutate(ord_value = parse_number(ord_value),
         result_name = case_when(str_detect(result_name, "Glucose") ~ "Glucose",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Lactate") ~ "Lactate",
                                 str_detect(result_name, "pH") ~ "Arterial Ph",
                                 str_detect(result_name, "CO2") ~ "Arterial PCO2",
                                 str_detect(result_name, "O2") ~ "Arterial PaO2",
                                 str_detect(result_name, "Bilirubin") ~ "Bilirubin",
                                 TRUE ~ result_name))  %>%
  group_by(empi, result_name, day) %>%
  arrange(desc(ord_value)) %>%
  filter(row_number() == 1)

labs_high_long <-
  labs %>%
  filter(result_name %in% high_lab_params) %>%
  select(empi, result_time, result_name, ord_value) %>%
  right_join(df_days) %>%
  filter(result_time %within% interval(day_start_dt, day_end_dt)) %>%
  mutate(ord_value = parse_number(ord_value),
         result_name = case_when(str_detect(result_name, "Glucose") ~ "Glucose",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Reactive") ~ "C-Reactive Protein",
                                 str_detect(result_name, "Lactate") ~ "Lactate",
                                 str_detect(result_name, "pH") ~ "Arterial Ph",
                                 str_detect(result_name, "CO2") ~ "Arterial PCO2",
                                 str_detect(result_name, "O2") ~ "Arterial PO2",
                                 str_detect(result_name, "Bilirubin") ~ "Bilirubin",
                                 TRUE ~ result_name))  %>%
  group_by(empi, result_name, day) %>%
  arrange(desc(ord_value)) %>%
  filter(row_number() == n())

labs_high_wide <-
  labs_high_long %>%
  pivot_wider(id_cols = c("empi","day","max_day","day_start_dt","day_end_dt"),
              names_from = result_name,
              values_from = ord_value,
              names_prefix = "daily_highest_"
              ) %>%
  clean_names() %>%
  arrange(empi, day)

labs_low_wide <-
  labs_low_long %>%
  pivot_wider(id_cols = c("empi","day","max_day","day_start_dt","day_end_dt"),
              names_from = result_name,
              values_from = ord_value,
              names_prefix = "daily_lowest_"
  ) %>%
  clean_names() %>%
  arrange(empi, day)

labs_wide <- full_join(labs_high_wide, labs_low_wide)

labs_high_long_clean <- labs_high_wide %>% 
  pivot_longer(cols = starts_with("daily_")) %>% 
  drop_na(value)

labs_low_long_clean <- labs_low_wide %>% 
  pivot_longer(cols = starts_with("daily_")) %>% 
  drop_na(value)

labs_long <- full_join(labs_high_long_clean, labs_low_long_clean) %>%
  ungroup() %>%
  mutate(lab = str_remove(name, "daily_"))

labs_long %>%
  count(day, lab) %>%
  ggplot(aes(day, n)) +
  geom_bar(stat="identity") +
  facet_wrap(~lab) +
  theme_bw()

write_rds(labs_long, "data/derived/labs_long.rds")
write_rds(labs_wide, "data/derived/labs_wide.rds")
write_rds(cohort, "data/derived/cohort_w_labs.rds")
