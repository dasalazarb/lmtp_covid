## ---------------------------
##
## Script name: Wrangle oxygen time
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
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)


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

## ---------------------------

## load data

file_date <- "2020-08-17"
file_path <- paste0("data/raw/",file_date)

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/hospitalized_cohort_with_queens.rds"))

# vital signs, with only the relevant respiratory + general vs kept
vitals <- read_rds(here::here(file_path, "wcm_vitals_filtered_rel.rds"))

ids_with_vitals <- vitals %>%
  filter(empi %in% cohort$empi) %>%
  distinct(empi) %>%
  pull(empi)

cohort <- cohort %>% filter(empi %in% ids_with_vitals)

# DNR/DNI status
dnr <-
  read_rds(here::here(paste0(file_path, "/dnr.rds"))) %>%
  mutate(empi = as.character(empi))

# set time zones
tz(dnr$procedure_dt_tm) <- "America/New_York"
tz(vitals$recordeddtm) <- "America/New_York"

## ---------------------------

## functions

# DAYS TEMPLATE - GENERAL - 24 HOURS RELATIVE TO ED

# make periods (days) in redcap data
create_periods = function(ed_adm_dt, end_dt, ...){
  seq.POSIXt(from = floor_date(ed_adm_dt, unit="days"), to = ceiling_date(end_dt, unit="days"), by = "days") %>% as_tibble()
}

## ---------------------------

# make a day (with start and end date-times for intervals) to merge with steroid data
df_days <-
  cohort %>%
  select(empi, ed_adm_dt, end_dt) %>%
  mutate(ed_adm_dt = as.Date(floor_date(ed_adm_dt, unit="days")), end_date = as.Date(floor_date(end_dt, unit="days")),
         day_start_dt = ed_adm_dt) %>%
  group_by(empi) %>%
  complete(day_start_dt = seq.Date(ed_adm_dt, end_date, by="day")) %>%
  # starting from hospitalizationday, make 24 hour periods and call the first "day 0"
  mutate(day = row_number() - 1,
         day_end_dt = day_start_dt + hours(24),
         max_day = max(day)) %>%
  select(empi, day_start_dt, day_end_dt, day, max_day) %>% 
  ungroup() %>%
  filter(day <= 60)

# add DNR/DNI in
dnr_day <-
  inner_join(dnr, df_days) %>%
  filter(procedure_dt_tm %within% interval(day_start_dt, day_end_dt)) %>%
  mutate(day_dnrdni = day,
         dnrdni = "yes",
         dnrdni_date = as.Date(procedure_dt_tm))  %>%
  group_by(empi) %>%
  arrange(empi, procedure_dt_tm) %>%
  filter(row_number() == 1) %>%
  select(empi, day_dnrdni, dnrdni, procedure_description, dnrdni_date)

o2_clean <- 
  cohort %>%
  select(empi, ed_adm_dt, end_dt) %>%
  left_join(vitals) %>%
  filter(observation_name %in%
           c("resp_non vent device",
             "resp_insp_gas_label_oxygen",
             "vs_resp_devicepacu"),
         recordeddtm %within% interval(ed_adm_dt, end_dt)) %>%
  mutate(any_o2 = case_when(
           observation_name == "resp_insp_gas_label_oxygen" &
             parse_number(valuetext) > 0 ~ "Yes",
           observation_name == "resp_non vent device" &
             ((valuetext != "0") |  !str_detect(tolower(valuetext), "room"))  ~ "Yes",
           observation_name == "vs_resp_devicepacu" &
             !str_detect(tolower(valuetext), "room")  ~ "Yes",
         TRUE ~ "No"
         )) %>%
  group_by(empi) %>%
  select(empi, recordeddtm, any_o2) %>%
  filter(any_o2 == "Yes")

o2_by_day <- df_days %>%
  inner_join(o2_clean) %>% # only keep pts with at least one o2 supp 
  filter(recordeddtm %within% interval(day_start_dt, day_end_dt)) %>% # only keep days with an o2 supp
  distinct(empi, day, .keep_all = T) %>%
  right_join(df_days) %>% # add the other days back in
  mutate(any_o2 = case_when(is.na(any_o2) ~ "No", TRUE ~ any_o2))


events_day <-
  cohort %>%
  select(empi, death, death_date, ed_adm_dt,
         intubation1_date, extubation1_date, intubation2_date, extubation2_date, end_dt) %>%
  mutate(across(contains("date"), ~as.Date(.x))) %>%
  full_join(o2_by_day) %>%
  full_join(dnr_day) %>%
  group_by(empi) %>%
  mutate(dnrdni = replace_na(dnrdni, "No"),
         death_this_day = case_when(death_date == day_start_dt ~ day),
         o2_this_day = case_when(any_o2 == "Yes" ~ day),
         intubation1_this_day = case_when(intubation1_date == day_start_dt ~ day),
         extubation1_this_day = case_when(extubation1_date %within% interval(day_start_dt, day_end_dt) ~ day),
         intubation2_this_day = case_when(intubation2_date %within% interval(day_start_dt, day_end_dt) ~ day),
         extubation2_this_day = case_when(extubation2_date %within% interval(day_start_dt, day_end_dt) ~ day)) %>%
  fill(intubation1_this_day, intubation2_this_day, extubation1_this_day, extubation2_this_day, end_dt,
       .direction = "downup") %>%
  mutate(
    max_day_deceased = case_when(death == "Yes" ~ max_day),
    max_day_discharged = case_when(death == "No" ~ max_day)
  ) %>%
  ungroup() %>%
  mutate(
    empi = fct_reorder(empi, max_day),
    intubation_status = case_when(day >= intubation1_this_day & day <= extubation1_this_day ~ "Intubated", 
                                  day >= intubation1_this_day & is.na(extubation1_this_day) ~ "Intubated", 
                                  day >= intubation2_this_day & day <= extubation2_this_day ~ "Intubated", 
                                  day >= intubation2_this_day & is.na(extubation2_this_day) ~ "Intubated",
                                  TRUE ~ "Not intubated")) %>%
  # there were a ~600 mismatches where pts were intubated and not recorded in VITALS as intubated
  # overrode this 
  mutate(any_o2 = case_when(intubation_status == "Intubated" ~ "Yes", TRUE ~ any_o2),
         o2_status = case_when(intubation_status == "Intubated" ~ "Intubated",
                               any_o2 == "Yes" ~ "Supp O2",
                               any_o2 == "No" ~ "None"))



## ---------------------------
# Compute weights for steroids

# keep only patients weights
weights <- vitals %>%
  filter(observation_name == "vs_weight_amb_kg_cal") %>%
  right_join(cohort %>% select(empi, sex, ed_adm_dt, ed_adm_dt)) %>%
  # regular and absolute value of the difference in time. positive is after time of meeting o2 crit
  mutate(
    diff_ed_adm_dt_weight_time = difftime(recordeddtm, ed_adm_dt, units="day"),
    abs_diff_ed_adm_dt_weight_time = abs( diff_ed_adm_dt_weight_time)) %>%
  group_by(empi) %>%
  arrange(empi, abs_diff_ed_adm_dt_weight_time) %>%
  filter(row_number() == 1) %>%
  mutate(weight_kg = parse_number(valuetext),
         weight_kg = case_when(weight_kg == 0 ~ NA_real_,
                               TRUE ~ weight_kg)) %>%
  select(empi, sex, weight_kg, weight_time = recordeddtm, diff_ed_adm_dt_weight_time) %>%
  group_by(sex) %>%
  # Impute weight for males and females using averages of other males and females
  mutate(
    imputed_weight = case_when(is.na(weight_kg) ~ "Yes", TRUE ~ "No"),
    weight_kg = case_when(is.na(weight_kg) ~ mean(weight_kg, na.rm=T),
                          TRUE ~ weight_kg)
  ) %>%
  ungroup()


# write rds of oxygen time
write_rds(cohort, "data/derived/cohort_w_vitals.rds")
write_rds(o2_by_day, "data/derived/o2_time.rds")
write_rds(dnr_day, "data/derived/dnrdni_time.rds")
write_rds(df_days, "data/derived/df_days.rds")
write_rds(events_day, "data/derived/events_day.rds")
write_rds(weights, "data/derived/weights.rds")

