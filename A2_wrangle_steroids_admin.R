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

# medications - steroid administration
mar <-
  read_rds(here::here(paste0(file_path, "/covid_datalake_wcm_mar.rds"))) %>%
  mutate(empi = as.character(empi))

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_vitals.rds"))
cohort <- cohort %>%
  filter(empi %in% mar$empi) # remove 5 more pts

events_day <- read_rds(here::here("data/derived/events_day.rds"))

weights <- read_rds(here::here("data/derived/weights.rds"))

df_days <- read_rds(here::here("data/derived/df_days.rds"))

o2_time <- read_rds(here::here("data/derived/o2_time.rds"))

# set time zones
#tz(labs$result_time) <- "EST"
tz(mar$start_dt_tm) <- "America/New_York"
tz(mar$end_dt_tm) <- "America/New_York"

# add weights into full data set
cohort <- cohort %>%
  left_join(weights)

# Begin steroid cleaning

# Steroids

## Steroid doses of interest (by time stamp)

# Doses were first converted to methlprednisolone, then if <10 mg, converted to mg/kg body weight. For patients missing body weight, the average weight (stratified by male or female) was used.
# 
# Methylprednisolone conversions:
#   
#   - prednisoLONE and predniSONE: `dose * (4/5)`
# 
# - hydrocortisone: `dose * (4/20)`
# 
# - dexamethasone: `dose * (4/0.75)`


# steroids of interest per email convo with Ed
steroid_meds <- c(
  "Dexamethasone Inj",
  "Dexamethasone Inj 10 mg/ml for oral use",
  "Dexamethasone Oral Liq",
  "Dexamethasone Oral Liq Peds ED",
  "Dexamethasone Tab",
  "Dexamethasone Tab Peds ED",
  #"Fludrocortisone Oral",
  "Hydrocortisone Oral",
  "Hydrocortisone Oral Susp",
  "Hydrocortisone Succinate Inj",
  "Methylprednisolone Oral",
  "Methylprednisolone Sod Suc Inj",
  "Methylprednisolone Sod Suc Inj Peds ED",
  "prednisoLONE Oral Liq",
  "prednisoLONE Oral Liq Peds ED",
  "predniSONE Oral",
  "predniSONE Oral Liq",
  "predniSONE oral Peds ED"
)

steroids_cohort_raw <-
  cohort %>%
  select(empi, weight_kg, ed_adm_dt, end_dt) %>%
  left_join(mar) %>%
  filter(medication_name %in% steroid_meds)

# for examining doses with Ed
steroid_dose_ed <-
  weights %>%
  select(empi, weight_kg) %>%
  left_join(mar) %>%
  # keep only meds in the steroids list
  filter(medication_name %in% steroid_meds) %>%
  count(generic_name, dose, route) %>%
  arrange(desc(n)) %>% 
  drop_na(dose)

# use the MAR from the data lake
steroids <-
  steroids_cohort_raw %>%
  mutate(dose=as.numeric(dose)) %>%
  mutate(mp_dose = case_when(generic_name == "methylPREDNISolone" ~ dose,
                             generic_name %in% c("prednisoLONE", "predniSONE") ~ dose * (4/5),
                             generic_name == "hydrocortisone" ~ dose * (4/20),
                             generic_name == "dexamethasone" ~ dose * (4/0.75)
                             # FIGURE OUT WHAT TO DO WITH FLUDROCORTISONE?
  ),
  # may group p and mp
  # critiques of hernan's work - confounding by indication impossible to capture
  mp_dose_kg = mp_dose / weight_kg
  ) %>%
  select(empi, medication_name, generic_name, route, dose, mp_dose, weight_kg, mp_dose_kg, start_dt_tm, end_dt_tm, ed_adm_dt, end_dt) %>%
  arrange(empi, start_dt_tm) %>%
  # keep only start times of medications that fall within their ED admission date and last followup date
  filter(start_dt_tm %within% interval(ed_adm_dt, end_dt)) 

# GEOM_DENSITY_RIDGES
# steroids %>%
#   filter(!is.infinite(mp_dose_kg)) %>%
#   mutate(id = paste(tolower(generic_name), tolower(route))) %>%
#   ggplot(aes(mp_dose_kg, id)) +
#   geom_density_ridges() +
#   labs(x = "Methylprednisolone equivalent dose / kg",y="",
#        title = "Distributions of methylprednsolone equivalents")

steroids %>%
  filter(!is.infinite(mp_dose_kg)) %>%
  mutate(id = paste(tolower(generic_name), tolower(route))) %>%
  ggplot(aes(mp_dose_kg)) +
  geom_histogram() +
  facet_wrap(~id, scales = "free", ncol = 3) +
  labs(x = "Methylprednisolone equivalent dose / kg",y="Number of instances in patient charts",
       title = "Instances of methylprednsolone equivalents")


## Determining steroids ≥ 0.5 mg per kg per day

# - Made day intervals based upon every 24 hour window since ED admission
# 
# - If at any point during a day window, a rolling 24 hour steroid dose became ≥ 0.5 mg/kg, that patient was marked as a "yes" for treatment for that day


#### STEROIDS ####

# get each patients cumulative dose and rolling 24 hour dose
steroids_cumsum <-
  steroids %>%
  group_by(empi) %>%
  arrange(start_dt_tm) %>%
  # total cumulative sum
  mutate(cumsum = cumsum(mp_dose_kg)) %>%
  # make a row for every hour (60 min*60 sec) from the earliest to latest steroid administration time
  complete(start_dt_tm = seq(min(start_dt_tm), max(start_dt_tm), 60*60), fill = list(mp_dose_kg=0)) %>%
  # do a rolling sum for every 24 rows (1 hour)
  mutate(cum_rolling24 = rollapplyr(mp_dose_kg, width=24, FUN = sum, partial = T)) %>%
  # remove the hours that had no doses administered
  drop_na(cumsum)



steroids_alldays <-
  # keep all the days, add the steroids back in
  df_days %>%
  left_join(steroids_cumsum) %>%
  # only keep values in the interval of the days between hospitalizationand last follow up
  filter(start_dt_tm %within% interval(day_start_dt, day_end_dt)) %>%
  group_by(empi, day) %>%
  # keep only the max rolling dose in that day
  filter(cum_rolling24 == max(cum_rolling24)) %>%
  # but some have equivalent (due to doses being staggered so that they maintain a rolling sum of X, so only keep first row after that)
  filter(row_number() == 1)

# people who received at least one dose of steroids ≥ .5 mg/kg body weight
steroid_count <-
  steroids_alldays %>%
  filter(cum_rolling24 >= .5) %>%
  group_by(empi) %>%
  filter(row_number() == 1)

steroids_over1mg_days <-
  steroids_alldays %>%
  mutate(steroid_1mg = ifelse(cum_rolling24 >= .5 , 1, 0)) %>%
  filter(cum_rolling24 >= .5) %>%
  group_by(empi) %>%
  arrange(empi, day) %>%
  mutate(n_steroids_1mg = cumsum(steroid_1mg),
         total_steroids_1mg = sum(steroid_1mg),
         any_steroids_1mg = ifelse(total_steroids_1mg > 0, 1, 0)) %>%
  filter(any_steroids_1mg == 1) %>%
  arrange(total_steroids_1mg) %>%
  ungroup() %>%
  mutate(day_steroids = day) %>%
  left_join(cohort %>% select(empi, mrn)) %>%
  select(mrn, empi, day, everything()) %>%
  group_by(empi) %>%
  fill(any_steroids_1mg, .direction="downup") %>%
  mutate(any_steroids_1mg = replace_na(any_steroids_1mg, 0),
         generic_name = replace_na(generic_name, "None"),
         # generic_name = fct_relabel(generic_name,
         #                            Methylprednisolone = Methylprednisolone,
         #                            Prednisone = predniSONE,
         #                            Dexamethasone = dexamethasone,
         #                            Hydrocortisone = hydrocortisone),
         generic_name = fct_relevel(generic_name, "None")
  ) %>%
  mutate(any_steroids_1mg = case_when(any_steroids_1mg == 1 ~ "Steroids ever", TRUE ~ "Steroids never")) %>%
  ungroup() %>%
  arrange(empi, day_start_dt)


steroids_over1mg_days %>%  
  mutate(drug_route = paste(tolower(generic_name), tolower(route))) %>%
  filter(cum_rolling24 > 5) %>%
  count(drug_route) %>%
  gt() %>%
  tab_header("Drug and route of instances where max cumulative rolling 24 hour value > 5 mg/kg per day")


steroids_over1mg_days %>%  
  mutate(drug_route = paste(tolower(generic_name), tolower(route))) %>%
  filter(cum_rolling24 > 15) %>%
  count(drug_route) %>%
  gt() %>%
  tab_header("Drug and route of patients where max cumulative rolling 24 hour value > 15 mg/kg per day")
# people who received at least one dose of steroids ≥ .5 mg/kg body weight
steroid_count <-
  steroids_alldays %>%
  filter(cum_rolling24 >= .5) %>%
  group_by(empi) %>%
  filter(row_number() == 1)

steroids_over1mg_days <-
  steroids_alldays %>%
  mutate(steroid_1mg = ifelse(cum_rolling24 >= .5 , 1, 0)) %>%
  filter(cum_rolling24 >= .5) %>%
  group_by(empi) %>%
  arrange(empi, day) %>%
  mutate(n_steroids_1mg = cumsum(steroid_1mg),
         total_steroids_1mg = sum(steroid_1mg),
         any_steroids_1mg = ifelse(total_steroids_1mg > 0, 1, 0)) %>%
  filter(any_steroids_1mg == 1) %>%
  arrange(total_steroids_1mg) %>%
  ungroup() %>%
  mutate(day_steroids = day) %>%
  left_join(cohort %>% select(empi, mrn)) %>%
  select(mrn, empi, day, everything()) %>%
  group_by(empi) %>%
  fill(any_steroids_1mg, .direction="downup") %>%
  mutate(any_steroids_1mg = replace_na(any_steroids_1mg, 0),
         generic_name = replace_na(generic_name, "None"),
         # generic_name = fct_relabel(generic_name,
         #                            Methylprednisolone = Methylprednisolone,
         #                            Prednisone = predniSONE,
         #                            Dexamethasone = dexamethasone,
         #                            Hydrocortisone = hydrocortisone),
         generic_name = fct_relevel(generic_name, "None")
  ) %>%
  mutate(any_steroids_1mg = case_when(any_steroids_1mg == 1 ~ "Steroids ever", TRUE ~ "Steroids never")) %>%
  ungroup()



events_day_w_steroids <-
  events_day %>%
  full_join(steroids_over1mg_days) %>%
  mutate(steroid_this_day = case_when(!is.na(steroid_1mg) ~ day),
         empi = fct_reorder(empi, max_day)) %>%
  group_by(empi) %>%
  fill(any_steroids_1mg, .direction = "downup") %>%
  mutate(any_steroids_1mg = replace_na(any_steroids_1mg, "Steroids never"),
         # o2_crit_this_day = case_when(day_start_dt == as.Date(o2_crit_time) ~ day),
         steroids_within_5d = case_when(day <= 5 & !is.na(steroid_this_day) ~ "Steroids within 5d"),
         day_steroids_after_5d = case_when(day > 5 & !is.na(steroid_this_day) & steroids_within_5d == "No steroids within 5d" ~ day
                                       )) %>%
  fill(steroids_within_5d,
       #o2_crit_this_day,
       .direction = "downup") %>%
  ungroup()  %>%
  mutate(death_yn = case_when(death == "Yes" ~ "Died in hospital",
                              death == "No" ~ "Discharged alive"),
         generic_name = replace_na(as.character(generic_name),""),
         generic_name = fct_relevel(factor(generic_name), "", after=Inf),
         steroids_within_5d = replace_na(steroids_within_5d, "No steroids within 5d")#,
         
         #day_from_o2 = day - o2_crit_this_day
         ) # %>%
  #group_by(empi) %>%
  #mutate(steroids_within_5d_hypoxia = case_when(day_from_o2 <= 5 & !is.na(steroid_this_day) ~ "Steroids within 5d hypoxia"),
  #       day_steroids_after_5d) %>%
  #fill(steroids_within_5d_hypoxia, .direction = "downup") %>%
  #mutate(steroids_within_5d_hypoxia = replace_na(steroids_within_5d_hypoxia, "No steroids within 5d hypoxia")) %>%
  #ungroup()

# events_day %>%
#   drop_na(death_yn) %>%
#   ggplot(aes(day, empi, group=empi, col=intubation_status)) +
#   geom_line(alpha=.5) +
#   geom_point(aes(death_this_day,empi), shape=4, size=.4, col="black") +
#   geom_point(aes(day_dnrdni,empi), shape=21, col="purple",
#              #fill="white",
#              size=.5) +
#   geom_point(aes(steroid_this_day, empi, shape=generic_name), col="red", 
#              #fill="white",
#              size=.5) +
#   scale_x_continuous(limits=c(0,30), minor_breaks=0:30, breaks=seq(0,30,2), expand=c(.01,.01)) +
#   labs(x="Days since hypoxia", y= "Individual Patient", col = "Intubation status", shape="Generic name",
#        title = "Steroids, Death, DNR, and Intubation statuses relative to hospitalization day",
#        subtitle = "Patients split by final mortality status") +
#   theme(#legend.position="top",
#     axis.text=element_text(size=6),
#     axis.text.y = element_blank()) +
#   scale_color_manual(values = c("dodgerblue","navy")) +
#   facet_wrap(~death_yn, scales = "free") +
#   scale_shape_manual(values=c(22,23,24,25,20,NA))

events_day_w_steroids %>%
  drop_na(death_yn) %>%
  ggplot(aes(day, empi, group=empi, col=intubation_status)) +
  geom_line(alpha=.5) +
  geom_point(aes(death_this_day,empi), shape=4, size=.4, col="black") +
  geom_point(aes(day_dnrdni,empi), shape=21, col="purple",
             #fill="white",
             size=.5) +
  geom_point(aes(steroid_this_day, empi, shape=generic_name), col="red", 
             #fill="white",
             size=.5) +
  scale_x_continuous(limits=c(0,30), minor_breaks=0:30, breaks=seq(0,30,2), expand=c(.01,.01)) +
  labs(x="Days since hypoxia", y= "Individual Patient", col = "Intubation\nStatus", shape = "Steroid\nAdministration",
       title = "Steroids, Death, DNR, Intubation statuses relative to hospitalization day", subtitle="Patients split by ever/never steroids") +
  theme(#legend.position="top",
    axis.text=element_text(size=6),
    axis.text.y = element_blank()) +
  scale_color_manual(values = c("dodgerblue","navy")) +
  facet_wrap(~any_steroids_1mg, scales = "free") +
  scale_shape_manual(values=c(22,23,24,25,20,NA))


steroids_patient_list <- events_day_w_steroids %>%
  group_by(empi) %>%
  summarize(steroid_list = paste(generic_name, collapse = ",")) %>%
  mutate(methylpred_ever = case_when(str_detect(steroid_list, "methyl") ~ "Yes",
                                     TRUE ~ "No"),
         hydro_ever = case_when(str_detect(steroid_list, "hydro") ~ "Yes",
                                TRUE ~ "No"),
         pred_ever = case_when(str_detect(steroid_list, "pred") ~ "Yes",
                               TRUE ~ "No"),
         dex_ever = case_when(str_detect(steroid_list, "dex") ~ "Yes",
                              TRUE ~ "No")
  ) %>%
  select(-steroid_list)

write_rds(steroids_over1mg_days, "data/derived/steroids_over_1mg_days.rds")
write_rds(events_day_w_steroids, "data/derived/events_day_w_steroids.rds")
write_rds(cohort, "data/derived/cohort_w_mar.rds")
