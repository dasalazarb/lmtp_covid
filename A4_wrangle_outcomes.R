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

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds"))

o2_time <- read_rds(here::here("data/derived/o2_time.rds"))

steroids_over1mg_days <- read_rds(here::here("data/derived/steroids_over_1mg_days.rds"))

#events_day_w_steroids <- read_rds(here::here("data/derived/events_day_w_steroids.rds"))

outcomes <- 
  cohort %>% 
  left_join(o2_time) %>%
  select(empi, ed_adm_dt, end_dt, death) %>%
  #mutate(#days_from_hypoxia = ceiling(time_length(difftime(end_dt, o2_crit_time), unit="day")),
         #event_death_28d_from_hypoxia = case_when(days_from_hypoxia <= 28 & death == "Yes" ~ 1,
         #                                         TRUE ~ 0),
         # if the patient doesn't die in 28 days, extend the days to 29 (assume no death)
         #days_from_hypoxia_28d_ext = case_when(days_from_hypoxia <= 28 & event_death_28d_from_hypoxia == 1 ~ days_from_hypoxia,
         #                                      TRUE ~ 28),
         # follow up time is censored at discharge
         #days_from_hypoxia_28d = case_when(days_from_hypoxia <= 28 ~ days_from_hypoxia,
         #                                  TRUE ~ 28)) %>%
  # add days from hospitalization to death
  mutate(days_from_hosp = ceiling(time_length(difftime(end_dt, ed_adm_dt), unit="day")),
         event_death_28d_from_hosp = case_when(days_from_hosp <= 28 & death == "Yes" ~ 1,
                                               TRUE ~ 0),
         # if the patient doesn't die in 28 days, extend the days to 29 (assume no death)
         days_from_hosp_28d_ext = case_when(days_from_hosp <= 28 & event_death_28d_from_hosp == 1 ~ days_from_hosp,
                                            TRUE ~ 28),
         # follow up time is censored at discharge
         days_from_hosp_28d = case_when(days_from_hosp <= 28 ~ days_from_hosp,
                                        TRUE ~ 28))

# for the flat file - don't use events_day because you can't rename day_start_dt for people who ddn't have steroids
first_steroids <-
  steroids_over1mg_days %>%
  select(empi, first_steroid_date = day_start_dt) %>%
  arrange(empi, first_steroid_date) %>%
  distinct(empi, .keep_all = T) 

  
outcome_exposure <-
  outcomes %>%
  left_join(first_steroids) %>%
  mutate(time_to_steroids = ceiling(time_length(difftime(first_steroid_date, ed_adm_dt), unit="day"))) %>%
         #time_to_steroids_hypoxia = ceiling(time_length(difftime(first_steroid_date, o2_crit_time), unit="day"))) %>%
# add point treatment exposure
   mutate(steroids_within_5d = case_when(round(time_length(difftime(first_steroid_date, ed_adm_dt), unit="day")) <= 5 ~ 1,
                                      TRUE ~ 0)) %>%
        # steroids_within_5d_hypoxia = case_when(is.na(o2_crit_time) ~ NA_real_,
        #                                round(time_length(difftime(first_steroid_date, o2_crit_time), unit="day")) <= 5 ~ 1,
        #                                TRUE ~ 0))  %>%
  # censor at the date of steroids when controls get steroids
  mutate(days_from_hosp_28d_cens_ster = case_when(steroids_within_5d == 0 & !is.na(first_steroid_date) ~
                                                   ceiling(time_length(difftime(first_steroid_date, ed_adm_dt), unit="day")),
                                                  TRUE ~ days_from_hosp_28d),
         # those patients wouldn't have an event bc they're censored
         event_death_28d_from_hosp_cens_ster = case_when(steroids_within_5d == 0 & !is.na(first_steroid_date) ~ 0,
                                                         TRUE ~ event_death_28d_from_hosp),
         # change the starting time for the treatment group to the day of steroids
         days_from_hosp_ster_28d_cens_ster = case_when(steroids_within_5d == 0 & 
                                                         ceiling(time_length(difftime(end_dt, first_steroid_date), unit="day")) <= 28 ~
                                                         ceiling(time_length(difftime(end_dt, first_steroid_date), unit="day")),
                                                       steroids_within_5d == 0 ~ 28,
                                                       TRUE ~ days_from_hosp_28d_cens_ster),
            # steroids_within_5d == 0 & !is.na(first_steroid_date) ~
            #                                         ceiling(time_length(difftime(first_steroid_date, ed_adm_dt), unit="day")),
            #                                       TRUE ~ days_from_hosp_28d),s
         days_from_hosp_28d_cens_ster = case_when(days_from_hosp_28d_cens_ster  > 28 ~ 28, TRUE ~ days_from_hosp_28d_cens_ster )) # ,
         
         #days_from_hypoxia_28d_cens_ster = case_when(steroids_within_5d_hypoxia == 0 & !is.na(first_steroid_date) ~
         #                                             ceiling(time_length(difftime(first_steroid_date, o2_crit_time), unit="day")),
         #                                         TRUE ~ days_from_hypoxia_28d),
         #days_from_hypoxia_28d_cens_ster = case_when(days_from_hypoxia_28d_cens_ster  > 28 ~ 28, TRUE ~ days_from_hypoxia_28d_cens_ster ),
        # 
        # event_death_28d_from_hypoxia_cens_ster = case_when(steroids_within_5d_hypoxia == 0 & !is.na(first_steroid_date) ~ 0,
        #                                                 TRUE ~ event_death_28d_from_hypoxia)) #%>%
  # filter(steroids_within_5d_hypoxia == "No steroids within 5d hypoxia",
  #        !is.na(first_steroid_date)) %>%
  # select(empi, empi, ed_adm_dt, o2_crit_time, steroids_within_5d, steroids_within_5d_hypoxia, first_steroid_date,
  #        days_from_hosp_28d, days_from_hypoxia_28d, days_from_hosp_28d_cens_ster:event_death_28d_from_hypoxia_cens_ster,
  #        event_death_28d_from_hosp, event_death_28d_from_hypoxia) %>%
  #  view()
  #        


outcome_exposure %>%
  filter(time_to_steroids <= 0) %>%
  select(empi, ed_adm_dt, first_steroid_date)


outcome_exposure %>%
  filter(steroids_within_5d == "No steroids within 5d", event_death_28d_from_hosp_cens_ster == 1,
         days_from_hosp_28d_cens_ster <= 5) %>%
  nrow()

outcome_exposure %>%
  filter(steroids_within_5d_hypoxia == "No steroids within 5d hypoxia", event_death_28d_from_hypoxia_cens_ster == 1,
         days_from_hypoxia_28d_cens_ster <= 5) %>%
  nrow()

outcome_exposure %>%
  count(time_to_steroids, steroids_within_5d, event_death_28d_from_hosp_cens_ster) %>%
  ggplot(aes(time_to_steroids, n, fill=factor(event_death_28d_from_hosp_cens_ster))) +
  geom_bar(stat="identity") +
  labs(fill = "Death (censored\nif given steroids\nafter 5 days)",
       x = "Days to first steroids administration") +
  theme_bw()


outcome_exposure %>%
  count(time_to_steroids_hypoxia, steroids_within_5d_hypoxia, event_death_28d_from_hypoxia) %>%
  drop_na(steroids_within_5d_hypoxia) %>%
  ggplot(aes(time_to_steroids_hypoxia, n)) +
  geom_bar(stat="identity") +
  labs(x = "Days to first steroids administration from hypoxia") +
  theme_bw() +
  # scale_fill_manual(values = c("black","brown1")) +
  # labs(fill = "Mortality")+
  ggtitle("Distribution of first steroids administration among patients who received steroids",
          "Hypoxic patients only, relative to time of hypoxia") +
  scale_x_continuous(breaks = seq(from = min(outcome_exposure$time_to_steroids_hypoxia, na.rm=T),
                                  to = max(outcome_exposure$time_to_steroids_hypoxia, na.rm=T), by = 5))

outcome_exposure %>%
  count(time_to_steroids, steroids_within_5d, event_death_28d_from_hosp) %>%
  ggplot(aes(time_to_steroids, n)) +
  geom_bar(stat="identity") +
  labs(x = "Days to first steroids administration from hospitalization") +
  theme_bw() +
  # scale_fill_manual(values = c("black","brown1")) +
  # labs(fill = "Mortality")+
  ggtitle("Distribution of first steroids administration among patients who received steroids",
          "Relative to time of hospitalization") +
  scale_x_continuous(breaks = seq(from = min(outcome_exposure$time_to_steroids, na.rm=T),
                                  to = max(outcome_exposure$time_to_steroids, na.rm=T), by = 5))

outcome_exposure %>%
  mutate(steroids_within_5d = case_when(steroids_within_5d == 1 ~ "Steroids within 5d of hospitalization",
                                        steroids_within_5d == 0 ~ "No steroids within 5d of hospitalization"
  )) %>%
  count(days_from_hosp_28d, steroids_within_5d, event_death_28d_from_hosp) %>%
  ggplot(aes(days_from_hosp_28d, n, fill=factor(event_death_28d_from_hosp))) +
  geom_bar(stat="identity") +
  labs(x = "Days to death or discharge") +
  theme_bw() +
  scale_fill_manual(values = c("black","brown1")) +
  labs(fill = "Mortality") +
  facet_wrap(~steroids_within_5d, scales="free")


outcome_exposure %>%
  mutate(steroids_within_5d_hypoxia = case_when(steroids_within_5d_hypoxia == 1 ~ "Steroids within 5d of hypoxia",
                                        steroids_within_5d_hypoxia == 0 ~ "No steroids within 5d of hypoxia"
  )) %>%
  drop_na(steroids_within_5d_hypoxia) %>%
  count(days_from_hypoxia_28d, steroids_within_5d_hypoxia, event_death_28d_from_hypoxia) %>%
  ggplot(aes(days_from_hypoxia_28d, n, fill=factor(event_death_28d_from_hypoxia))) +
  geom_bar(stat="identity") +
  labs(x = "Days to death or discharge") +
  theme_bw() +
  scale_fill_manual(values = c("black","brown1")) +
  labs(fill = "Mortality") +
  facet_wrap(~steroids_within_5d_hypoxia, scales="free")

write_rds(outcome_exposure, "data/derived/outcome_exposure.rds")
