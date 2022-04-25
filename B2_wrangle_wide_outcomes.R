## Code created by Katherine Hoffman - kah2797@med.cornell.edu
library(tidyverse)
library(lubridate)

## load data

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds"))

outcome_exposure <- read_rds(here::here("data/derived/outcome_exposure.rds"))%>%
  filter(empi %in% cohort$empi)

aki <- read_rds(here::here("data/derived/kh_aki.rds"))
tz(aki$aki_time) <- "America/New_York"

events_day_with_steroids <- readr::read_rds(here::here("data/derived/events_day_w_steroids.rds")) %>%
  filter(empi %in% cohort$empi)

dat_full <-  readr::read_rds(here::here("data/derived/dat_full.rds"))


outcomes <- 
  cohort %>% 
  left_join(aki) %>% 
  drop_na(max_aki_stage) %>% # missing data lake data
  select(empi, ed_adm_dt, aki_time, max_aki_stage, end_dt, death) %>% 
  # add days from hospitalization to death
  mutate(event_aki = case_when(max_aki_stage > 0 & (aki_time %within% interval(ed_adm_dt, end_dt)) ~ 1,
                                    TRUE ~ 0),
         # add in AKI time to the end date (should be before death)
         end_dt = pmin(end_dt, aki_time, na.rm=T),
         days_to_aki_death_discharge = ceiling(time_length(difftime(end_dt, ed_adm_dt), unit="day")),
         event_aki_28d_from_hosp = case_when(days_to_aki_death_discharge <= 28 & event_aki == 1 ~ 1,
                                               TRUE ~ 0),
         event_death_28d_from_hosp = case_when(days_to_aki_death_discharge <= 28 & death == "Yes" &
                                                 event_aki_28d_from_hosp == 0 ~ 1,
                                             TRUE ~ 0),
         # follow up time is censored at discharge
         days_to_aki_death_discharge_28d = case_when(days_to_aki_death_discharge <= 28 ~ days_to_aki_death_discharge,
                                        TRUE ~ 28)) %>% 
  select(id = empi,
         fu = days_to_aki_death_discharge,
         event = event_aki_28d_from_hosp,
         cr = event_death_28d_from_hosp) %>%
  filter(fu > 0) # note that final cohort is 3,300

max_fu_day <- 28
intubation <-
  outcomes %>%
  select(id, fu, event, cr) %>% 
  left_join(events_day_with_steroids %>% select(id = empi, day, o2_status)) %>%
  filter(day <= fu, day < max_fu_day) %>%
  filter(!(day == fu & (event == 1 | cr == 1))) %>%
  mutate(I = case_when(o2_status == "Intubated" ~ 2,
                       o2_status == "Supp O2" ~ 1,
                       o2_status == "None" ~ 0)) %>%
  mutate(day = str_pad(day, 2, pad="0")) %>% # make the cens match Nick's examples
  pivot_wider(id_cols = c("id","fu"), 
              names_from = day,
              values_from = I, 
              names_prefix = "I_")

max_fu_day <- 28
padded_days <- str_pad(0:(max_fu_day-1),2,pad="0")

replace_obs <- as.list(rep(0,max_fu_day)) # replace with all 0s
names(replace_obs) <- paste0("C_",padded_days)


cens <-
  outcomes %>%
  group_by(id) %>%
  mutate(day = fu) %>%
  complete(day = full_seq(1:day, 1)) %>%
  fill(fu, event, cr, .direction = "up") %>%
  mutate(C = case_when(event == 0 & cr == 0 & day == fu ~ 0, # we don't have data after patients are discharged, so back their censoring time up one day
                       TRUE ~ 1)) %>%
  mutate(day = str_pad(day - 1, 2, pad="0")) %>% # make the cens match Nick's examples
  pivot_wider(id_cols = c("id"), 
              names_from = day,
              values_from = C, 
              names_prefix = "C_") %>%
  replace_na(replace_obs) %>%
  ungroup() %>%
  select(id, one_of(names(replace_obs)))

# outcome wide format
outcome <-
  outcomes %>%
  group_by(id) %>%
  mutate(day = 1) %>%
  complete(day = full_seq(1:max_fu_day, 1)) %>%
  fill(fu, event, cr, .direction = "downup") %>%
  mutate(
    Y = case_when(
      cr == 1 ~ 0,
      day < fu ~ 0, # patients who died should not have the event before their last follow up day
      event == 0 & day == max_fu_day ~ 0, # patients who are censored on final day have an outcome observed
      event == 0 & day >= fu ~ NA_real_, # as soon as patients are censored, they should have NAs for outcome
      event == 1 & day >= fu ~ 1  # once a patient has an event, they always have an outcome status of 1
    )) %>%
  mutate(day = str_pad(day,2,pad="0")) %>%
  pivot_wider(id_cols = c("id","event"), 
              names_from = day,
              values_from = Y, 
              names_prefix = "Y_") %>%
  ungroup()


# outcome wide format
comp_risk <-
  outcomes %>%
  group_by(id) %>%
  mutate(day = 1) %>%
  complete(day = full_seq(1:max_fu_day, 1)) %>%
  fill(fu, event, cr, .direction = "downup") %>%
  mutate(
    cr = case_when(
      day < fu ~ 0, # patients who died should not have the comp risk before their last follow up day
      cr == 0 & day == max_fu_day ~ 0, # patients who are censored on final day have a comp risk of 0 observed
      cr == 0 & day >= fu ~ NA_real_, # as soon as patients are censored or have outcome, they should have NAs for comp risk
      cr == 1 & day >= fu ~ 1  # once a patient has competing risk, they always have an outcome status of 1
    )) %>%
  mutate(day = str_pad(day,2,pad="0")) %>%
  pivot_wider(names_from = "day",
              values_from = "cr", 
              names_prefix = "CR_") %>%
  ungroup()

i <- sample(dim(outcome)[1],size = 1);outcomes[i,];outcome[i,1:15] %>% select(-event);cens[i,1:14];intubation[i,1:15] %>% select(-fu);comp_risk[i,1:15] %>% select(-fu)

for (i in 1:dim(outcome)[1]) {
  if (sum(is.na(outcome[i,1:15] %>% select(-event))) == 0 & sum(is.na(intubation[i,1:15] %>% select(-fu))) > 0) {
    print(outcomes[i,])
    print(outcome[i,1:15] %>% select(-event))
    print(cens[i,1:14])
    print(intubation[i,1:15] %>% select(-fu))
    print(comp_risk[i,1:15] %>% select(-fu))
  }
  
}

dat_final <- 
  outcomes %>% # contains fu, event, cr columns
  full_join(cens) %>%
  full_join(outcome) %>%
  left_join(comp_risk) %>%
  left_join(intubation) %>%
  left_join(dat_full %>% rename(id = empi))  


for (i in 0:26){
  ip <- str_pad(i, 2, pad="0")
  ip1 <- str_pad(i+1, 2, pad="0")
  vars_for_na <- dat_final %>%
    select(ends_with(ip1),
           -paste0("CR_",ip1),
           -paste0("Y_",ip1),
           -paste0("C_",ip1)) %>%
    names()
  dat_final <-
    dat_final %>%
    mutate(across(vars_for_na, ~case_when(.data[[paste0("C_", ip)]] == 0 ~ NA_real_,
                                          .data[[paste0("Y_", ip1)]] == 1 ~ NA_real_,
                                          .data[[paste0("CR_", ip1)]] == 1 ~ NA_real_,
                                          TRUE ~ .x))) 
}

saveRDS(dat_final, "data/derived/dat_final.rds")


