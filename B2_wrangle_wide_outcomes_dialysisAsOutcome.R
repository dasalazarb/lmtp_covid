## Code created by Katherine Hoffman - kah2797@med.cornell.edu
## and edited by me
library(tidyverse)
library(lubridate)

## load data

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds")); head(cohort); dim(cohort)

outcome_exposure <- read_rds(here::here("data/derived/outcome_exposure.rds"))%>%
  filter(empi %in% cohort$empi); head(outcome_exposure); dim(outcome_exposure)

## old dialysis data
dialysis.1 <- read_rds(here::here("data/raw/2020-08-17/covid_datalake_wcm_procedures.rds")); head(dialysis.1); dim(dialysis.1)

## new dialysis data
dialysis.2 <- readxl::read_excel(here::here("data/raw/2020-08-17/RITM0529623_hoffman.xlsx")) %>% 
  mutate(procedure_description = "dialysis", procedure_dt_tm=as.Date(dialysis_date)) %>% 
  select(-dialysis_date); head(dialysis.2); dim(dialysis.2)

# tz(dialysis.1$procedure_dt_tm) <- "America/New_York"
# tz(dialysis.2$procedure_dt_tm) <- "America/New_York"

## merge old and new data
dialysis <- dialysis.2 %>% 
  mutate(empi=as.character(empi)) %>% 
  ## bind rows of old data
  bind_rows(dialysis.1 %>% 
              filter(str_detect(tolower(procedure_description), "dialysis")) %>% 
              mutate(procedure_description = "dialysis") %>% 
              select(empi, procedure_dt_tm, procedure_description) ) %>% 
  ## remove dates before 2020-03-03 -> useful for the slice_min
  filter(procedure_dt_tm >= "2020-03-03") %>% 
  ## when was the first dialysis procedure for each patient?
  group_by(empi) %>% 
  slice_min(procedure_dt_tm) %>% 
  ungroup() %>% 
  ## maintain unique records per patient
  distinct()

## there are 232 patients with dialysis procedure
head(dialysis, 5); dim(dialysis)

## between 2020-03-03 and 2020-07-01, there are 201 patients
sum(dialysis$procedure_dt_tm >= as.Date("2020-03-03") & dialysis$procedure_dt_tm <= as.Date("2020-07-01"))

## However, # common patients in cohort and dialysis dataset
length(intersect(cohort$empi, dialysis$empi)) ## 52

# tz(dialysis$procedure_dt_tm) <- "America/New_York"
# tz(dialysis$dialysis_date) <- "America/New_York"

events_day_with_steroids <- readr::read_rds(here::here("data/derived/events_day_w_steroids.rds")) %>%
  filter(empi %in% cohort$empi)

dat_full <-  readr::read_rds(here::here("data/derived/dat_full.rds"))


outcomes <- 
  cohort %>% 
  left_join(dialysis, by="empi") %>% # -> 52 patients with dialysis info
  select(empi, ed_adm_dt, end_dt, procedure_dt_tm, procedure_description, death) %>% 
  # add days from hospitalization to death
  mutate(event_dialysis = case_when(str_detect(tolower(procedure_description), "dialysis") ~ 1, 
                                    TRUE ~ 0), ## -> 14 event_dialysis
         # add in dialysis time to the end date (should be before death)
         end_dt = pmin(end_dt, procedure_dt_tm, na.rm=T),
         days_to_dialysis_death_discharge = ceiling(time_length(difftime(end_dt, ed_adm_dt), unit="day")),
         event_dialysis_28d_from_hosp = case_when(days_to_dialysis_death_discharge <= 28 & event_dialysis == 1 ~ 1,
                                               TRUE ~ 0),
         event_death_28d_from_hosp = case_when(days_to_dialysis_death_discharge <= 28 & death == "Yes" &
                                                 event_dialysis_28d_from_hosp == 0 ~ 1,
                                               TRUE ~ 0),
         # follow up time is censored at discharge
         days_to_dialysis_death_discharge_28d = case_when(days_to_dialysis_death_discharge <= 28 ~ days_to_dialysis_death_discharge,
                                        TRUE ~ 28)) %>% 
  select(id = empi,
         fu = days_to_dialysis_death_discharge,
         event = event_dialysis_28d_from_hosp,
         cr = event_death_28d_from_hosp) %>% 
  filter(fu > 0) %>% 
  distinct(); outcomes; dim(outcomes) # note that final cohort is 3,300
sum(table(outcomes$fu, outcomes$event)[,2])

max_fu_day <- 28

padded_days <- str_pad(0:(max_fu_day-1),2,pad="0")

replace_obs <- as.list(rep(0,max_fu_day)) # replace with all 0s
names(replace_obs) <- paste0("C_",padded_days)


cens <- 
  outcomes %>%
  group_by(id) %>%
  mutate(day = fu) %>%
  complete(day = tidyr::full_seq(1:day, 1)) %>% 
  tidyr::fill(fu, event, .direction = "up") %>% 
  mutate(C = case_when(event == 0 & day == fu ~ 0, # we don't have data after patients are discharged, so back their censoring time up one day
                       TRUE ~ 1)) %>% 
  mutate(day = str_pad(day - 1, 2, pad="0")) %>% # make the cens match Nick's examples
  pivot_wider(id_cols = c("id"), 
              names_from = day, 
              values_from = C, 
              names_prefix = "C_") %>% 
  replace_na(replace_obs) %>% 
  ungroup() %>% 
  select(id, one_of(names(replace_obs)))

intubation <-
  outcomes %>%
  select(id, fu, event) %>%
  left_join(events_day_with_steroids %>% select(id = empi, day, o2_status)) %>%
  filter(day <= fu, day < max_fu_day) %>%
  filter(!(day == fu & event == 1)) %>%
  mutate(I = case_when(o2_status == "Intubated" ~ 2,
                       o2_status == "Supp O2" ~ 1,
                       o2_status == "None" ~ 0)) %>%
  mutate(day = str_pad(day, 2, pad="0")) %>% # make the cens match Nick's examples
  pivot_wider(id_cols = c("id","fu"),
              names_from = day,
              values_from = I,
              names_prefix = "I_")

# i <- sample(dim(outcome)[1],size = 1);outcomes[i,];outcome[i,1:14];cens[i,1:14];intubation[i,1:14]

# outcome wide format
outcome <-
  outcomes %>%
  group_by(id) %>%
  mutate(day = 1) %>% 
  complete(day = full_seq(1:max_fu_day, 1)) %>% 
  fill(fu, event, .direction = "downup") %>% 
  mutate(
    Y = case_when(
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
  # replace(is.na(.), 0) %>% 
  ungroup()
outcome %>% 
  slice_sample(n = 10)


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
comp_risk %>% 
  slice_sample(n=10)


## using the cens pattern to correct the intubation pattern
# intubation[,3:dim(intubation)[2]] <- (intubation[,3:dim(intubation)[2]] + 1) * cens[,2:dim(cens)[2]]

# intubation[,3:dim(intubation)[2]][intubation[,3:dim(intubation)[2]] == 0] <- NA

# intubation[,3:dim(intubation)[2]] <- intubation[,3:dim(intubation)[2]] - 1

i <- sample(which(outcome$event == 0),size = 1);outcomes[i,];outcome[i,1:15] %>% select(-event);cens[i,1:14];intubation[i,1:15] %>% select(-fu);comp_risk[i,1:16] %>% select(-fu, -event)


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
    mutate(across(all_of(vars_for_na), ~case_when(.data[[paste0("C_", ip)]] == 0 ~ NA_real_,
                                          .data[[paste0("Y_", ip1)]] == 1 ~ NA_real_,
                                          .data[[paste0("CR_", ip1)]] == 1 ~ NA_real_,
                                          TRUE ~ .x))) 
}

table(dat_final %>% filter(fu<=14) %>% select(event))

saveRDS(dat_final, "data/derived/dat_final_deathAsCompRisk.rds")


