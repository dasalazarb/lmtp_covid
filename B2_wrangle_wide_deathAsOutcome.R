## Code created by Katherine Hoffman - kah2797@med.cornell.edu
library(tidyverse)
library(lubridate)

## load data

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds")); head(cohort)

outcome_exposure <- read_rds(here::here("data/derived/outcome_exposure.rds"))%>%
  filter(empi %in% cohort$empi); head(outcome_exposure)

events_day_with_steroids <- readr::read_rds(here::here("data/derived/events_day_w_steroids.rds")) %>%
  filter(empi %in% cohort$empi); head(events_day_with_steroids)

dat_full <-  readr::read_rds(here::here("data/derived/dat_full.rds")); head(dat_full)


outcomes <- 
  cohort %>% 
  dplyr::select(empi, ed_adm_dt, end_dt, death) %>% #aki_time, max_aki_stage,
  # add days from hospitalization to death
  dplyr::mutate(
         # add in AKI time to the end date (should be before death)
         days_to_death_or_discharge = ceiling(time_length(difftime(end_dt, ed_adm_dt), unit="day")),
         event_death_28d_from_hosp = case_when(days_to_death_or_discharge <= 28 & death == "Yes" ~ 1,
                                             TRUE ~ 0)
         ) %>% 
  dplyr::select(id = empi,
         fu = days_to_death_or_discharge,
         event = event_death_28d_from_hosp) %>%
  filter(fu > 0) %>% 
  distinct(); dim(outcomes) # note that final cohort is 3,300
outcomes

max_fu_day <- 28
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

# max_fu_day <- 28
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
outcome

# ## using the cens pattern to correct the intubation pattern
# intubation[,3:dim(intubation)[2]] <- (intubation[,3:dim(intubation)[2]] + 1) * cens[,2:dim(cens)[2]]
# 
# intubation[,3:dim(intubation)[2]][intubation[,3:dim(intubation)[2]] == 0] <- NA
# 
# intubation[,3:dim(intubation)[2]] <- intubation[,3:dim(intubation)[2]] - 1

i <- sample(which(outcome$event == 0),size = 1);outcomes[i,];outcome[i,1:15] %>% select(-event);cens[i,1:14];intubation[i,1:15] %>% select(-fu)

dat_final <- 
  outcomes %>% # contains fu, event, cr columns
  full_join(cens) %>% 
  full_join(outcome) %>% 
  left_join(intubation) %>% 
  left_join(dat_full %>% rename(id = empi))


for (i in 0:26){
  ip <- str_pad(i, 2, pad="0")
  ip1 <- str_pad(i+1, 2, pad="0")
  vars_for_na <- dat_final %>% 
    select(ends_with(ip1),
           -paste0("Y_",ip1),
           -paste0("C_",ip1)) %>% 
    names()
  dat_final <-
    dat_final %>%
    mutate(across(all_of(vars_for_na), ~ case_when(.data[[paste0("C_", ip)]] == 0 ~ NA_real_,
                                          .data[[paste0("Y_", ip1)]] == 1 ~ NA_real_,
                                          TRUE ~ .x))) 
}

saveRDS(dat_final, "data/derived/dat_final_deathAsOutcome.rds")

# write.csv(dat_final %>% select_if(function(x) any(is.na(x))), "review_Y_L.csv")
# 
dat_final <- readRDS("data/derived/dat_final_deathAsOutcome.rds")


cens <- dat_final %>%
  select(id, starts_with("C_"))

intubation <- dat_final %>%
  select(id, starts_with("I_"))

outcome <- dat_final %>%
  select(id, starts_with("Y_"))

outcomes <- dat_final %>%
  select(id, fu, event)

i <- sample(dim(dat_final)[1],size = 1);outcomes[i,];outcome[i,1:15];cens[i,1:15];intubation[i,1:15]
i <- sample(which(dat_final$fu == 1),size = 1);outcomes[i,];outcome[i,1:15];cens[i,1:15];intubation[i,1:15]

# out_mtp

# ### Without mpt
# progressr::with_progress(
#   debug(out_NULL <-
#     lmtp_sdr(
#       dat_lmtp,
#       trt = a,
#       outcome = y,
#       # comp_risk = cr,
#       baseline = bs,
#       time_vary = tv,
#       cens = censoring,
#       shift = NULL,
#       outcome_type = "survival",
#       learners_outcome = lrnrs,
#       # learners_trt = lrnrs,
#       folds = folds,
#       .SL_folds = SL_folds,
#       # .trim = trim,
#       k=k,
#       intervention_type = "mtp"
#     ))
# )
# # out_NULL
# # 
# ggsurvplot(
#   fit = survfit(Surv(fu, event) ~ 1, 
#                 data = dat_lmtp %>% mutate(fu = ifelse(fu >=15, 15, fu))), 
#   xlab = "Days", 
#   ylab = "Overall survival probability",
#   risk.table = TRUE, break.time.by = 1)
# fu_ <- dat_lmtp %>% filter(fu <= 15)
# table(fu_$fu, fu_$event)

# write.csv(x = dat_lmtp %>% 
#   select(fu, event, starts_with("C_")) %>% 
#   group_by(fu) %>% 
#   arrange(fu, event) %>% 
#   ungroup() %>% 
#   distinct() %>% 
#   filter(fu >=2), file = "review_C_.csv")
