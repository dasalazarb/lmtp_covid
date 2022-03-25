
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

set.seed(7)

my_letters <- function(length.out) {
  a <- rep(letters, length.out = length.out)
  grp <- cumsum(a == "a")
  vapply(seq_along(a), 
         function(x) paste(rep(a[x], grp[x]), collapse = ""),
         character(1L))
}
my_letters <- my_letters(60)

## ---------------------------

# CHANGE DEPENDING ON EVENT
follow_up_var <- "days_from_hosp_28d"
event_var <- "event_death_28d_from_hosp"

max_fu_day <- 28
padded_days <- str_pad(0:(max_fu_day-1),2,pad="0")

## load data

# read in all cohort after determining in steroids_create_analytic_file-1.Rmd
cohort <- read_rds(here::here("data/derived/cohort_w_labs.rds"))

steroids_over1mg_days <- read_rds(here::here("data/derived/steroids_over_1mg_days.rds"))%>%
  filter(empi %in% cohort$empi)

df_days <- read_rds(here::here("data/derived/df_days.rds")) %>%
  filter(empi %in% cohort$empi)

labs_long <- read_rds(here::here("data/derived/labs_long.rds")) %>%
  filter(!str_detect(lab, "troponin|fibrinogen|creatine|lowest_neutrophil|lowest_absolute_neutrophil|lowest_glucose|highest_platelet"))

vitals_long <- read_rds(here::here("data/derived/vitals_long.rds"))

max_days <- df_days %>% distinct(empi, max_day)

# impute W covariates
dat_wide_adj_temp <-
  cohort %>%
  select(empi,
         red_cap_source,
         age, sex, race, ethnicity, bmi, smoking, 
         cad, home_o2_yn, dm, htn, cva, cirrhosis, ckd_or_esrd, asthma, copd, active_cancer, 
         immunosuppressed, ild, hiv, hypoxia_ed, hypoxia_ed_method
  ) %>%
  mutate(bmi_miss = case_when(is.na(bmi) ~ 1, TRUE ~ 0),
         home_o2_miss = case_when(is.na(home_o2_yn) ~ 1, TRUE ~ 0), 
         hypoxia_ed = case_when(hypoxia_ed == "No" ~ 0, hypoxia_ed == "Yes" ~ 1),
         race_miss = case_when(is.na(race) ~ 1, TRUE ~ 0),
         ethnicity_miss = case_when(is.na(ethnicity) ~ 1, TRUE ~ 0),
         race = case_when(is.na(race) ~ "Missing", TRUE ~ race),
         ethnicity = case_when(is.na(ethnicity) ~ "Missing", TRUE ~ ethnicity),
         smoking = case_when(is.na(smoking) ~ "No", TRUE ~ smoking),
         #steroids_within_5d = case_when(steroids_within_5d == "Steroids within 5d" ~ 1,
         #                              TRUE ~ 0),
         sex = case_when(sex=="Male"~1, TRUE ~0),
         red_cap_source = case_when(red_cap_source == "QUEENS" ~ 0,
                                    red_cap_source == "EAST" ~ 1),
         across(cad:hiv, ~case_when(.x == "No" ~ 0, TRUE ~ 1))) %>%
  fastDummies::dummy_columns(select_columns = c("race", "ethnicity","hypoxia_ed_method","smoking")) %>%
  clean_names() %>%
  select(-race, -ethnicity, -hypoxia_ed_method, -smoking)


dat_wide_adj_mice <- mice::mice(dat_wide_adj_temp %>% select(-empi))
dat_wide_adj <- mice::complete(dat_wide_adj_mice) %>%
  bind_cols(dat_wide_adj_temp %>% select(empi))
names(dat_wide_adj)


# wide day 1 labs for mice
labs_temp <-
  max_days %>%
  left_join(labs_long) %>%
  expand(empi, lab, day) %>%
  filter(day == 0) %>%
  left_join(labs_long %>% select(empi, day, lab, value)) %>%
  drop_na(lab) %>%
  pivot_wider(names_from = lab, values_from = value, names_prefix = "lab_") %>%
  select(-day)

# wide day 1 vs for mice
vs_temp <-
  max_days %>%
  left_join(vitals_long) %>%
  expand(empi, vs, day) %>%
  filter(day == 0) %>%
  left_join(vitals_long %>% select(empi, day, vs, value))  %>%
  drop_na(vs) %>%
  pivot_wider(names_from = vs, values_from = value, names_prefix = "vs_")%>%
  select(-day)

day0_for_mice <-
  full_join(dat_wide_adj, labs_temp) %>%
  full_join(vs_temp) 

day0_mice_fit <- mice::mice(day0_for_mice %>% select(-empi))
day0_mice <- mice::complete(day0_mice_fit) %>%
  bind_cols(day0_for_mice %>% select(empi))
names(day0_mice)

# get mice values for day 1 back in long format
day0_labs_long <- 
  day0_mice %>%
  select(empi, starts_with("lab_")) %>%
  pivot_longer(cols = starts_with("lab")) %>%
  mutate(lab = str_remove_all(name, "lab_"), day = 0) %>%
  select(empi, day, lab, value)

# get mice values for day 1 back in long format
day0_vs_long <- 
  day0_mice %>%
  select(empi, starts_with("vs_")) %>%
  pivot_longer(cols = starts_with("vs")) %>%
  mutate(vs = str_remove_all(name, "vs_"), day = 0) %>%
  select(empi, day, vs, value)

# merge the mice values with the original values
labs_long_mod <-
  labs_long %>%
  select(empi, day, lab, value) %>%
  filter(day != 0) %>%
  full_join(day0_labs_long)

# merge the mice values with the original values
vitals_long_mod <-
  vitals_long %>%
  select(empi, day, vs, value) %>%
  filter(day != 0) %>%
  full_join(day0_vs_long)

labs_long_locf <-
  labs_long_mod %>%
  expand(empi, lab, day) %>% # get all days 
  left_join(labs_long_mod %>% select(empi, day, lab, value)) %>% # re-add values after expanding (col not kept)
  group_by(empi, lab) %>%
  fill(value, .direction = "down") %>% # carry observations forward
  full_join(max_days) %>%
  filter(day <= max_day) %>% # keep only values before max day (goes all the way to day =60)
  filter(day <= 27) # for now only looking at first 28 days
  
vitals_long_locf <-
  vitals_long_mod %>%
  expand(empi, vs, day) %>% # get all days 
  left_join(vitals_long_mod %>% select(empi, day, vs, value)) %>% # re-add values after expanding (col not kept)
  group_by(empi, vs) %>%
  fill(value, .direction = "down") %>% # carry observations forward
  full_join(max_days) %>%
  filter(day <= max_day) %>% # keep only values before max day (goes all the way to day =60)
  filter(day <= 27) # for now only looking at first 28 days
  
steroids_wide <-
  steroids_over1mg_days  %>%
  mutate(steroids_yn = 1) %>% # everyone got steroids
  select(empi, day, steroids_yn, max_day) %>%
  full_join(df_days %>% select(empi, day, max_day)) %>%
  group_by(empi) %>%
  complete(empi, day = full_seq(0:max_day, 1)) %>%
  filter(day <= 27) %>% # for now only looking at first 28 days
  mutate(steroids_yn = case_when(is.na(steroids_yn) ~ 0,
                                 TRUE ~ steroids_yn),
         day = str_pad(day, 2, pad="0")) %>%
  pivot_wider(id_cols = empi, names_from = day, values_from = steroids_yn,
              names_prefix = "A_") %>%
  select(empi, sort(names(.)))

# make letter alt keys for labs and vitals so easier to read the full wide data set
lab_name_key <-
  tibble(lab = unique(labs_long$lab),
         alt_id = letters[1:length(unique(labs_long$lab))])

vs_name_key <-
  tibble(vs = unique(vitals_long$vs),
         alt_id = my_letters[(length(unique(labs_long$lab))+1):(length(unique(labs_long$lab))+length(unique(vitals_long$vs)))])

### LOCF VERISONS OF LABS AND VITALS
# make labs and vs in L_a_1, etc. format
labs_widest_locf <-
  labs_long_locf %>%
  full_join(lab_name_key) %>%
  mutate(day = str_pad(as.character(day),2,pad="0")) %>%
  mutate(lab_day = paste0(alt_id, "_", day)) %>%
  pivot_wider(id_cols = empi, names_from = lab_day,
              values_from = value, names_prefix = "L_") %>% # imputed version of value
  select(sort(names(.))) 

# vitals_long %>% filter(day <= 28) %>% distinct(empi, day) ## 16,538
vitals_widest_locf <-
  vitals_long_locf %>%
  full_join(vs_name_key) %>%
  mutate(day = str_pad(as.character(day),2,pad="0")) %>%
  mutate(vs_day = paste0(alt_id, "_", day)) %>%
  pivot_wider(id_cols = empi, names_from = vs_day,
              values_from = value, names_prefix = "L_") %>%
  select(sort(names(.))) 

replace_miss <- as.list(rep(0,max_fu_day)) # replace with all 0s
names(replace_miss) <- paste0("L_miss_",padded_days)

# indicator variable for whether d-dimer was measured
miss_labs_wide <-
  labs_long %>%
  filter(lab == "highest_d_dimer") %>% # only keep patients with a d-dimer
  select(empi, day) %>%
  filter(day <= 27) %>%
  mutate(value = 1, # change their "value" to 1 indicating they had a d-dimer
         day = str_pad(day, 2, pad = "0")) %>%
  pivot_wider(names_from = day, names_prefix = "L_miss_")  %>% # wide format
  right_join(cohort %>% select(empi)) %>% # add all the other patients back in (that didn't have d-dimer)
  replace_na(replace_miss) %>% # replace with 0 (no d-dimer)
  ungroup() %>%
  select(sort(names(.)))


# make full data set LOCF VERSION
dat_full_locf <-
  dat_wide_adj %>%
  left_join(labs_widest_locf) %>%
  left_join(miss_labs_wide) %>%
  left_join(vitals_widest_locf) %>%
  left_join(steroids_wide) 


# SAVE DATA FILES
saveRDS(steroids_wide, here::here("data/derived/steroids_wide.rds"))
saveRDS(vitals_widest_locf, here::here("data/derived/vitals_widest_locf.rds"))
saveRDS(labs_widest_locf, here::here("data/derived/labs_widest_locf.rds"))
saveRDS(miss_labs_wide, here::here("data/derived/miss_labs_wide.rds"))
saveRDS(dat_wide_adj, here::here("data/derived/dat_wide_adj.rds"))
saveRDS(dat_full_locf, here::here("data/derived/dat_full.rds"))
saveRDS(vs_name_key, here::here("data/derived/vs_name_key.rds"))
saveRDS(lab_name_key, here::here("data/derived/lab_name_key.rds"))

