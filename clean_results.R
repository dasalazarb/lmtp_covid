## Code created by Katherine Hoffman - kah2797@med.cornell.edu
library(lmtp) 
library(tidyverse)
library(survival)
library(gt)
library(wesanderson)
library(viridis)

set.seed(7)

simple_estimator_contrast <- function(result, ref){
    shift <- result[["theta"]]
    ref <- ref[["theta"]]
    theta <- shift - ref
    return(tibble(theta, shift, ref))
  }

dat_lmtp <- readr::read_rds(here::here("data/dat_final.rds")) 
  

date <- "2021-12-23"
file_list <- list.files(here::here("results",date))
unadj_file_list <- file_list[which(str_detect(file_list, "notv"))]
file_list <- file_list[which(!str_detect(file_list, "notv"))]
file_list <- file_list[which(str_detect(file_list, "full"))]
file_list <- file_list[which(parse_number(file_list) %in% 1:4)]

results_files <- map(file_list, ~read_rds(here::here("results",date,.x)))
name_file <- function(file, cohort){
 #  cohort <- ifelse(str_detect(file, "full"), "full", "hypoxic")
  shift <- str_remove(str_remove(file[["shift"]], "\\("), "\\)")
  paste0(file[["estimator"]],"_shift_",shift)#,"_",cohort)
}
names(results_files) <- map2(results_files, file_list, ~name_file(.x, .y))
names(results_files)

unadj_results_files <- map(unadj_file_list, ~read_rds(here::here("results",date,.x)))
names(unadj_results_files) <- map2(unadj_results_files, file_list, ~name_file(.x, .y))

results_files

time_points <- length(results_files[[1]]$weights_m[[1]])
folds <- length(results_files[[1]]$weights_m)
num_lrnrs <- results_files[[1]]$weights_m[[1]][1] %>% as.data.frame() %>% nrow()

tmle <- lmtp_contrast(results_files$TMLE_shift_mtp,
              ref = results_files$TMLE_shift_NULL
              )

sdr <- lmtp_contrast(results_files$SDR_shift_mtp,
              ref = results_files$SDR_shift_NULL)


cis <- tibble(estimator = c("TMLE","SDR"),
              shift_low = c(results_files$TMLE_shift_mtp$low, results_files$SDR_shift_mtp$low),
              shift_high = c(results_files$TMLE_shift_mtp$high, results_files$SDR_shift_mtp$high),
              ref_low = c(results_files$TMLE_shift_NULL$low, results_files$SDR_shift_NULL$low),
              ref_high = c(results_files$TMLE_shift_NULL$high, results_files$SDR_shift_NULL$high)) %>%
  mutate(across(shift_low:ref_high, ~round(.x, 3)))
              

bind_rows(tmle$vals, sdr$vals) %>%
  mutate(across(everything(), ~round(.x, 3)),
         p.value = case_when(p.value == 0 ~ "<0.001",
                             TRUE ~ as.character(p.value))) %>%
  mutate(estimator = c("TMLE","SDR")) %>%
  left_join(cis) %>%
  mutate(across(c(shift, ref, shift_low, shift_high, ref_low, ref_high), ~1-.x),
         across(c(theta, conf.low, conf.high), ~-.x),
         across(where(is.numeric), ~.x*100)) %>%
  select(estimator, ref, ref_low, ref_high, shift, shift_low, shift_high, theta, conf.high, conf.low, p.value) %>%
  gt() %>%
  cols_merge(c(theta, conf.high, conf.low), pattern = "{1} ({2}, {3})") %>%
  cols_merge(c(ref, ref_high, ref_low), pattern = "{1} ({2}, {3})") %>%
  cols_merge(c(shift, shift_high, shift_low), pattern = "{1} ({2}, {3})") %>%
  cols_label(ref = "No Intervention", shift = "Delay intubation",
             estimator = "Estimator", theta = "Estimated difference", p.value = "p-value") %>%
  tab_spanner("Estimated 14-day AKI Rate", columns = c(ref, shift)) 


bind_rows(tmle$vals, sdr$vals) %>%
  mutate(across(everything(), ~round(.x, 3)),
         p.value = case_when(p.value == 0 ~ "<0.001",
                             TRUE ~ as.character(p.value))) %>%
  mutate(estimator = c("TMLE","SDR")) %>%
  # select(estimator, theta, shift, ref, conf.low, conf.high) %>%
  mutate(across(c(shift, ref), ~1-.x),
         across(c(theta, conf.low, conf.high), ~-.x),
         across(where(is.numeric), ~.x*100)) %>%
  select(estimator, ref, shift, theta, conf.high, conf.low, p.value) %>%
  gt() %>%
  cols_merge(c(theta, conf.high, conf.low), pattern = "{1} ({2}, {3})") %>%
  cols_label(ref = "No Intervention", shift = "Delay intubation",
             estimator = "Estimator", theta = "Estimated difference", p.value = "p-value") %>%
  tab_spanner("Estimated 14-day AKI Rate", columns = c(ref, shift)) %>%
  as_latex()

# find high density ratios ------------------------------------------------------

get_drs <- function(result_file){
  dr <- as.data.frame(result_file[["density_ratios"]]) 
  colnames(dr) <- paste0("density_ratio_",0:26)
  dr_merg <- bind_cols(dat_lmtp, dr)
  return(dr_merg)
}

extract_high_drs <- function(result_file, cutoff){
  get_drs(result_file) %>%
    select(id, starts_with("density_ratio")) %>%
    pivot_longer( starts_with("density_ratio")) %>%
    group_by(id) %>%
    mutate(max_density_ratio = max(value)) %>%
    mutate(day = parse_number(name)) %>%
    filter(max_density_ratio > cutoff) 
}

plot_high_drs <- function(result_file, cutoff){
  extract_high_drs(result_file, cutoff) %>%
    ggplot(aes(day, value, col=id)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(x = "time", y="density_ratio")
}

# why is the frequency of zeros getting lower then heigher?
get_drs(results_files$TMLE_shift_NULL) %>%
  select(id, starts_with("density_ratio")) %>%
  pivot_longer(starts_with("density_ratio")) %>%
  filter(value == 0) %>%
  count(name) %>%
  mutate(time = parse_number(name),
         freq = n / nrow(dat_lmtp)) %>%
  ggplot(aes(time, freq)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,1,by=.1))


# why is the frequency of zeros getting lower then heigher?
get_drs(results_files$TMLE_shift_NULL) %>%
  select(id, starts_with("density_ratio")) %>%
  pivot_longer(starts_with("density_ratio")) %>%
  filter(value > 0) %>%
  count(name) %>%
  mutate(time = parse_number(name)#,
         #freq = n / nrow(dat_lmtp)
         ) %>%
  ggplot(aes(time, n, label=n)) +
  geom_bar(stat="identity") +
  theme_bw() +
  geom_label(nudge_y = 200) +
  labs(x = "day", y = "n", title = "Non-zero density ratios under null shift")
  #scale_y_continuous(breaks=seq(0,1,by=.1))

ggsave("nonzero_rats.pdf", width=15, height=4)

fit <- survfit(Surv(fu, event) ~ 1, data=dat_lmtp)
library(survminer)
ggsurvplot(fit, xlim=c(0,28), break.x.by = 2)

int_times <- dat_lmtp %>%
  select(id, starts_with("I_")) %>% 
  pivot_longer(starts_with("I_")) %>%
  filter(value == 2) %>%
  distinct(id, .keep_all = T) %>%
  mutate(time_to_int = parse_number(name))

dat_lmtp %>% 
  select(id, event, fu) %>%
  left_join(int_times) %>%
  drop_na(time_to_int) %>%
  mutate(diff_aki_int = fu - time_to_int) %>%
  #count(diff_aki_int)
  ggplot(aes(diff_aki_int)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=seq(0,28,by=2), limits=c(0,28)) +
  theme_bw() +
  labs(x = "Days to AKI after intubation", y = "Number of patients")

get_drs(results_files$TMLE_shift_mtp_steroids_after_hypoxia) %>%
  select(id, starts_with("density_ratio")) %>%
  pivot_longer( starts_with("density_ratio")) %>%
  filter(value < .01) %>%
  count(name) %>%
  mutate(time = parse_number(name),
         freq_null_shift = n / nrow(dat_lmtp)) %>%
  select(time, starts_with("freq"))


get_drs(results_files$TMLE_shift_mtp) %>%
  select(id, starts_with("density_ratio")) %>%
  pivot_longer( starts_with("density_ratio")) %>%
  mutate(time = parse_number(name)) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~time, scales="free") +
  labs(title = "Distribution of raw ratios -- Intubation policy") +
  theme_bw()

rawratio_plots <- imap(results_files, function(x,y){
  get_drs(x) %>%
    select(id, starts_with("density_ratio")) %>%
    pivot_longer( starts_with("density_ratio")) %>%
    mutate(time = parse_number(name)) %>%
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~time, scales="free") +
    labs(title = paste0("Distribution of raw ratios -- ", y)) +
    theme_bw()
})


dt <- "2021-12-02"
imap(rawratio_plots, function(x,y){
  ggsave(filename = here::here(paste0("plots/rawratio_",y,".pdf")), x)
})

# look at weights ------------------------------------------------------

weights_list <-
  as.list(paste0(
    rep(c("TMLE","SDR"), each = 2),
    "_shift_",
          rep(c("mtp", "NULL"), 2)))
names(weights_list) <- unlist(weights_list)

n_est <- results_files$SDR_shift_mtp$weights_m[[1]] %>% length()
mygrid <- expand.grid(1:folds, 1:n_est)

# wrangle weights
weights_by_time <- function(estimator_res, weight, max_time, folds, n_est) {
  weights <- paste0("weights_",weight)
  temp_weights <- estimator_res[[weights]]
  # a few models are completely dropped, fix this
  fixed_weights <- map2(mygrid$Var1, mygrid$Var2, function(x,y){
    out <- temp_weights[[x]][[y]]
    if(nrow(out) != num_lrnrs){
      miss_lrnrs <- num_lrnrs - nrow(out)
      out <- bind_rows(out, data.frame(lrnrs = rep(NA, miss_lrnrs), weights = rep(NA, miss_lrnrs)))
    }
    return(out)
  })
  fixed_weights %>%
    bind_rows() %>%
    mutate(time = rep(0:(max_time-1), each = folds*n_est),
           fold = rep(1:folds, each = max_time*n_est)) %>%
    group_by(time, lrnrs) %>%
    summarise(mean_weight = mean(weights))
}

weights_by_time(results_files[[1]], "m", max_time = time_points, folds = folds, n_est = num_lrnrs)
weights_by_time(results_files[[4]], "m", max_time = time_points, folds = folds, n_est = num_lrnrs)


weights_by_time(results_files[[1]], "m", max_time = time_points, folds = folds, n_est = num_lrnrs)%>%
  ggplot(aes(time, mean_weight)) +
  geom_point(size=.8) +
  facet_wrap(~lrnrs) +
  theme_bw() 


pal <- c("#7a95b0", "#9e0202")

mplots <- imap_dfr(results_files, function(x,y){
  weights_by_time(x, "m", max_time = time_points, folds = folds, n_est = num_lrnrs) %>%
    drop_na(lrnrs) %>%
    mutate(est_shift = y)
}
) %>%
  group_by(lrnrs) %>%
  mutate(any_weight_.1 = case_when(mean_weight > .15 ~ 1)) %>%
  fill(any_weight_.1, .direction = "downup") %>%
  ungroup() %>%
  filter(any_weight_.1 == 1)  %>%
  mutate(lrnrs = case_when(lrnrs == "Lrnr_bartMachine_50_250_FALSE_0.95_2_2_0.9_3_1000_0.5" ~ "BART",
                           lrnrs == "Lrnr_earth_2_1_backward_0_1_0_0" ~ "MARS 1",
                           lrnrs == "Lrnr_glmnet_NULL_deviance_10_0_100_TRUE" ~ "Ridge",
                           lrnrs == "Lrnr_glmnet_NULL_deviance_10_0.5_100_TRUE" ~ "Elastic Net",
                           lrnrs == "Lrnr_glmnet_NULL_deviance_10_1_100_TRUE" ~ "LASSO",
                           lrnrs == "Lrnr_mean" ~ "Mean"
                           ),
         lrnrs = fct_relevel(factor(lrnrs), "Mean","LASSO","Elastic Net", "Ridge"),
         estimator = case_when(str_detect(est_shift,"SDR") ~ "SDR",
                               str_detect(est_shift,"TMLE") ~ "TMLE"),
         shift = case_when(str_detect(est_shift,"NULL") ~ "No intervention",
                               str_detect(est_shift,"mtp") ~ "Delay intubation")
         ) %>%
  ggplot(aes(time, mean_weight, col=shift)) +
  geom_point(size=.8) +
  facet_grid(estimator~lrnrs) +
  theme_bw() +
  labs(x = "Study Day", y="Weight", col = "Treatment Policy",
       title = bquote(italic(m[t]))) +
  theme(text=element_text(family="Times", size=11),
        strip.background = element_blank(),
        #  axis.title = element_text(face="bold"),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA)) +
  # scale_color_manual(values=pal) +
  scale_color_nejm()
mplots
#ggsave("m_est_draft.pdf", mplots, width=10, height=5)

rplots <- imap_dfr(results_files, function(x,y){
  weights_by_time(x, "r", max_time = time_points, folds = folds, n_est = num_lrnrs) %>%
    drop_na(lrnrs) %>%
    mutate(est_shift = y)
}
) %>%
  group_by(lrnrs) %>%
  mutate(any_weight_.1 = case_when(mean_weight > .15 ~ 1)) %>%
  fill(any_weight_.1, .direction = "downup") %>%
  ungroup() %>%
  filter(any_weight_.1 == 1)  %>%
  mutate(lrnrs = case_when(lrnrs == "Lrnr_bartMachine_50_250_FALSE_0.95_2_2_0.9_3_1000_0.5" ~ "BART",
                           lrnrs == "Lrnr_earth_2_1_backward_0_1_0_0" ~ "MARS 1",
                           lrnrs == "Lrnr_earth_2_2_backward_0_1_0_0" ~ "MARS 2",
                           lrnrs == "Lrnr_earth_2_3_backward_0_1_0_0" ~ "MARS 3",
                           lrnrs == "Lrnr_earth_3_1_backward_0_1_0_0" ~ "MARS 4",
                           lrnrs == "Lrnr_earth_3_2_backward_0_1_0_0" ~ "MARS 5",
                           lrnrs == "Lrnr_glmnet_NULL_deviance_10_0.5_100_TRUE" ~ "Elastic Net",
                           lrnrs == "Lrnr_glmnet_NULL_deviance_10_1_100_TRUE" ~ "LASSO",
                           TRUE ~ lrnrs
                          # lrnrs == "Lrnr_mean" ~ "Mean"
  ),
  lrnrs = fct_relevel(factor(lrnrs), "LASSO","Elastic Net"),
  estimator = case_when(str_detect(est_shift,"SDR") ~ "SDR",
                        str_detect(est_shift,"TMLE") ~ "TMLE"),
  shift = case_when(str_detect(est_shift,"NULL") ~ "No intervention",
                    str_detect(est_shift,"mtp") ~ "Delay intubation")
  ) %>%
  filter(estimator == "TMLE") %>%
  ggplot(aes(time, mean_weight, col=shift)) +
  geom_point(size=.8) +
  facet_wrap(~lrnrs, nrow=1) +
  theme_bw() +
  labs(x = "Study Day", y="Weight", col = "Treatment Policy",
       title = bquote(italic(r[t]))) +
  theme(text=element_text(family="Times", size=11),
        strip.background = element_blank(),
        #  axis.title = element_text(face="bold"),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA)) +
 # scale_color_manual(values=pal) +
  scale_color_nejm()
rplots

#ggsave("r_est_draft.pdf", rplots, width=12, height=5)



library(patchwork)
totalplots <- rplots + mplots + plot_layout(nrow = 2, guides = "collect", height=c(1,2)) +
  plot_annotation(tag_levels = "A", ) &
  theme(legend.position='bottom', plot.tag = element_text(face = "bold")) 
totalplots
ggsave("output/figure_weights.pdf", totalplots, width=6.7, height=5.4)



# SAVE r and m plots
dt <- "10132021"
dir.create(file.path(here::here("output/ivan_mtgs/"), dt), showWarnings = FALSE)

imap(mplots, function(x,y){
  ggsave(filename = here::here(paste0("output/ivan_mtgs/", dt,"/mweights_",y,".pdf")), x)
})

imap(rplots, function(x,y){
  ggsave(filename = here::here(paste0("output/ivan_mtgs/", dt, "/rweights_",y,".pdf")), x)
})

eif_df <- data.frame("tmle_eif" = results_files[[2]]$eif,
                      "sdr_eif" = results_files[[5]]$eif)

eif_df %>%
  ggplot(aes(tmle_eif, sdr_eif)) +
  geom_point() +
  theme_bw()

plot_eifs <- function(rf1, rf2){
  
}

tmle_m <- results_files$TMLE_shift_mtp$outcome_reg %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>%
  mutate(time = as.numeric(str_remove(name, "V"))) %>%
  select(time, tmle_value = value)

sdr_m <- results_files$SDR_shift_mtp$outcome_reg %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>%
  mutate(time = as.numeric(str_remove(name, "V"))) %>%
  select(sdr_value = value)

ms <- bind_cols(tmle_m, sdr_m)

ms %>%
  filter(time != 15) %>%
  ggplot(aes(tmle_value, sdr_value)) +
  geom_point(size=.5) +
  theme_bw() +
  facet_wrap(~time)

# making a table 1 ------------------------------------------------------
library(gtsummary)
cohort <- read_rds(here::here("data/derived/hospitalized_cohort_with_queens.rds"))

o2_time <- read_rds(here::here("data/derived/o2_time.rds"))
o2_time <- o2_time %>%
  mutate(vent = case_when(first_g6l_type %in% c("Ventilator", "Intubated via ETT
") ~ "Severe hypoxia with immediate IMV", TRUE ~ "Severe hypoxia without IMV initially"))

cohort %>% 
  left_join(o2_time) %>%
  select(age, sex, race, dm, cad, hiv, hf, copd, hepatitis, ckd_or_esrd, vent,
         intubation1) %>%
  mutate(vent = case_when(is.na(vent) ~ "Never severely hypoxic",
                          T ~ vent)) %>%
  labelled::set_variable_labels(
    age = "Age", sex = "Sex", race = "Race",
    dm = "DM", cad = "CAD", hiv = "HIV", copd = "COPD", hf = "HF",
    hepatitis = "Hepatitis", ckd_or_esrd = "Renal disease",
    vent = "Hypoxia criteria", intubation1 = "Ever intubated"
  ) %>%
  tbl_summary() %>%
  bold_labels()


cohort %>% 
  left_join(o2_time) %>%
  select(age, sex, race, dm, cad, hiv, hf, copd, hepatitis, ckd_or_esrd, vent,
         intubation1) %>%
  mutate(vent = case_when(is.na(vent) ~ "Never severely hypoxic",
                          T ~ vent)) %>%
  labelled::set_variable_labels(
    age = "Age", sex = "Sex", race = "Race",
    dm = "DM", cad = "CAD", hiv = "HIV", copd = "COPD", hf = "HF",
    hepatitis = "Hepatitis", ckd_or_esrd = "Renal disease",
    vent = "Hypoxia criteria", intubation1 = "Ever intubated"
  ) %>%
  tbl_summary(by = vent) %>%
  bold_labels()


