## Code created by Katherine Hoffman - kah2797@med.cornell.edu
options(java.parameters = "-Xmx20000m")

# remotes::install_local(here::here("lmtp_1.0.0.5001.tar.gz"))

library(sl3)
library(lmtp)
library(tidyverse)
library(future)
library(earth)

progressr::handlers(global = TRUE)
set.seed(7)

# task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# plan(multicore)

dat_lmtp <- read_rds(here::here("data/derived/dat_final_deathAsOutcome.rds")) %>%
  filter(hypoxia_ed == 1) %>% # our cohort is only people who were hypoxic initially
  mutate(I_00 = ifelse(I_00 == 0, 1, I_00)) %>% 
  replace(is.na(.), 0) %>% 
  select(-ckd_or_esrd, -hypoxia_ed, ## 1 class
         -ild, -hiv, ## few observation in one class
         -hypoxia_ed_method_na, -hypoxia_ed_method_none, ## few obs in one class
         -ethnicity_missing, ##ethnicity_missing == ethnicity_miss
         ## checking these ones:
         -cirrhosis,
         -smoking_active_smoker, -hypoxia_ed_method_venti_mask, -hypoxia_ed_method_niv_bipap_cpap, 
         -hypoxia_ed_method_high_flow_nasal_cannula); dim(dat_lmtp) 


# dat_slice <- cbind(dat_lmtp[,c(bs,"event")], dat_lmtp[,tv[[1]]])
# dat_slice.2 <- dat_slice %>%
#   select(#ckd_or_esrd, cirrhosis, ild, hiv, hypoxia_ed_method_high_flow_nasal_cannula, hypoxia_ed_method_niv_bipap_cpap
#     #hypoxia_ed_method_venti_mask, smoking_active_smoker
#          event, age, sex, bmi,red_cap_source, cad, home_o2_yn,
#         dm, htn, cva, asthma, copd, active_cancer, immunosuppressed,
#         bmi_miss, home_o2_miss, race_miss, ethnicity_miss, race_asian,
#         race_black, race_missing, race_other, race_white, ethnicity_hispanic_or_latino_or_spanish_origin,
#         ethnicity_not_hispanic_or_latino_or_spanish_origin,ethnicity_not_hispanic_or_latino_or_spanish_origin,
#         hypoxia_ed_method_mechanical_ventilation,
#         hypoxia_ed_method_nasal_cannula, hypoxia_ed_method_non_rebreather,
#         smoking_former_smoker, smoking_no); dim(dat_slice.2)
# 
# 
# for (i in bs) {
#   if (!(grepl("age", i) | grepl("bmi", i))) {
#     print(i); print(table(dat_slice[,i][[1]], dat_slice[,"event"]$event))
#   }
# }
# 
# fit <- glmnet::glmnet(x = as.matrix(dat_slice.2[, c("age", "sex", "bmi")]), y = as.matrix(dat_slice.2[,"event"]$event), family = binomial(), alpha = 1)
# 
# summary(glm(formula = event ~ ., family = binomial(), data = dat_slice.2))
# 
# cvfit <- glmnet::cv.glmnet(x = as.matrix(dat_slice[, !grepl("event", colnames(dat_slice))]),
#                            y = as.matrix(dat_slice[,"event"]), family = binomial(),
#                            alpha = 1, nfolds = 5); coef(cvfit, s = "lambda.min")

trim <- .995
folds <- 5
SL_folds <- 5
k <- 2

lrn_rf <- Lrnr_randomForest$new()
lrn_glmfast <- Lrnr_glm_fast$new()
lrn_lasso <- Lrnr_glmnet$new(alpha = 1, stratify_cv = TRUE)
lrn_glm <- Lrnr_glm$new(stratify_cv = TRUE)
lrn_mean <- Lrnr_mean$new()
lrn_bart <- Lrnr_bartMachine$new()
lrn_earth <- Lrnr_earth$new()
lrn_rpart <- Lrnr_rpart$new()
lrnr_lgb <- Lrnr_lightgbm$new()
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_enet <- Lrnr_glmnet$new(alpha = 0.5)

learners_simple <- unlist(list(
  # lrn_rf, 
  # lrn_glm,
  # lrn_lasso#,
  lrnr_lgb,
  lrn_mean
), recursive = TRUE)

lrnrs <- make_learner(Stack, learners_simple)

# set parameters of outcome, trt, and adjustment vars
outcome_day <- 14
padded_days <- str_pad(0:(outcome_day-1), 2, pad = "0")
padded_days_out <- str_pad(1:outcome_day, 2, pad = "0")

a <-  paste0("I_", padded_days) # treatment 
bs <- dat_lmtp %>% # baseline covariates
  select(-id, -fu, -event,
         -starts_with("L_"), -starts_with("C_"),
         -starts_with("Y_"), -starts_with("A_"),
         -starts_with("I_"), -starts_with("CR_"),
         -starts_with("H_")) %>% names()
y <- paste0("Y_",padded_days_out) # outcome (AKI)
# cr <- paste0("CR_",padded_days_out) # competing risk (death)
censoring <- paste0("C_",padded_days) # observed at next time

used_letters <- dat_lmtp %>% # letters for time varying
  select(starts_with("L_"),
         starts_with("A_"),
         -starts_with("L_NA_"),
         -ends_with(paste0("_",outcome_day))) %>%
  names() 
tv <- map(0:(outcome_day - 1), function(x) { # list for time varying covars
  used_letters[str_detect(used_letters, str_pad(x, 2, pad="0"))]
})


mtp <- function(data, trt) {
  tau <- readr::parse_number(trt) # extract time point
  trt_prev <- paste0("I_", stringr::str_pad(tau - 1, 2, "left", "0")) # get the col name of previous trt
  
  if(trt == "I_00") {
    data[data[[trt]] == 2, trt] <- 1
  }
  
  else {
    data[which(data[[trt]] == 2 & data[[trt_prev]] != 2), trt] <- 1
  }
  
  return(data[[trt]])
}

### With mpt
progressr::with_progress(
  out_mtp <-
    lmtp_sdr(
      dat_lmtp,
      trt = a,
      outcome = y,
      # comp_risk = cr,
      baseline = bs,
      time_vary = tv,
      cens = censoring,
      shift = mtp,
      outcome_type = "survival",
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      # .trim = trim,
      k=k,
      intervention_type = "mtp"
    )
)
# out_mtp

# ### Without mpt
# progressr::with_progress(
#   out_NULL <-
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
#       learners_trt = lrnrs,
#       folds = folds,
#       .SL_folds = SL_folds,
#       # .trim = trim,
#       k=k,
#       intervention_type = "mtp"
#     )
# )
# out_NULL
# 
