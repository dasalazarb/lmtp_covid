  options(java.parameters = "-Xmx20000m")

# remotes::install_local(here::here("lmtp_1.0.0.5001.tar.gz"))

library(sl3)
library(lmtp)
library(tidyverse)
library(future)
library(earth)

progressr::handlers(global = TRUE)
set.seed(7)

task <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

plan(multicore)

dat_lmtp <- read_rds(here::here("data/dat_final.rds")) %>%
  filter(hypoxia_ed == 1)

trim <- .995
folds <- 10
SL_folds <- 10
k <- 2

results_folder <- Sys.Date()
trim_num <- str_split(as.character(trim),"\\.")[[1]][2]

folder_to_save <- paste0("/home/kah2797/comprisks/results/", results_folder)
dir.create(file.path(folder_to_save), showWarnings = FALSE)
file_to_save <- paste0(task, "_tv_locf_", trim_num, "_k", k, "_f", folds, "_notv.rds")

mars_grid_params <- list(
  degree = c(2,3),
  penalty = c(1,2,3)
)

mars_grid <- expand.grid(mars_grid_params, KEEP.OUT.ATTRS = FALSE)
mars_learners <- apply(mars_grid, MARGIN = 1, function(tuning_params) {
  do.call(Lrnr_earth$new, as.list(tuning_params))
})


lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_enet <- Lrnr_glmnet$new(alpha = 0.5)
lrn_bart <- Lrnr_bartMachine$new()
lrn_mean <- Lrnr_mean$new()


learners <- unlist(list(
  mars_learners,
  lrn_lasso,
  lrn_ridge,
  lrn_enet,
  lrn_bart,
  lrn_mean
),
recursive = TRUE
)


lrnrs <- make_learner(Stack, learners)

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
cr <- paste0("CR_",padded_days_out) # competing risk (death)
censoring <- paste0("C_",padded_days) # observed at next time

used_letters <- dat_lmtp %>% # letters for time varying
  select(starts_with("L_"),
         starts_with("A_"),
         starts_with("H_"),
         -starts_with("L_NA_"),
         -ends_with(paste0("_",outcome_day))) %>%
  names() 
tv <- map(0:(outcome_day - 1), function(x) { # list for time varying covars
  used_letters[str_detect(used_letters, str_pad(x, 2, pad="0"))]
})


mtp <- function(data, trt) {
  tau <- readr::parse_number(trt) # extract time point
  trt_prev <- paste0("I_", stringr::str_pad(tau - 1, 2, "left", "0")) # get the col name of previous trt
  
  trt_prev

  if(trt == "I_00"){
    data[[trt]]
  }

  else if(data[["H_00"]] == 0 && data[[trt]] == 1 && data[[trt_prev]] == 0) {
    data[[trt]] <- 0
    data[[trt]]
  }

  else{
    data[[trt]]
  }

}

if (task == 1) {
progressr::with_progress(
  out <-
    lmtp_sdr(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = mtp,
      outcome_type = "survival",
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = 5,
      .SL_folds = 5,
      .trim = trim,
      k=k,
      intervention_type = "mtp"
    )
)
}

if (task == 2) {
  out <-
    lmtp_sdr(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = NULL,
      outcome_type = "survival",
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = 5,
      .SL_folds = 5,
      .trim = trim,
      k=k,
      intervention_type = "static"
    )
}

if (task == 3) {
  
  out <-
    lmtp_tmle(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = mtp,
      outcome_type = "survival",
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      .trim = trim,
      k=k,
      intervention_type = "dynamic"
    )
}

if (task == 4) {
  out <-
    lmtp_tmle(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = NULL,
      outcome_type = "survival",
      learners_outcome = lrnrs,
      learners_trt = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      .trim = trim,
      k=k,
      intervention_type = "static"
    )
}

if (task == 5) {
  
  out <-
    lmtp_sub(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = mtp,
      outcome_type = "survival",
      learners = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      .trim = trim,
      k=k,
      intervention_type = "dynamic"
    )
}

if (task == 6) {
  out <-
    lmtp_sub(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = NULL,
      outcome_type = "survival",
      learners = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      .trim = trim,
      k=k,
      intervention_type = "static"
    )
}
  



if (task == 7) {
  out <-
    lmtp_ipw(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = mtp,
      outcome_type = "survival",
      learners = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      k=k
    )
}


if (task == 8) {
  out <-
    lmtp_ipw(
      dat_lmtp,
      trt = a,
      outcome = y,
      comp_risk = cr,
      baseline = bs,
      cens = censoring,
      shift = NULL,
      outcome_type = "survival",
      learners = lrnrs,
      folds = folds,
      .SL_folds = SL_folds,
      k=k
    )
}


# save in new cohort folder
saveRDS(out, paste0(folder_to_save, "/", file_to_save))


