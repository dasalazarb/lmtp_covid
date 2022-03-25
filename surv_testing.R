library(dplyr)
# library(lmtp)
devtools::load_all()

sim_point_surv2 <-
  sim_point_surv %>%
  mutate(CR.1 = ifelse(Y.1 == 0, rbinom(2, 1, .5), 0),
         CR.2 = ifelse(Y.2 == 0, rbinom(2, 1, .5), 0),
         CR.3 = ifelse(Y.3 == 0, rbinom(2, 1, .5), 0),
         CR.4 = ifelse(Y.4 == 0, rbinom(2, 1, .5), 0),
         CR.5 = ifelse(Y.5 == 0, rbinom(2, 1, .5), 0),
         CR.6 = ifelse(Y.6 == 0, rbinom(2, 1, .5), 0),
         )

# Example 5.1
# Time-to-event analysis with a binary time-invariant exposure. Interested in
# the effect of treatment being given to all observations on the cumulative
# incidence of the outcome.
# For a survival problem, the outcome argument now takes a vector of outcomes
# if an observation experiences the event prior to the end of follow-up, all future
# outcome nodes should be set to 1 (i.e., last observation carried forward).
A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
CR <- paste0("CR.", 1:6)
W <- c("W1", "W2")

debug(lmtp_sdr)
debug(cf_sdr)
debug(at_risk)
debug(check_at_risk)

lmtp_sdr(
  sim_point_surv2, A, Y,
  comp_risk = CR, baseline = W, cens = C, folds = 2,
  shift = static_binary_on, outcome_type = "survival"
)


lmtp_sdr(
  sim_point_surv, A, Y,
  baseline = W, cens = C, folds = 2,
  shift = static_binary_on, outcome_type = "survival"
)
