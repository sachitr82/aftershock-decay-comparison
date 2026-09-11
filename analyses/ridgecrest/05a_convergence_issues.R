#===============================================================================
# Diagnose unconverged magnitude-threshold fits
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and M0 sensitivity fits
#-------------------------------------------------------------------------------

library(dplyr)
library(here)

fit_dir <- here(
  "outputs", "ridgecrest", "sensitivity", "magnitude_threshold", "fits")

jobs <- data.frame(M0 = c(3.0, 3.0, 3.3, 3.3, 3.3),
                   kernel = c("mse", "rate_state", "ou", "mse", "rate_state"),
                  file = c("fit_M0_3_mse.rds",
                           "fit_M0_3_rate_state.rds",
                           "fit_M0_3p3_ou.rds",
                           "fit_M0_3p3_mse.rds",
                           "fit_M0_3p3_rate_state.rds"),
                  stringsAsFactors = FALSE)

fits <- lapply(file.path(fit_dir, jobs$file), readRDS)

#-------------------------------------------------------------------------------
# Formal convergence status
#-------------------------------------------------------------------------------

status <- bind_rows(lapply(seq_along(fits), function(i) {
  x <- fits[[i]]
  
  data.frame(M0 = jobs$M0[i],
                  kernel = jobs$kernel[i],
                  fit_status = x$fit_status,
                  converged = x$converged,
                  hit_max = x$hit_max,
                  n_iter = x$n_iter,
                  lml = x$lml)
}))

status

#-------------------------------------------------------------------------------
# Last 10 outer-iteration log-posterior modes
#-------------------------------------------------------------------------------

last10_logpost <- bind_rows(lapply(seq_along(fits), function(i) {
  track <- fits[[i]]$fit$bru_iinla$track
  z <- track[track$effect == "log.posterior.mode", c("iteration", "mode")]
  z <- tail(z, 10)
  
  data.frame(M0 = jobs$M0[i],
             kernel = jobs$kernel[i],
             iteration = z$iteration,
             log_posterior_mode = z$mode)
}))

last10_logpost

#-------------------------------------------------------------------------------
# Check for period-two cycling
#-------------------------------------------------------------------------------

cycle_summary <- bind_rows(lapply(seq_along(fits), function(i) {
  track <- fits[[i]]$fit$bru_iinla$track
  z <- track[track$effect == "log.posterior.mode", "mode"]
  z <- tail(z, 10)
  
  data.frame(M0 = jobs$M0[i],
             kernel = jobs$kernel[i],
             n_unique = length(unique(round(z, 6))),
             range = diff(range(z)),
             max_lag1_change = max(abs(diff(z))),
             max_lag2_change = max(abs(z[3:length(z)] - z[1:(length(z) - 2)])))
}))

cycle_summary

#-------------------------------------------------------------------------------
# Last 10 parameter modes
#-------------------------------------------------------------------------------

last10_parameters <- bind_rows(lapply(seq_along(fits), function(i) {
  track <- fits[[i]]$fit$bru_iinla$track
  last_iter <- max(track$iteration)
  
  z <- track[track$effect != "log.posterior.mode" & 
               track$iteration > last_iter - 10,
             c("effect", "iteration", "mode", "sd", "new_linearisation")]
  
  data.frame(M0 = jobs$M0[i], kernel = jobs$kernel[i], z)
}))

last10_parameters

#-------------------------------------------------------------------------------
# Terminal outer-convergence messages
#-------------------------------------------------------------------------------

terminal_logs <- lapply(seq_along(fits), function(i) {
  log_i <- as.character(inlabru::bru_log(fits[[i]]$fit))
  
  keep <- grepl(
    "Max deviation from previous|Convergence criterion met|Maximum iterations",
    log_i)
  
  data.frame(M0 = jobs$M0[i], kernel = jobs$kernel[i], 
             message = tail(log_i[keep], 6))
})

terminal_logs <- bind_rows(terminal_logs)
terminal_logs