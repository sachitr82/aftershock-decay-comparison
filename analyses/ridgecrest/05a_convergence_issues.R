#===============================================================================
# Diagnose unconverged magnitude-threshold fits
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and M0 sensitivity fits
#-------------------------------------------------------------------------------

library(dplyr)
library(here)

source(here("analyses", "ridgecrest", "00_design.R"))

jobs <- data.frame(M0 = c(3.0, 3.0, 3.3, 3.3, 3.3),
                   form = c("mse", "rate_state", "ou", "mse", "rate_state"),
                  file = c("fit_M0_3_mse.rds",
                           "fit_M0_3_rate_state.rds",
                           "fit_M0_3p3_ou.rds",
                           "fit_M0_3p3_mse.rds",
                           "fit_M0_3p3_rate_state.rds"),
                  stringsAsFactors = FALSE)

fits <- lapply(file.path(magnitude_fit_dir, jobs$file), readRDS)

#-------------------------------------------------------------------------------
# Formal convergence status
#-------------------------------------------------------------------------------

status <- bind_rows(lapply(seq_along(fits), function(i) {
  x <- fits[[i]]
  d <- x$diagnostics
  
  data.frame(M0 = jobs$M0[i],
             form = jobs$form[i],
             fit_status = d$fit_status,
             converged = d$converged,
             hit_max = d$hit_max,
             n_iter = d$n_iter,
             min_sd = d$min_sd,
             degenerate = d$degenerate,
             usable = x$usable,
             lml = x$lml)
}))

#-------------------------------------------------------------------------------
# Check for period-two cycling
#-------------------------------------------------------------------------------

cycle_summary <- bind_rows(lapply(seq_along(fits), function(i) {
  track <- fits[[i]]$fit$bru_iinla$track
  z <- track[track$effect == "log.posterior.mode", "mode"]
  z <- tail(z, 10)
  
  data.frame(M0 = jobs$M0[i],
             form = jobs$form[i],
             n_unique = length(unique(round(z, 6))),
             range = diff(range(z)),
             max_lag1_change = max(abs(diff(z))),
             max_lag2_change = max(abs(z[3:length(z)] - z[1:(length(z) - 2)])))
}))

#-------------------------------------------------------------------------------
# Last 10 parameter modes
#-------------------------------------------------------------------------------

last10_parameters <- bind_rows(lapply(seq_along(fits), function(i) {
  track <- fits[[i]]$fit$bru_iinla$track
  last_iter <- max(track$iteration)
  
  z <- track[track$effect != "log.posterior.mode" & 
               track$iteration > last_iter - 10,
             c("effect", "iteration", "mode", "sd", "new_linearisation")]
  
  data.frame(M0 = jobs$M0[i], form = jobs$form[i], z)
}))

#-------------------------------------------------------------------------------
# Terminal outer-convergence messages
#-------------------------------------------------------------------------------

terminal_logs <- lapply(seq_along(fits), function(i) {
  log_i <- as.character(inlabru::bru_log(fits[[i]]$fit))
  
  keep <- grepl(
    "Max deviation from previous|Convergence criterion met|Maximum iterations",
    log_i)
  
  data.frame(M0 = jobs$M0[i], form = jobs$form[i], 
             message = tail(log_i[keep], 6))
})

terminal_logs <- bind_rows(terminal_logs)

write.csv(last10_parameters,
          file.path(magnitude_table_dir,
                    "convergence_parameter_modes.csv"),
          row.names = FALSE)

write.csv(status,
          file.path(magnitude_table_dir, "convergence_status.csv"),
          row.names = FALSE)

write.csv(cycle_summary,
          file.path(magnitude_table_dir, "convergence_cycle_summary.csv"),
          row.names = FALSE)

write.csv(terminal_logs,
          file.path(magnitude_table_dir, "convergence_terminal_logs.csv"),
          row.names = FALSE)

message("Saved convergence diagnostics to: ", magnitude_table_dir)
message("Finished magnitude-threshold convergence diagnostics.")