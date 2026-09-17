#===============================================================================
# Fit Ridgecrest catalogue
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ETAS.inlabru)
library(here)
library(inlabru)

source(here("analyses", "ridgecrest", "00_design.R"))
source(here("src", "fit_helpers", "fit_diagnostics.R"))
source(here("src", "fit_helpers", "fit_temporal_etas.R"))

#-------------------------------------------------------------------------------
# Load and prepare final Ridgecrest catalogue
#-------------------------------------------------------------------------------

cat_rc <- readRDS(here("data", "derived", "ridgecrest_temporal_catalogue.rds"))

# Mainshock defines t = 0
mainshock_i <- which.max(cat_rc$mag)
mainshock_datetime <- cat_rc$datetime[mainshock_i]

# Convert to form expected by Temporal.ETAS()
cat_rc <- cat_rc %>%
  mutate(ts = as.numeric(difftime(datetime, mainshock_datetime, units = "days")),
         magnitudes = mag) %>%
  arrange(ts)

cat_rc$idx.p <- seq_len(nrow(cat_rc))

#-------------------------------------------------------------------------------
# Ridgecrest observation window
#-------------------------------------------------------------------------------

M0 <- 2.5

window_start <- as.POSIXct("2016-01-01 00:00:00", tz = "UTC")
window_end <- as.POSIXct("2025-12-31 23:59:59", tz = "UTC")


T_fit_start <- as.numeric(
  difftime(window_start, mainshock_datetime, units = "days"))

T_fit_end <- as.numeric(
  difftime(window_end, mainshock_datetime, units = "days"))

#-------------------------------------------------------------------------------
# Catalogue checks
#-------------------------------------------------------------------------------

stopifnot(is.data.frame(cat_rc))
stopifnot(all(c("ts", "magnitudes", "idx.p") %in% names(cat_rc)))
stopifnot(all(is.finite(cat_rc$ts)))
stopifnot(all(is.finite(cat_rc$magnitudes)))
stopifnot(all(cat_rc$magnitudes >= M0))
stopifnot(all(diff(cat_rc$ts) >= 0))
stopifnot(any(cat_rc$ts == 0))
stopifnot(min(cat_rc$ts) >= T_fit_start)
stopifnot(max(cat_rc$ts) <= T_fit_end)

cat("Ridgecrest catalogue:", nrow(cat_rc), "events\n")
cat("Observation window:", T_fit_start, "to", T_fit_end, "days\n")
cat("Event-time range:", range(cat_rc$ts), "days\n")
cat("Magnitude range:", range(cat_rc$magnitudes), "\n")

#-------------------------------------------------------------------------------
# Fixed fitting design
#-------------------------------------------------------------------------------

candidate_forms <- c("ou", "mse", "rate_state")

stopifnot(temporal_binning$N.max == 14)
stopifnot(fit_control$rel_tol == 0.1)
stopifnot(fit_control$max_iter == 100)

fit_seeds <- c(ou = 400001, mse = 400002, rate_state = 400003)

#-------------------------------------------------------------------------------
# Log marginal likelihood
#-------------------------------------------------------------------------------

extract_lml <- function(fit) {
  stopifnot(!is.null(fit$mlik))
  
  i <- grep("integration", rownames(fit$mlik), ignore.case = TRUE)
  stopifnot(length(i) == 1)
  
  as.numeric(fit$mlik[i, 1])
}

#-------------------------------------------------------------------------------
# Fit OU, MSE and rate-state
#-------------------------------------------------------------------------------

fit_manifest <- vector("list", length(candidate_forms))

for (k in seq_along(candidate_forms)) {
  
  fitted_i <- candidate_forms[k]
  
  outfile_i <- file.path(baseline_fit_dir, paste0("fit_", fitted_i, ".rds"))
  
  message("============================================================")
  message("Ridgecrest | Fit: ", fitted_i)
  message("============================================================")
  
  if (file.exists(outfile_i)) {
    
    message("Already exists: ", basename(outfile_i))
    output_i <- readRDS(outfile_i)
  } else {
    
    # fit each candidate form
    result_i <- fit_temporal_etas(catalogue = cat_rc,
                                  fitted_form = fitted_i,
                                  binning = temporal_binning,
                                  fit_control = fit_control,
                                  M0 = M0,
                                  T1 = T_fit_start,
                                  T2 = T_fit_end,
                                  compute_model_criteria = TRUE,
                                  seed = fit_seeds[[fitted_i]])
    
    fit_i <- result_i$fit
    link_i <- result_i$link.functions
    runtime_i <- result_i$runtime_minutes
    diag_i <- result_i$diagnostics
    usable_i <- result_i$usable
    lml_i <- extract_lml(fit_i)
    
    output_i <- list(fit = fit_i,
                     link.functions = link_i,
                     fitted_form = fitted_i,
                     M0 = M0,
                     T1 = T_fit_start,
                     T2 = T_fit_end,
                     n_fit = nrow(cat_rc),
                     binning_parameters = temporal_binning,
                     rel_tol = fit_control$rel_tol,
                     max_iter = fit_control$max_iter,
                     fit_seed = fit_seeds[[fitted_i]],
                     runtime_minutes = runtime_i,
                     diagnostics = diag_i,
                     usable = usable_i,
                     lml = lml_i)
    
    # Save atomically via a temporary file to avoid incomplete fit objects
    tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i), 
                      fileext = ".rds")
    
    saveRDS(output_i, tmp_i)
    
    if (!file.rename(tmp_i, outfile_i)) {
      stop("Could not move temporary fit file to final output.")
    }
    
    message("Saved ", basename(outfile_i), 
            " | runtime = ", round(runtime_i, 2), " min", 
            " | LML = ", round(lml_i, 2))
  }
  
  diag_i <- output_i$diagnostics
  usable_i <- output_i$usable
  
  fit_manifest[[k]] <- data.frame(fitted_form = fitted_i,
                                  n_fit = output_i$n_fit,
                                  lml = output_i$lml,
                                  runtime_minutes = output_i$runtime_minutes,
                                  fit_status = diag_i$fit_status,
                                  converged = diag_i$converged,
                                  hit_max = diag_i$hit_max,
                                  inla_failure = diag_i$inla_failure,
                                  nan_inf_logl = diag_i$nan_inf_logl,
                                  vb_aborted = diag_i$vb_aborted,
                                  n_iter = diag_i$n_iter,
                                  min_sd = diag_i$min_sd,
                                  degenerate = diag_i$degenerate,
                                  usable = usable_i,
                                  file = basename(outfile_i))
}
#-------------------------------------------------------------------------------
# Save Ridgecrest fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- bind_rows(fit_manifest)

fit_manifest <- fit_manifest %>% arrange(desc(lml))

write.csv(fit_manifest,
          file.path(baseline_fit_summary_dir, "fit_manifest.csv"),
          row.names = FALSE)

message("Saved baseline fits to: ", baseline_fit_dir)
message("Saved baseline fit manifest to: ", baseline_fit_summary_dir)
message("Finished Ridgecrest baseline fits.")