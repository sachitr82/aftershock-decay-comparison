#===============================================================================
# Fit Ridgecrest catalogue
#===============================================================================

library(dplyr)
library(ETAS.inlabru)
library(here)
library(inlabru)

source(here("analyses", "simulation", "00_design.R"))

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

candidate_kernels <- c("ou", "mse", "rate_state")

stopifnot(temporal_binning$N.max == 14)
stopifnot(fit_control$rel_tol == 0.1)
stopifnot(fit_control$max_iter == 100)

ridgecrest_fit_dir <- here("outputs", "ridgecrest", "baseline")
dir.create(ridgecrest_fit_dir, recursive = TRUE, showWarnings = FALSE)

fit_seeds <- c(ou = 400001, mse = 400002, rate_state = 400003)

#-------------------------------------------------------------------------------
# Fit diagnostics
#-------------------------------------------------------------------------------

fit_diagnostics <- function(fit) {
  log_i <- as.character(inlabru::bru_log(fit))
  
  inla_failure <- any(grepl(
    "Problem in inla|Giving up and returning last successfully obtained result|inla-program exited with an error|maximum number of tries has been reached|Newton-Raphson optimizer did not converge",
    log_i))
  
  nan_inf_logl <- any(grepl("NAN/INF values in logl", log_i, fixed = TRUE))
  vb_aborted <- any(grepl("max_correction|vb.correction.*aborted", log_i))
  
  converged <- !inla_failure &&
    any(grepl("Convergence criterion met", log_i, fixed = TRUE))
  
  hit_max <- !inla_failure &&
    any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
  
  fit_status <- case_when(inla_failure ~ "inla_failure",
                          converged ~ "converged",
                          hit_max ~ "hit_max",
                          TRUE ~ "unknown")
  
  n_iter <- max(fit$bru_iinla$track$iteration, na.rm = TRUE)
  
  list(fit_status = fit_status,
       converged = converged,
       hit_max = hit_max,
       inla_failure = inla_failure,
       nan_inf_logl = nan_inf_logl,
       vb_aborted = vb_aborted,
       n_iter = n_iter)
}

#-------------------------------------------------------------------------------
# Posterior diagnostics
#-------------------------------------------------------------------------------

posterior_diagnostics <- function(fit) {
  x <- fit$summary.fixed
  
  data.frame(latent_parameter = rownames(x),
             mean = x$mean,
             sd = x$sd,
             q025 = x$`0.025quant`,
             median = x$`0.5quant`,
             q975 = x$`0.975quant`,
             degenerate = !is.finite(x$sd) | x$sd < 1e-6,
             row.names = NULL)
}

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

fit_manifest <- vector("list", length(candidate_kernels))

for (k in seq_along(candidate_kernels)) {
  
  fitted_i <- candidate_kernels[k]
  
  outfile_i <- file.path(ridgecrest_fit_dir, paste0("fit_", fitted_i, ".rds"))
  
  message("============================================================")
  message("Ridgecrest | Fit: ", fitted_i)
  message("============================================================")
  
  if (file.exists(outfile_i)) {
    
    message("Already exists: ", basename(outfile_i))
    output_i <- readRDS(outfile_i)
  } else {
    
    set.seed(fit_seeds[[fitted_i]])
    
    link_i <- make_links_P0(fitted_i)
    
    bru_i <- make_bru_options_P0(fitted_i, rel_tol = fit_control$rel_tol,
                                 max_iter = fit_control$max_iter)
    
    bru_i$control.compute <- list(config = TRUE, mlik = TRUE)
    
    start_i <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = cat_rc,
                                         M0 = M0,
                                         T1 = T_fit_start,
                                         T2 = T_fit_end,
                                         link.functions = link_i,
                                         coef.t. = temporal_binning$coef.t,
                                         delta.t. = temporal_binning$delta.t,
                                         N.max. = temporal_binning$N.max,
                                         bru.opt = bru_i,
                                         kernel = fitted_i)
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
    
    diag_i <- fit_diagnostics(fit_i)
    posterior_diag_i <- posterior_diagnostics(fit_i)
    lml_i <- extract_lml(fit_i)
    
    output_i <- list(fit = fit_i,
                     link.functions = link_i,
                     fitted_kernel = fitted_i,
                     M0 = M0,
                     T1 = T_fit_start,
                     T2 = T_fit_end,
                     n_fit = nrow(cat_rc),
                     binning_parameters = temporal_binning,
                     rel_tol = fit_control$rel_tol,
                     max_iter = fit_control$max_iter,
                     fit_seed = fit_seeds[[fitted_i]],
                     runtime_minutes = runtime_i,
                     fit_status = diag_i$fit_status,
                     converged = diag_i$converged,
                     hit_max = diag_i$hit_max,
                     inla_failure = diag_i$inla_failure,
                     nan_inf_logl = diag_i$nan_inf_logl,
                     vb_aborted = diag_i$vb_aborted,
                     n_iter = diag_i$n_iter,
                     posterior_diagnostics = posterior_diag_i,
                     lml = lml_i)
    
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
  
  posterior_diag_i <- posterior_diagnostics(output_i$fit)
  
  n_degenerate_i <- sum(posterior_diag_i$degenerate)
  posterior_valid_i <- n_degenerate_i == 0
  
  if (!posterior_valid_i) {
    warning(fitted_i, " fit contains ", 
            n_degenerate_i, " near-degenerate latent posterior marginal(s).")
  }
  
  fit_manifest[[k]] <- data.frame(fitted_kernel = fitted_i,
                                  n_fit = output_i$n_fit,
                                  lml = output_i$lml,
                                  runtime_minutes = output_i$runtime_minutes,
                                  fit_status = output_i$fit_status,
                                  converged = output_i$converged,
                                  posterior_valid = posterior_valid_i,
                                  hit_max = output_i$hit_max,
                                  inla_failure = output_i$inla_failure,
                                  nan_inf_logl = output_i$nan_inf_logl,
                                  vb_aborted = output_i$vb_aborted,
                                  n_iter = output_i$n_iter,
                                  n_degenerate = n_degenerate_i,
                                  file = basename(outfile_i))
}
#-------------------------------------------------------------------------------
# Save Ridgecrest fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- bind_rows(fit_manifest)

fit_manifest <- fit_manifest %>% arrange(desc(lml))

print(fit_manifest)

write.csv(fit_manifest, file.path(ridgecrest_fit_dir, "fit_manifest.csv"),
          row.names = FALSE)

message("Finished Ridgecrest baseline fits.")