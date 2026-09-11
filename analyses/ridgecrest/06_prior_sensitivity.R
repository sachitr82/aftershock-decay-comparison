#===============================================================================
# Ridgecrest prior sensitivity (narrower rho, wider t_a,B)
#===============================================================================

library(dplyr)
library(ETAS.inlabru)
library(here)
library(inlabru)
source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Alternative prior transformations
#-------------------------------------------------------------------------------

make_links_PS <- function(kernel) {
  links <- make_links_P0(kernel)
  
  if (kernel == "mse") {
    links$rho <- \(x) loggaus_t(x, log(1.5), 0.7)
  }
  
  if (kernel == "rate_state") {
    links$B <- \(x) logitgaus_t(x, 7.5, sqrt(2))
    links$ta <- \(x) loggaus_t(x, log(200), 0.5 * sqrt(2))
  }
  
  links
}

#-------------------------------------------------------------------------------
# Preserve baseline physical initial values under alternative priors
#-------------------------------------------------------------------------------

make_bru_PS <- function(kernel) {
  opt <- make_bru_options_P0(
    kernel, rel_tol = fit_control$rel_tol, max_iter = fit_control$max_iter)
  
  if (kernel == "mse") {
    opt$bru_initial$th.rho <- inv_loggaus_t(initials$mse["rho"], log(1.5), 0.7)
  }
  
  if (kernel == "rate_state") {
    opt$bru_initial$th.B <- inv_logitgaus_t(initials$rate_state["B"], 7.5, sqrt(2))
    opt$bru_initial$th.ta <- inv_loggaus_t(initials$rate_state["ta"], log(200),
                                           0.5 * sqrt(2))
  }
  
  opt
}
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
# Log marginal likelihood
#-------------------------------------------------------------------------------

extract_lml <- function(fit) {
  stopifnot(!is.null(fit$mlik))
  
  i <- grep("integration", rownames(fit$mlik), ignore.case = TRUE)
  stopifnot(length(i) == 1)
  
  as.numeric(fit$mlik[i, 1])
}
#-------------------------------------------------------------------------------
# Baseline Ridgecrest catalogue and observation window
#-------------------------------------------------------------------------------

cat_rc <- readRDS(here("data", "derived", "ridgecrest_temporal_catalogue.rds"))

mainshock_datetime <- cat_rc$datetime[which.max(cat_rc$mag)]

cat_rc <- cat_rc %>%
  mutate(ts = as.numeric(difftime(datetime, mainshock_datetime, units = "days")),
         magnitudes = mag) %>%
  arrange(ts)

cat_rc$idx.p <- seq_len(nrow(cat_rc))

M0 <- 2.5

T_fit_start <- as.numeric(difftime(as.POSIXct("2016-01-01 00:00:00", tz = "UTC"),
                                   mainshock_datetime, units = "days"))

T_fit_end <- as.numeric(difftime(as.POSIXct("2025-12-31 23:59:59", tz = "UTC"),
                                 mainshock_datetime, units = "days"))

stopifnot(temporal_binning$N.max == 14, fit_control$rel_tol == 0.1,
          fit_control$max_iter == 100)

#-------------------------------------------------------------------------------
# Prior-sensitivity fits
#-------------------------------------------------------------------------------

out_dir <- here("outputs", "ridgecrest", "sensitivity", "prior")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
jobs <- data.frame(scenario = c("mse_narrow_rho", "rs_broad"), 
                   kernel = c("mse", "rate_state"), seed = c(510001, 510002),
                   stringsAsFactors = FALSE)
results <- vector("list", nrow(jobs))

for (i in seq_len(nrow(jobs))) {
  scenario_i <- jobs$scenario[i]
  kernel_i <- jobs$kernel[i]
  outfile_i <- file.path(out_dir, paste0("fit_", scenario_i, ".rds"))
  
  if (file.exists(outfile_i)) {
    obj_i <- readRDS(outfile_i)
  
  } else {
    set.seed(jobs$seed[i])
    links_i <- make_links_PS(kernel_i)
    bru_i <- make_bru_PS(kernel_i)
    bru_i$control.compute <- list(config = TRUE, mlik = TRUE)
    
    start_i <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = cat_rc, 
                                         M0 = M0, 
                                         T1 = T_fit_start,
                                         T2 = T_fit_end,
                                         link.functions = links_i,
                                         coef.t. = temporal_binning$coef.t, 
                                         delta.t. = temporal_binning$delta.t, 
                                         N.max. = temporal_binning$N.max,
                                         bru.opt = bru_i, kernel = kernel_i)
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
    
    obj_i <- list(fit = fit_i, 
                  link.functions = links_i,
                  scenario = scenario_i,
                  kernel = kernel_i,
                  seed = jobs$seed[i],
                  runtime_minutes = runtime_i,
                  lml = extract_lml(fit_i))
    saveRDS(obj_i, outfile_i)
  }
  
  d_i <- fit_diagnostics(obj_i$fit)
  
  results[[i]] <- data.frame(scenario = scenario_i, 
                             kernel = kernel_i, 
                             lml = obj_i$lml, 
                             runtime_minutes = obj_i$runtime_minutes,
                             d_i)
}

#-------------------------------------------------------------------------------
# Fit summary
#-------------------------------------------------------------------------------

fit_summary <- bind_rows(results)
write.csv(fit_summary, file.path(out_dir, "fit_summary.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Compare with baseline model evidence
#-------------------------------------------------------------------------------

base_dir <- here("outputs", "ridgecrest", "baseline")
lml_ou <- readRDS(file.path(base_dir, "fit_ou.rds"))$lml
lml_mse <- readRDS(file.path(base_dir, "fit_mse.rds"))$lml
lml_rs <- readRDS(file.path(base_dir, "fit_rate_state.rds"))$lml
lml_mse_ps <- fit_summary$lml[fit_summary$scenario == "mse_narrow_rho"]
lml_rs_ps <- fit_summary$lml[fit_summary$scenario == "rs_broad"]
comparison <- data.frame(comparison = c("OU - MSE", "OU - RS"), 
                         baseline_delta_lml = c(lml_ou - lml_mse, lml_ou - lml_rs), 
                         sensitivity_delta_lml = c(lml_ou - lml_mse_ps, lml_ou - lml_rs_ps))
comparison$change <- comparison$sensitivity_delta_lml - comparison$baseline_delta_lml
write.csv(comparison, file.path(out_dir, "prior_sensitivity_lml.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Results
#-------------------------------------------------------------------------------

fit_summary
comparison