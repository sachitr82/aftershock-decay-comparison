#===============================================================================
# Ridgecrest prior sensitivity (narrower rho; wider t_a, B)
#===============================================================================

library(dplyr)
library(ETAS.inlabru)
library(here)
library(inlabru)
source(here("analyses", "ridgecrest", "00_design.R"))
source(here("src", "fit_helpers", "fit_diagnostics.R"))

#-------------------------------------------------------------------------------
# Alternative prior transformations
#-------------------------------------------------------------------------------

make_links_PS <- function(form) {
  links <- make_links_P0(form)
  
  if (form == "mse") {
    links$rho <- \(x) loggaus_t(x, log(1.5), 0.7)
  }
  
  if (form == "rate_state") {
    links$B <- \(x) logitgaus_t(x, 7.5, sqrt(2))
    links$ta <- \(x) loggaus_t(x, log(200), 0.5 * sqrt(2))
  }
  
  links
}

#-------------------------------------------------------------------------------
# Preserve baseline physical initial values under alternative priors
#-------------------------------------------------------------------------------

make_bru_PS <- function(form) {
  opt <- make_bru_options_P0(
    form, rel_tol = fit_control$rel_tol, max_iter = fit_control$max_iter)
  
  if (form == "mse") {
    opt$bru_initial$th.rho <- inv_loggaus_t(initials$mse["rho"], log(1.5), 0.7)
  }
  
  if (form == "rate_state") {
    opt$bru_initial$th.B <- inv_logitgaus_t(initials$rate_state["B"], 7.5, sqrt(2))
    opt$bru_initial$th.ta <- inv_loggaus_t(initials$rate_state["ta"], log(200),
                                           0.5 * sqrt(2))
  }
  
  opt
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

# Mainshock defines t = 0
mainshock_datetime <- cat_rc$datetime[which.max(cat_rc$mag)]

# Convert to form expected by Temporal.ETAS()
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

# Verify numerical design
stopifnot(temporal_binning$N.max == 14, fit_control$rel_tol == 0.1,
          fit_control$max_iter == 100)

#-------------------------------------------------------------------------------
# Prior-sensitivity fits
#-------------------------------------------------------------------------------

jobs <- data.frame(scenario = c("mse_narrow_rho", "rs_broad"), 
                   form = c("mse", "rate_state"), seed = c(510001, 510002),
                   stringsAsFactors = FALSE)
results <- vector("list", nrow(jobs))

for (i in seq_len(nrow(jobs))) {
  scenario_i <- jobs$scenario[i]
  form_i <- jobs$form[i]
  outfile_i <- file.path(prior_fit_dir, paste0("fit_", scenario_i, ".rds"))
  
  if (file.exists(outfile_i)) {
    obj_i <- readRDS(outfile_i)
  
  } else {
    set.seed(jobs$seed[i])
    
    # Two additional fits based on the two additional priors
    links_i <- make_links_PS(form_i)
    bru_i <- make_bru_PS(form_i)
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
                                         bru.opt = bru_i, form = form_i)
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
    
    diag_i <- fit_diagnostics(fit_i)
    usable_i <- fit_is_usable(diag_i)
    lml_i <- extract_lml(fit_i)
    
    obj_i <- list(fit = fit_i,
                  link.functions = links_i,
                  scenario = scenario_i,
                  fitted_form = form_i,
                  fit_seed = jobs$seed[i],
                  runtime_minutes = runtime_i,
                  diagnostics = diag_i,
                  usable = usable_i,
                  lml = lml_i)
    
    # Save atomically via a temporary file to avoid incomplete fit objects
    tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i),
                      fileext = ".rds")
    
    saveRDS(obj_i, tmp_i)
    
    if (!file.rename(tmp_i, outfile_i)) {
      stop("Could not move temporary fit file to final output.")
    }
  }
  
  d_i <- obj_i$diagnostics
  
  results[[i]] <- data.frame(scenario = scenario_i,
                  fitted_form = form_i,
                  lml = obj_i$lml,
                  runtime_minutes = obj_i$runtime_minutes,
                  fit_status = d_i$fit_status,
                  converged = d_i$converged,
                  hit_max = d_i$hit_max,
                  inla_failure = d_i$inla_failure,
                  nan_inf_logl = d_i$nan_inf_logl,
                  vb_aborted = d_i$vb_aborted,
                  n_iter = d_i$n_iter,
                  min_sd = d_i$min_sd,
                  degenerate = d_i$degenerate,
                  usable = obj_i$usable)
}

#-------------------------------------------------------------------------------
# Fit summary
#-------------------------------------------------------------------------------

fit_summary <- bind_rows(results)
write.csv(fit_summary, file.path(prior_table_dir, "fit_summary.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Compare with baseline model evidence
#-------------------------------------------------------------------------------

stopifnot(all(fit_summary$usable))
stopifnot(all(is.finite(fit_summary$lml)))

lml_ou <- readRDS(file.path(baseline_fit_dir, "fit_ou.rds"))$lml
lml_mse <- readRDS(file.path(baseline_fit_dir, "fit_mse.rds"))$lml
lml_rs <- readRDS(file.path(baseline_fit_dir, "fit_rate_state.rds"))$lml

lml_mse_ps <- fit_summary$lml[fit_summary$scenario == "mse_narrow_rho"]
lml_rs_ps <- fit_summary$lml[fit_summary$scenario == "rs_broad"]
comparison <- data.frame(comparison = c("OU - MSE", "OU - RS"), 
                         baseline_delta_lml = c(lml_ou - lml_mse, lml_ou - lml_rs), 
                         sensitivity_delta_lml = c(lml_ou - lml_mse_ps, lml_ou - lml_rs_ps))
comparison$change <- comparison$sensitivity_delta_lml - comparison$baseline_delta_lml
write.csv(comparison, file.path(prior_table_dir, "prior_sensitivity_lml.csv"), row.names = FALSE)

message("Saved prior-sensitivity fits to: ", prior_fit_dir)
message("Saved prior-sensitivity tables to: ", prior_table_dir)
message("Finished Ridgecrest prior sensitivity.")