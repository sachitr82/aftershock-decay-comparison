#===============================================================================
# rel_tol convergence check
# NOTE: On an Apple M4 Mac, the full run takes approximately 15 minutes.
#===============================================================================

#-------------------------------------------------------------------------------
# Load package and design choices
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(future)
library(here)
library(inlabru)

source(here("analyses", "simulation", "00_design.R"))
source(here("src", "fit_helpers", "fit_diagnostics.R"))
source(here("src", "fit_helpers", "fit_temporal_etas.R"))

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Load fixed pilot catalogues (from 01)
#-------------------------------------------------------------------------------

pilot_files <- c(ou = file.path(pilot_catalogue_dir, "ou_pilot.rds"),
                 mse = file.path(pilot_catalogue_dir, "mse_pilot.rds"),
                 rate_state = file.path(pilot_catalogue_dir, "rate_state_pilot.rds"))

pilots <- lapply(pilot_files, readRDS)

#-------------------------------------------------------------------------------
# Candidate fitted decay forms
#-------------------------------------------------------------------------------

candidate_forms <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Targeted rel_tol refinement with temporal binning held fixed
#-------------------------------------------------------------------------------

rel_tol_schemes <- c(RT01 = 0.1, RT005 = 0.05, RT001 = 0.01)

#-------------------------------------------------------------------------------
# Fit index: one correctly specified fit per truth under each rel_tol (9 fits)
#-------------------------------------------------------------------------------

rel_tol_index <- do.call(rbind, lapply(candidate_forms, function(form) {
  data.frame(truth_form = form,
             fitted_form = form,
             tolerance = names(rel_tol_schemes),
             stringsAsFactors = FALSE)
}))

rownames(rel_tol_index) <- NULL

stopifnot(nrow(rel_tol_index) == 9,
          all(rel_tol_index$truth_form == rel_tol_index$fitted_form))

#-------------------------------------------------------------------------------
# Storage for rel_tol-fit manifest
#-------------------------------------------------------------------------------

rel_tol_manifest <- vector("list", nrow(rel_tol_index))

#-------------------------------------------------------------------------------
# Run rel_tol fits
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(rel_tol_index))) {
  
  # Identify generating form, fitted form and convergence tolerance
  truth_i <- rel_tol_index$truth_form[i]
  fitted_i <- rel_tol_index$fitted_form[i]
  tolerance_i <- rel_tol_index$tolerance[i]
  rel_tol_i <- rel_tol_schemes[[tolerance_i]]
  
  # Prepare the fixed pilot catalogue for fitting
  catalogue_i <- prepare_temporal_catalogue(pilots[[truth_i]]$catalogue,
                                            T1 = T_fit_start,
                                            T2 = T_fit_end)
  
  # Use common fitting controls with the candidate convergence tolerance
  fit_control_i <- fit_control
  fit_control_i$rel_tol <- rel_tol_i
  
  outfile_i <- file.path(
    rel_tol_fit_dir,
    paste0("truth_", truth_i, "_fit_", fitted_i, "_", tolerance_i, ".rds")
  )
  
  message("Truth: ", truth_i, " | Fit: ", fitted_i, " | rel_tol: ", rel_tol_i)
  
  if (file.exists(outfile_i)) {
    
    # Reuse completed fit if present
    message("Already exists: ", basename(outfile_i), " — skipping fit")
    
    obj_i <- readRDS(outfile_i)
    
    diag_i <- fit_diagnostics(obj_i$fit)
    usable_i <- fit_is_usable(diag_i)
    runtime_i <- obj_i$runtime_minutes
    
  } else {
    
    # Fit correctly specified ETAS model under the candidate tolerance
    result_i <- fit_temporal_etas(catalogue = catalogue_i,
                                  fitted_form = fitted_i,
                                  binning = temporal_binning,
                                  fit_control = fit_control_i,
                                  M0 = M0,
                                  T1 = T_fit_start,
                                  T2 = T_fit_end,
                                  compute_model_criteria = FALSE)
    
    fit_i <- result_i$fit
    link_i <- result_i$link.functions
    runtime_i <- result_i$runtime_minutes
    diag_i <- result_i$diagnostics
    usable_i <- result_i$usable
    
    # Save fitted object and diagnostics
    saveRDS(list(fit = fit_i,
                 link.functions = link_i,
                 truth_form = truth_i,
                 fitted_form = fitted_i,
                 tolerance = tolerance_i,
                 rel_tol = rel_tol_i,
                 binning_parameters = temporal_binning,
                 runtime_minutes = runtime_i,
                 diagnostics = diag_i,
                 usable = usable_i), outfile_i)
    
    message("Saved ", basename(outfile_i),
            " | runtime = ", round(runtime_i, 2), " min")
  }
  
  # Record one-row summary for the rel_tol manifest
  rel_tol_manifest[[i]] <- data.frame(truth_form = truth_i,
                                      fitted_form = fitted_i,
                                      tolerance = tolerance_i,
                                      rel_tol = rel_tol_i,
                                      coef_t = temporal_binning$coef.t,
                                      delta_t = temporal_binning$delta.t,
                                      N_max = temporal_binning$N.max,
                                      n_fit = nrow(catalogue_i),
                                      runtime_minutes = runtime_i,
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
# Save fit manifest
#-------------------------------------------------------------------------------

rel_tol_manifest <- do.call(rbind, rel_tol_manifest)
rownames(rel_tol_manifest) <- NULL

manifest_file <- file.path(rel_tol_fit_dir, "rel_tol_manifest.csv")

write.csv(rel_tol_manifest, manifest_file, row.names = FALSE)

stopifnot(nrow(rel_tol_manifest) == 9)

message("Saved ", manifest_file)
message("Finished all 9 rel_tol-check fits.")