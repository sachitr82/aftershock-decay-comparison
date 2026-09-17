#===============================================================================
# Temporal binning check - choosing max number of bins for fitting
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
# Targeted N.max refinement with coef.t and delta.t held fixed
#-------------------------------------------------------------------------------

binning_schemes <- list(N8  = list(coef.t = 1, delta.t = 0.1, N.max = 8),
                        N10 = list(coef.t = 1, delta.t = 0.1, N.max = 10),
                        N12 = list(coef.t = 1, delta.t = 0.1, N.max = 12),
                        N14 = list(coef.t = 1, delta.t = 0.1, N.max = 14),
                        N16 = list(coef.t = 1, delta.t = 0.1, N.max = 16))

#-------------------------------------------------------------------------------
# Fit index: one correctly specified fit per truth under five N.max values (15 fits)
#-------------------------------------------------------------------------------

binning_index <- do.call(rbind, lapply(candidate_forms, function(form) {
  data.frame(truth_form = form,
             fitted_form = form, 
             binning = names(binning_schemes),
             stringsAsFactors = FALSE)
}))

rownames(binning_index) <- NULL

stopifnot(nrow(binning_index) == 15,
          all(binning_index$truth_form == binning_index$fitted_form))

#-------------------------------------------------------------------------------
# Storage for binning-fit manifest
#-------------------------------------------------------------------------------

binning_manifest <- vector("list", nrow(binning_index))

#-------------------------------------------------------------------------------
# Run binning fits
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(binning_index))) {
  
  # Identify generating form, fitted form and binning scheme
  truth_i <- binning_index$truth_form[i]
  fitted_i <- binning_index$fitted_form[i]
  binning_i <- binning_index$binning[i]
  
  # Prepare the fixed pilot catalogue for fitting
  catalogue_i <- prepare_temporal_catalogue(pilots[[truth_i]]$catalogue,
                                            T1 = T_fit_start,
                                            T2 = T_fit_end)
  
  # Extract candidate temporal-binning settings
  bin_i <- binning_schemes[[binning_i]]
  
  outfile_i <- file.path(binning_fit_dir, paste0("truth_", truth_i, "_fit_", 
                                             fitted_i, "_", binning_i, ".rds"))
  
  message("Truth: ", truth_i, " | Fit: ", fitted_i, " | Binning: ", binning_i)
  
  if (file.exists(outfile_i)) {
    
    # Reuse completed fit if present
    message("Already exists: ", basename(outfile_i), " — skipping fit")
    
    obj_i <- readRDS(outfile_i)
    
    diag_i <- fit_diagnostics(obj_i$fit)
    usable_i <- fit_is_usable(diag_i)
    runtime_i <- obj_i$runtime_minutes
    
  } else {
    
    # Fit correctly specified ETAS model under the candidate binning scheme
    result_i <- fit_temporal_etas(catalogue = catalogue_i,
                                  fitted_form = fitted_i,
                                  binning = bin_i,
                                  fit_control = fit_control,
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
                 binning = binning_i,
                 binning_parameters = bin_i,
                 runtime_minutes = runtime_i,
                 diagnostics = diag_i,
                 usable = usable_i), outfile_i)
    
    message("Saved ", basename(outfile_i), " | runtime = ",
            round(runtime_i, 2), " min")
  }
  
  # Record one-row summary for the binning manifest
  binning_manifest[[i]] <- data.frame(truth_form = truth_i,
                                      fitted_form = fitted_i,
                                      binning = binning_i,
                                      coef_t = bin_i$coef.t,
                                      delta_t = bin_i$delta.t,
                                      N_max = bin_i$N.max,
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

binning_manifest <- do.call(rbind, binning_manifest)
rownames(binning_manifest) <- NULL

manifest_file <- file.path(binning_fit_dir, "binning_manifest.csv")

write.csv(binning_manifest, manifest_file, row.names = FALSE)

stopifnot(nrow(binning_manifest) == 15)

message("Saved ", manifest_file)
message("Finished all 15 binning-check fits.")