#===============================================================================
# rel_tol convergence check
#===============================================================================

#-------------------------------------------------------------------------------
# Load package and design choices
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(future)
library(here)
library(inlabru)

source(here::here("analyses", "simulation", "00_design.R"))

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Load fixed pilot catalogues (from 01)
#-------------------------------------------------------------------------------

pilot_files <- c(ou = file.path(pilot_dir, "ou_pilot.rds"),
                 mse = file.path(pilot_dir, "mse_pilot.rds"),
                 rate_state = file.path(pilot_dir, "rate_state_pilot.rds"))

pilots <- lapply(pilot_files, readRDS)

#-------------------------------------------------------------------------------
# Candidate fitted kernels
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Targeted rel_tol refinement
#-------------------------------------------------------------------------------

rel_tol_schemes <- c(RT01 = 0.1, RT005 = 0.05, RT001 = 0.01)

#-------------------------------------------------------------------------------
# Fit index: one correctly specified fit per truth under each rel_tol (9 fits)
#-------------------------------------------------------------------------------

rel_tol_index <- do.call(rbind, lapply(candidate_kernels, function(kernel) {
  data.frame(truth_kernel = kernel, fitted_kernel = kernel,
             tolerance = names(rel_tol_schemes), stringsAsFactors = FALSE)
}))

rownames(rel_tol_index) <- NULL

stopifnot(nrow(rel_tol_index) == 9,
          all(rel_tol_index$truth_kernel == rel_tol_index$fitted_kernel))

#-------------------------------------------------------------------------------
# Output directory
#-------------------------------------------------------------------------------

rel_tol_dir <- file.path(fit_dir, "rel_tol_check")
dir.create(rel_tol_dir, recursive = TRUE, showWarnings = FALSE)

rel_tol_manifest <- vector("list", nrow(rel_tol_index))

#-------------------------------------------------------------------------------
# Run rel_tol fits
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(rel_tol_index))) {
  
  truth_i <- rel_tol_index$truth_kernel[i]
  fitted_i <- rel_tol_index$fitted_kernel[i]
  tolerance_i <- rel_tol_index$tolerance[i]
  rel_tol_i <- rel_tol_schemes[[tolerance_i]]
  
  catalogue_i <- pilots[[truth_i]]$catalogue
  catalogue_i <- catalogue_i[catalogue_i$ts >= T_fit_start &
                               catalogue_i$ts <= T_fit_end, ]
  catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  outfile_i <- file.path(rel_tol_dir, paste0("truth_", truth_i, "_fit_", 
                                             fitted_i, "_", tolerance_i, ".rds"))
  
  message("Truth: ", truth_i, " | Fit: ", fitted_i, " | rel_tol: ", rel_tol_i)
  
  if (file.exists(outfile_i)) {
    
    message("Already exists: ", basename(outfile_i), " — skipping fit")
    
    obj_i <- readRDS(outfile_i)
  
    runtime_i <- obj_i$runtime_minutes
    converged_i <- obj_i$converged
    hit_max_i <- obj_i$hit_max
    n_iter_i <- obj_i$n_iter
  } else {
    
    link_i <- make_links_P0(fitted_i)
    bru_i <- make_bru_options_P0(fitted_i, rel_tol = rel_tol_i, 
                                 max_iter = 100)
    
    start_i <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(
      total.data = catalogue_i, 
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
    
    log_i <- as.character(inlabru::bru_log(fit_i))
    converged_i <- any(grepl("Convergence criterion met", log_i, fixed = TRUE))
    hit_max_i <- any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
    n_iter_i <- max(fit_i$bru_iinla$track$iteration, na.rm = TRUE)
    
    saveRDS(list(fit = fit_i,
                 link.functions = link_i,
                 truth_kernel = truth_i,
                 fitted_kernel = fitted_i,
                 tolerance = tolerance_i,
                 rel_tol = rel_tol_i,
                 binning_parameters = temporal_binning,
                 runtime_minutes = runtime_i,
                 converged = converged_i,
                 hit_max = hit_max_i,
                 n_iter = n_iter_i),
                 outfile_i)
              
    message("Saved ", basename(outfile_i),
            " | runtime = ", round(runtime_i, 2), " min")
  }
  
  rel_tol_manifest[[i]] <- data.frame(truth_kernel = truth_i,
                                      fitted_kernel = fitted_i,
                                      tolerance = tolerance_i,
                                      rel_tol = rel_tol_i,
                                      N_max = temporal_binning$N.max,
                                      n_fit = nrow(catalogue_i),
                                      runtime_minutes = runtime_i,
                                      converged = converged_i,
                                      hit_max = hit_max_i,
                                      n_iter = n_iter_i,
                                      file = basename(outfile_i))
}

#-------------------------------------------------------------------------------
# Save fit manifest
#-------------------------------------------------------------------------------

rel_tol_manifest <- do.call(rbind, rel_tol_manifest)
rownames(rel_tol_manifest) <- NULL

write.csv(rel_tol_manifest, file.path(rel_tol_dir, "rel_tol_manifest.csv"),
          row.names = FALSE)

stopifnot(nrow(rel_tol_manifest) == 9)

message("Finished all 9 rel_tol-check fits.")