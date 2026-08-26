#===============================================================================
# Temporal binning check
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

pilot_files <- c(
  ou = file.path(pilot_dir, "ou_pilot.rds"),
  mse = file.path(pilot_dir, "mse_pilot.rds"),
  rate_state = file.path(pilot_dir, "rate_state_pilot.rds"))

pilots <- lapply(pilot_files, readRDS)

#-------------------------------------------------------------------------------
# Candidate fitted kernels
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Targeted N.max refinement
#-------------------------------------------------------------------------------

binning_schemes <- list(
  N8  = list(coef.t = 1, delta.t = 0.1, N.max = 8),
  N10 = list(coef.t = 1, delta.t = 0.1, N.max = 10),
  N12 = list(coef.t = 1, delta.t = 0.1, N.max = 12),
  N14 = list(coef.t = 1, delta.t = 0.1, N.max = 14),
  N16 = list(coef.t = 1, delta.t = 0.1, N.max = 16)
)

#-------------------------------------------------------------------------------
# Fit index: one correctly specified fit per truth under N8, N10 and N12 (15 fits)
#-------------------------------------------------------------------------------

binning_index <- do.call(rbind, lapply(candidate_kernels, function(kernel) {
  data.frame(
    truth_kernel = kernel, fitted_kernel = kernel,
    binning = names(binning_schemes), stringsAsFactors = FALSE)
}))

rownames(binning_index) <- NULL

stopifnot(nrow(binning_index) == 15,
          all(binning_index$truth_kernel == binning_index$fitted_kernel))

#-------------------------------------------------------------------------------
# Output directory
#-------------------------------------------------------------------------------

binning_dir <- file.path(fit_dir, "binning_check")
dir.create(binning_dir, recursive = TRUE, showWarnings = FALSE)

binning_manifest <- vector("list", nrow(binning_index))

#-------------------------------------------------------------------------------
# Run binning fits
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(binning_index))) {
  
  truth_i <- binning_index$truth_kernel[i]
  fitted_i <- binning_index$fitted_kernel[i]
  binning_i <- binning_index$binning[i]
  
  catalogue_i <- pilots[[truth_i]]$catalogue
  catalogue_i <- catalogue_i[catalogue_i$ts >= T_fit_start & 
                               catalogue_i$ts <= T_fit_end,]
  catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  bin_i <- binning_schemes[[binning_i]]
  
  outfile_i <- file.path(binning_dir, paste0("truth_", truth_i, "_fit_", 
                                             fitted_i, "_", binning_i, ".rds"))
  
  message("Truth: ", truth_i, " | Fit: ", fitted_i, " | Binning: ", binning_i)
  
  if (file.exists(outfile_i)) {
    
    message("Already exists: ", basename(outfile_i), " — skipping fit")
    
    obj_i <- readRDS(outfile_i)
    
    runtime_i <- obj_i$runtime_minutes
    converged_i <- obj_i$converged
    hit_max_i <- obj_i$hit_max
    n_iter_i <- obj_i$n_iter
  } else{
  
    link_i <- make_links_P0(fitted_i)
    bru_i <- make_bru_options_P0(fitted_i, rel_tol = 0.1, max_iter = 100)
    
    start_i <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(
      total.data = catalogue_i, M0 = M0, T1 = T_fit_start, T2 = T_fit_end,
      link.functions = link_i, coef.t. = bin_i$coef.t, delta.t. = bin_i$delta.t,
      N.max. = bin_i$N.max, bru.opt = bru_i, kernel = fitted_i)
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
    log_i <- as.character(inlabru::bru_log(fit_i))
    converged_i <- any(grepl("Convergence criterion met", log_i, fixed = TRUE))
    hit_max_i <- any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
    n_iter_i <- max(fit_i$bru_iinla$track$iteration, na.rm = TRUE)
    
    saveRDS(list(
      fit = fit_i,
      link.functions=link_i,
      truth_kernel = truth_i,
      fitted_kernel = fitted_i,
      binning = binning_i,
      binning_parameters = bin_i,
      runtime_minutes = runtime_i,
      converged = converged_i,
      hit_max = hit_max_i,
      n_iter = n_iter_i), outfile_i)
    
    message("Saved ", basename(outfile_i), " | runtime = ", round(runtime_i, 2),
            " min")
  }
  binning_manifest[[i]] <- data.frame(
    truth_kernel = truth_i, fitted_kernel = fitted_i, binning = binning_i,
    coef_t = bin_i$coef.t, delta_t = bin_i$delta.t, N_max = bin_i$N.max,
    n_fit = nrow(catalogue_i),
    runtime_minutes = runtime_i, converged = converged_i,
    hit_max = hit_max_i, n_iter = n_iter_i,
    file = basename(outfile_i))

}

#-------------------------------------------------------------------------------
# Save fit manifest
#-------------------------------------------------------------------------------

binning_manifest <- do.call(rbind, binning_manifest)
rownames(binning_manifest) <- NULL

write.csv(binning_manifest, file.path(binning_dir, "binning_manifest.csv"), row.names = FALSE)

stopifnot(nrow(binning_manifest) == 15)

message("Finished all 15 binning-check fits.")