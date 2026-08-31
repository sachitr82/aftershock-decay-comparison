#===============================================================================
# Fit final synthetic catalogue ensemble
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(future)
library(future.apply)
library(here)
library(inlabru)

source(here::here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Fixed production design
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

stopifnot(temporal_binning$N.max == 14,
          fit_control$rel_tol == 0.1,
          fit_control$max_iter == 100
)

#-------------------------------------------------------------------------------
# Load final simulation manifest (from 03)
#-------------------------------------------------------------------------------

simulation_manifest <- read.csv(file.path(catalogue_dir, "simulation_manifest.csv"),
                                stringsAsFactors = FALSE)

simulation_manifest <- simulation_manifest[
  order(match(simulation_manifest$kernel, candidate_kernels), 
        simulation_manifest$rep), ]

rownames(simulation_manifest) <- NULL
simulation_manifest$catalogue_file <- file.path(catalogue_dir, simulation_manifest$file)

stopifnot(
  nrow(simulation_manifest) == 3 * n_rep,
  all(table(simulation_manifest$kernel) == n_rep),
  !anyDuplicated(simulation_manifest[, c("kernel", "rep")]),
  all(file.exists(simulation_manifest$catalogue_file))
)

n_fits <- nrow(simulation_manifest) * length(candidate_kernels)

message(nrow(simulation_manifest), " catalogues | ",
        length(candidate_kernels), " fitted kernels | ", 
        n_fits, " total fits")

#-------------------------------------------------------------------------------
# Output directory
#-------------------------------------------------------------------------------

main_fit_dir <- file.path(fit_dir, "main_simulation")
dir.create(main_fit_dir, recursive = TRUE, showWarnings = FALSE)

for (kernel in candidate_kernels) {
  dir.create(file.path(main_fit_dir, paste0("truth_", kernel)),
             recursive = TRUE, showWarnings = FALSE)
}

#-------------------------------------------------------------------------------
# Fixed reproducible fit seeds
#-------------------------------------------------------------------------------

fit_seed_base <- 1000000

fit_index <- do.call(rbind, lapply(seq_len(nrow(simulation_manifest)), function(i) {
  
  data.frame(
    truth_kernel = simulation_manifest$kernel[i],
    rep = simulation_manifest$rep[i],
    catalogue_seed = simulation_manifest$seed[i],
    fitted_kernel = candidate_kernels,
    stringsAsFactors = FALSE
  )
}))

rownames(fit_index) <- NULL
fit_index$fit_seed <- fit_seed_base + seq_len(nrow(fit_index))

stopifnot(
  nrow(fit_index) == n_fits,
  !anyDuplicated(fit_index$fit_seed),
  !anyDuplicated(fit_index[, c("truth_kernel", "rep", "fitted_kernel")])
)

write.csv(fit_index, file.path(main_fit_dir, "fit_index.csv"), row.names = FALSE)
#-------------------------------------------------------------------------------
# RS and MSE fits first
#-------------------------------------------------------------------------------

fit_order <- c("mse", "rate_state", "ou")
execution_order <- order(match(simulation_manifest$kernel, c("mse", "rate_state", "ou")), simulation_manifest$rep)

#-------------------------------------------------------------------------------
# Record production settings
#-------------------------------------------------------------------------------

saveRDS(
  list(candidate_kernels = candidate_kernels, temporal_binning = temporal_binning,
       fit_control = fit_control, fit_seed_base = fit_seed_base,
       initials = initials, prior_baseline = prior_baseline),
  file.path(main_fit_dir, "run_design.rds")
)

writeLines(capture.output(sessionInfo()), file.path(main_fit_dir, "sessionInfo.txt"))

#-------------------------------------------------------------------------------
# Fit diagnostics
#-------------------------------------------------------------------------------

fit_diagnostics <- function(fit) {
  
  log_i <- as.character(inlabru::bru_log(fit))
  
  inla_failure <- any(grepl(paste("Problem in inla",
                                  "Giving up and returning last successfully obtained result",
                                  "inla-program exited with an error",
                                  "maximum number of tries has been reached",
                                  "Newton-Raphson optimizer did not converge", sep = "|"), log_i))
  
  nan_inf_logl <- any(grepl("NAN/INF values in logl", log_i, fixed = TRUE))
  vb_aborted <- any(grepl("max_correction|vb.correction.*aborted", log_i))
  
  converged <- !inla_failure && any(grepl("Convergence criterion met", log_i, fixed = TRUE))
  hit_max <- !inla_failure && any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
  
  fit_status <- dplyr::case_when(inla_failure ~ "inla_failure",
                                 converged ~ "converged",
                                 hit_max ~ "hit_max",
                                 TRUE ~ "unknown")
  
  n_iter <- max(fit$bru_iinla$track$iteration, na.rm = TRUE)
  
  list(fit_status = fit_status, converged = converged, hit_max = hit_max,
       inla_failure = inla_failure, nan_inf_logl = nan_inf_logl,
       vb_aborted = vb_aborted, n_iter = n_iter)
}
#-------------------------------------------------------------------------------
# Fit all three candidate kernels to one catalogue
#-------------------------------------------------------------------------------

fit_catalogue <- function(i) {
  
  # One numerical thread per worker; parallelism occurs between catalogues
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  
  INLA::inla.setOption(num.threads = "1:1")
  
  truth_i <- simulation_manifest$kernel[i]
  rep_i <- simulation_manifest$rep[i]
  catalogue_seed_i <- simulation_manifest$seed[i]
  
  #-----------------------------------------------------------------------------
  # Load catalogue once
  #-------------------------------------------------------------------------------
  
  obj_cat <- readRDS(simulation_manifest$catalogue_file[i])
  
  stopifnot(
    obj_cat$truth_kernel == truth_i,
    obj_cat$rep == rep_i,
    obj_cat$seed == catalogue_seed_i,
    isTRUE(all.equal(obj_cat$truth_parameters, truths[[truth_i]])),
    isTRUE(all.equal(obj_cat$design$M0, M0)),
    isTRUE(all.equal(obj_cat$design$Mmax, Mmax)),
    isTRUE(all.equal(obj_cat$design$T_fit_start, T_fit_start)),
    isTRUE(all.equal(obj_cat$design$T_fit_end, T_fit_end)))
  
  catalogue_i <- obj_cat$catalogue
  catalogue_i <- catalogue_i[catalogue_i$ts >= T_fit_start &
                               catalogue_i$ts <= T_fit_end, ]
  catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  fit_rows <- vector("list", length(candidate_kernels))
  
  #-----------------------------------------------------------------------------
  # Fit OU, MSE and rate-state
  #-----------------------------------------------------------------------------
  
  for (k in seq_along(fit_order)) {
    
    fitted_i <- fit_order[k]
    
    job_i <- fit_index[
      fit_index$truth_kernel == truth_i & 
        fit_index$rep == rep_i &
        fit_index$fitted_kernel == fitted_i, ]
    
    stopifnot(nrow(job_i) == 1)
    
    fit_seed_i <- job_i$fit_seed
    
    outfile_i <- file.path(main_fit_dir, paste0("truth_", truth_i),
                           sprintf("rep_%04d_fit_%s.rds", rep_i, fitted_i))
    
    message("Truth: ", truth_i, " | Rep: ", rep_i, " | Fit: ", fitted_i)
    
    fit_rows[[k]] <- tryCatch({
      
      #-------------------------------------------------------------------------
      # Reuse completed fit
      #-------------------------------------------------------------------------
      
      if (file.exists(outfile_i)) {
        
        existing_i <- readRDS(outfile_i)
        
        stopifnot(
          existing_i$truth_kernel == truth_i,
          existing_i$fitted_kernel == fitted_i,
          existing_i$rep == rep_i,
          existing_i$fit_seed == fit_seed_i,
          isTRUE(all.equal(existing_i$binning_parameters, temporal_binning)),
          isTRUE(all.equal(existing_i$rel_tol, fit_control$rel_tol)),
          existing_i$max_iter == fit_control$max_iter
        )
        
        message("Already exists: ", basename(outfile_i), " — skipping fit")
        diag_i <- fit_diagnostics(existing_i$fit)
        min_sd_i <- min(existing_i$fit$summary.fixed$sd, na.rm = TRUE)
        
        row_i <- data.frame(
          truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
          catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
          n_fit = existing_i$n_fit, runtime_minutes = existing_i$runtime_minutes,
          fit_status = diag_i$fit_status, converged = diag_i$converged,
          hit_max = diag_i$hit_max, inla_failure = diag_i$inla_failure,
          nan_inf_logl = diag_i$nan_inf_logl, vb_aborted = diag_i$vb_aborted,
          n_iter = diag_i$n_iter, min_sd = min_sd_i, degenerate = isTRUE(min_sd_i < 1e-6),
          status = "existing", error = NA_character_, file = basename(outfile_i))
        
      } else {
        
        #-----------------------------------------------------------------------
        # New fit
        #-----------------------------------------------------------------------
        
        set.seed(fit_seed_i)
        
        link_i <- make_links_P0(fitted_i)
        bru_i <- make_bru_options_P0(fitted_i, rel_tol = fit_control$rel_tol,
                                     max_iter = fit_control$max_iter)
        bru_i$control.compute <- list(config = TRUE, dic = TRUE, waic = TRUE
                                      , mlik = TRUE)
        
        start_i <- Sys.time()
        
        fit_i <- ETAS.inlabru::Temporal.ETAS(
          total.data = catalogue_i, M0 = M0, T1 = T_fit_start, T2 = T_fit_end,
          link.functions = link_i, coef.t. = temporal_binning$coef.t,
          delta.t. = temporal_binning$delta.t, N.max. = temporal_binning$N.max,
          bru.opt = bru_i, kernel = fitted_i
        )
        
        runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
        
        diag_i <- fit_diagnostics(fit_i)
        
        fit_status_i <- diag_i$fit_status
        converged_i <- diag_i$converged
        hit_max_i <- diag_i$hit_max
        inla_failure_i <- diag_i$inla_failure
        nan_inf_i <- diag_i$nan_inf_logl
        vb_aborted_i <- diag_i$vb_aborted
        n_iter_i <- diag_i$n_iter
        min_sd_i <- min(fit_i$summary.fixed$sd, na.rm = TRUE)
        
        missing_i <- setdiff(c("dic", "waic", "mlik"), names(fit_i))
        
        if (length(missing_i) > 0) {
          stop("Missing requested INLA output: ", paste(missing_i, collapse = ", "))
        }
        
        #-----------------------------------------------------------------------
        # Save completed fit
        #-----------------------------------------------------------------------
        
        output_i <- list(
          fit = fit_i, link.functions = link_i,
          truth_kernel = truth_i, truth_parameters = obj_cat$truth_parameters,
          fitted_kernel = fitted_i, rep = rep_i,
          catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
          n_fit = nrow(catalogue_i), binning_parameters = temporal_binning,
          rel_tol = fit_control$rel_tol, max_iter = fit_control$max_iter,
          runtime_minutes = runtime_i, fit_status = fit_status_i,
          converged = converged_i, hit_max = hit_max_i,
          inla_failure = inla_failure_i, nan_inf_logl = nan_inf_i,
          vb_aborted = vb_aborted_i, n_iter = n_iter_i)
        
        # Temporary file prevents an interrupted save appearing as a completed fit
        tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i), fileext = ".rds")
        saveRDS(output_i, tmp_i)
        
        if (!file.rename(tmp_i, outfile_i)) {
          stop("Could not move temporary fit file to final output.")
        }
        
        message("Saved ", basename(outfile_i), 
                " | runtime = ", round(runtime_i, 2), " min")
        
        row_i <- data.frame(
          truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
          catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
          n_fit = nrow(catalogue_i), runtime_minutes = runtime_i,
          fit_status = fit_status_i, converged = converged_i, hit_max = hit_max_i,
          inla_failure = inla_failure_i, nan_inf_logl = nan_inf_i,
          vb_aborted = vb_aborted_i, n_iter = n_iter_i, 
          min_sd = min_sd_i, degenerate = isTRUE(min_sd_i < 1e-6), 
          status = "fitted", error = NA_character_,
          file = basename(outfile_i))
        
        rm(fit_i, output_i)
      }
      
      gc(verbose = FALSE)
      row_i
      
    }, error = function(e) {
      
      gc(verbose = FALSE)
      
      error_i <- conditionMessage(e)
      inla_failure_i <- grepl("Problem in inla|inla-program exited|Newton-Raphson optimizer did not converge",
                              error_i)
      
      data.frame(
        truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
        catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
        n_fit = nrow(catalogue_i), runtime_minutes = NA_real_,
        fit_status = ifelse(inla_failure_i, "inla_failure", "error"),
        converged = FALSE, hit_max = NA, inla_failure = inla_failure_i,
        nan_inf_logl = NA, vb_aborted = NA, n_iter = NA_integer_, 
        min_sd = NA_real_, degenerate = NA,
        status = "error", error = error_i,
        file = basename(outfile_i))
    })
  }
  
  do.call(rbind, fit_rows)
}

#-------------------------------------------------------------------------------
# Parallel execution
#-------------------------------------------------------------------------------

# Six parallel workers, with one numerical thread per worker
n_workers <- min(6, future::availableCores())

message("Running ", n_fits, " fits using ", n_workers, " parallel workers.")

future::plan(future::multisession, workers = n_workers)

fit_rows <- future.apply::future_lapply(
  execution_order, fit_catalogue,
  future.seed = 12345,
  future.packages = c("ETAS.inlabru", "inlabru", "INLA"),
  future.scheduling = 25
)

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Save fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- do.call(rbind, fit_rows)
rownames(fit_manifest) <- NULL

fit_manifest <- fit_manifest[
  order(match(fit_manifest$truth_kernel, candidate_kernels),
        fit_manifest$rep,
        match(fit_manifest$fitted_kernel, candidate_kernels)), ]

write.csv(fit_manifest, file.path(main_fit_dir, "fit_manifest.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Problem fits
#-------------------------------------------------------------------------------

problem_fits <- fit_manifest[
  fit_manifest$status == "error" |
    fit_manifest$fit_status %in% c("hit_max", "inla_failure", "unknown") |
    fit_manifest$hit_max %in% TRUE |
    fit_manifest$degenerate %in% TRUE |
    fit_manifest$vb_aborted %in% TRUE |
    fit_manifest$nan_inf_logl %in% TRUE, ]

write.csv(problem_fits, file.path(main_fit_dir, "problem_fits.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Final completeness checks
#-------------------------------------------------------------------------------

print(table(fit_manifest$status))
print(table(fit_manifest$converged, useNA = "ifany"))
print(table(fit_manifest$fitted_kernel, fit_manifest$degenerate, useNA = "ifany"))
message("Usable fits: ", sum(fit_manifest$converged & !fit_manifest$hit_max &
                               !fit_manifest$degenerate & 
                               !fit_manifest$vb_aborted &
                               !fit_manifest$inla_failure & 
                               !fit_manifest$nan_inf_logl, 
                             na.rm = TRUE), "/", n_fits)
stopifnot(
  nrow(fit_manifest) == n_fits,
  !anyDuplicated(fit_manifest[, c("truth_kernel", "rep", "fitted_kernel")])
)

message("Finished main simulation fitting: ",
        sum(fit_manifest$status != "error"), "/", n_fits, " fits available.")