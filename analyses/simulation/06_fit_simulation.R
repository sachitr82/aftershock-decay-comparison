#===============================================================================
# Fit final synthetic catalogue ensemble.
# WARNING: The full production run writes approximately 10–12 GB of fitted
# model objects to disk. On an Apple M4 Mac using 6 parallel workers, the full
# run took approximately 20 hours. Runtime will vary with available hardware.
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(future)
library(future.apply)
library(here)
library(inlabru)

source(here("analyses", "simulation", "00_design.R"))
source(here("src", "fit_helpers", "fit_diagnostics.R"))
source(here("src", "fit_helpers", "fit_temporal_etas.R"))

#-------------------------------------------------------------------------------
# Fixed production design
#-------------------------------------------------------------------------------

candidate_forms <- c("ou", "mse", "rate_state")

stopifnot(temporal_binning$N.max == 14,
          fit_control$rel_tol == 0.1,
          fit_control$max_iter == 100)

#-------------------------------------------------------------------------------
# Load final simulation manifest (from 05)
#-------------------------------------------------------------------------------

simulation_manifest <- read.csv(file.path(catalogue_dir, "simulation_manifest.csv"),
                                stringsAsFactors = FALSE)

simulation_manifest <- simulation_manifest[
  order(match(simulation_manifest$form, candidate_forms), 
        simulation_manifest$rep), ]

rownames(simulation_manifest) <- NULL
simulation_manifest$catalogue_file <- file.path(catalogue_dir, simulation_manifest$file)

stopifnot(
  nrow(simulation_manifest) == 3 * n_rep,
  all(table(simulation_manifest$form) == n_rep),
  !anyDuplicated(simulation_manifest[, c("form", "rep")]),
  all(file.exists(simulation_manifest$catalogue_file)))

n_fits <- nrow(simulation_manifest) * length(candidate_forms)

message(nrow(simulation_manifest), " catalogues | ",
        length(candidate_forms), " fitted forms | ", 
        n_fits, " total fits")

#-------------------------------------------------------------------------------
# Output directories for fits
#-------------------------------------------------------------------------------

for (form in candidate_forms) {
  dir.create(file.path(fit_dir, paste0("truth_", form)),
             recursive = TRUE, showWarnings = FALSE)
}

#-------------------------------------------------------------------------------
# Fixed reproducible fit seeds, assign one to each fit/catalogue combination
#-------------------------------------------------------------------------------

fit_seed_base <- 1000000

fit_index <- do.call(rbind, lapply(seq_len(nrow(simulation_manifest)), function(i) {
  
  data.frame(
    truth_form = simulation_manifest$form[i],
    rep = simulation_manifest$rep[i],
    catalogue_seed = simulation_manifest$seed[i],
    fitted_form = candidate_forms,
    stringsAsFactors = FALSE
  )
}))

rownames(fit_index) <- NULL
fit_index$fit_seed <- fit_seed_base + seq_len(nrow(fit_index))

stopifnot(
  nrow(fit_index) == n_fits,
  !anyDuplicated(fit_index$fit_seed),
  !anyDuplicated(fit_index[, c("truth_form", "rep", "fitted_form")])
)

write.csv(fit_index, file.path(fit_dir, "fit_index.csv"), row.names = FALSE)
#-------------------------------------------------------------------------------
# RS and MSE fits first
#-------------------------------------------------------------------------------

fit_order <- c("mse", "rate_state", "ou")
execution_order <- order(match(simulation_manifest$form, fit_order),
                         simulation_manifest$rep)
#-------------------------------------------------------------------------------
# Record production settings
#-------------------------------------------------------------------------------

saveRDS(list(candidate_forms = candidate_forms, 
             temporal_binning = temporal_binning,
             fit_control = fit_control,
             fit_seed_base = fit_seed_base,
             initials = initials, 
             prior_baseline = prior_baseline),
        file.path(fit_dir, "run_design.rds"))

#-------------------------------------------------------------------------------
# Fit all three candidate forms to one catalogue
#-------------------------------------------------------------------------------

fit_catalogue <- function(i) {
  
  # One numerical thread per worker; parallelism occurs between catalogues
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  
  INLA::inla.setOption(num.threads = "1:1")
  
  truth_i <- simulation_manifest$form[i]
  rep_i <- simulation_manifest$rep[i]
  catalogue_seed_i <- simulation_manifest$seed[i]
  
  #-----------------------------------------------------------------------------
  # Load catalogue once
  #-------------------------------------------------------------------------------
  
  obj_cat <- readRDS(simulation_manifest$catalogue_file[i])
  
  stopifnot(
    obj_cat$truth_form == truth_i,
    obj_cat$rep == rep_i,
    obj_cat$seed == catalogue_seed_i,
    isTRUE(all.equal(obj_cat$truth_parameters, truths[[truth_i]])),
    isTRUE(all.equal(obj_cat$design$M0, M0)),
    isTRUE(all.equal(obj_cat$design$Mmax, Mmax)),
    isTRUE(all.equal(obj_cat$design$T_fit_start, T_fit_start)),
    isTRUE(all.equal(obj_cat$design$T_fit_end, T_fit_end)))
  
  catalogue_i <- prepare_temporal_catalogue(obj_cat$catalogue,
                                            T1 = T_fit_start,
                                            T2 = T_fit_end)
  
  fit_rows <- vector("list", length(candidate_forms))
  
  #-----------------------------------------------------------------------------
  # Fit OU, MSE and rate-state
  #-----------------------------------------------------------------------------
  
  for (k in seq_along(fit_order)) {
    
    fitted_i <- fit_order[k]
    
    job_i <- fit_index[
      fit_index$truth_form == truth_i & 
        fit_index$rep == rep_i &
        fit_index$fitted_form == fitted_i, ]
    
    stopifnot(nrow(job_i) == 1)
    
    fit_seed_i <- job_i$fit_seed
    
    outfile_i <- file.path(fit_dir, paste0("truth_", truth_i),
                           sprintf("rep_%04d_fit_%s.rds", rep_i, fitted_i))
    
    message("Truth: ", truth_i, " | Rep: ", rep_i, " | Fit: ", fitted_i)
    
    fit_rows[[k]] <- tryCatch({
      
      #-------------------------------------------------------------------------
      # Reuse completed fit
      #-------------------------------------------------------------------------
      
      if (file.exists(outfile_i)) {
        
        existing_i <- readRDS(outfile_i)
        
        stopifnot(
          existing_i$truth_form == truth_i,
          existing_i$fitted_form == fitted_i,
          existing_i$rep == rep_i,
          existing_i$fit_seed == fit_seed_i,
          isTRUE(all.equal(existing_i$binning_parameters, temporal_binning)),
          isTRUE(all.equal(existing_i$rel_tol, fit_control$rel_tol)),
          existing_i$max_iter == fit_control$max_iter
        )
        
        message("Already exists: ", basename(outfile_i), " — skipping fit")
        
        diag_i <- fit_diagnostics(existing_i$fit)
        usable_i <- fit_is_usable(diag_i)
        
        row_i <- data.frame(truth_form = truth_i,
                            fitted_form = fitted_i,
                            rep = rep_i,
                            catalogue_seed = catalogue_seed_i,
                            fit_seed = fit_seed_i,
                            n_fit = existing_i$n_fit, 
                            runtime_minutes = existing_i$runtime_minutes,
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
                            status = "existing",
                            error = NA_character_, 
                            file = basename(outfile_i))
        rm(existing_i)                  
      } else {
        
        #-----------------------------------------------------------------------
        # New fit
        #-----------------------------------------------------------------------
        
        result_i <- fit_temporal_etas(catalogue = catalogue_i,
                                      fitted_form = fitted_i,
                                      binning = temporal_binning,
                                      fit_control = fit_control,
                                      M0 = M0,
                                      T1 = T_fit_start,
                                      T2 = T_fit_end,
                                      compute_model_criteria = TRUE,
                                      seed = fit_seed_i)
        
        fit_i <- result_i$fit
        link_i <- result_i$link.functions
        runtime_i <- result_i$runtime_minutes
        diag_i <- result_i$diagnostics
        usable_i <- result_i$usable
        
        missing_i <- setdiff(c("dic", "mlik"), names(fit_i))
        
        if (length(missing_i) > 0) {
          stop("Missing requested INLA output: ", paste(missing_i, collapse = ", "))
        }
        
        #-----------------------------------------------------------------------
        # Save completed fit
        #-----------------------------------------------------------------------
        
        output_i <- list(fit = fit_i,
                         link.functions = link_i,
                         truth_form = truth_i,
                         truth_parameters = obj_cat$truth_parameters,
                         fitted_form = fitted_i,
                         rep = rep_i,
                         catalogue_seed = catalogue_seed_i,
                         fit_seed = fit_seed_i,
                         n_fit = nrow(catalogue_i),
                         binning_parameters = temporal_binning,
                         rel_tol = fit_control$rel_tol,
                         max_iter = fit_control$max_iter,
                         runtime_minutes = runtime_i,
                         diagnostics = diag_i,
                         usable = usable_i)
        
        # Temporary file prevents an interrupted save appearing as a completed fit
        tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i), fileext = ".rds")
        saveRDS(output_i, tmp_i)
        
        if (!file.rename(tmp_i, outfile_i)) {
          stop("Could not move temporary fit file to final output.")
        }
        
        message("Saved ", basename(outfile_i), 
                " | runtime = ", round(runtime_i, 2), " min")
        
        row_i <- data.frame(truth_form = truth_i, 
                            fitted_form = fitted_i,
                            rep = rep_i,
                            catalogue_seed = catalogue_seed_i,
                            fit_seed = fit_seed_i,
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
                            status = "fitted", 
                            error = NA_character_, 
                            file = basename(outfile_i))
        
        rm(fit_i, result_i, output_i)
      }
      
      gc(verbose = FALSE)
      row_i
      
    }, error = function(e) {
      
      gc(verbose = FALSE)
      
      error_i <- conditionMessage(e)
      inla_failure_i <- grepl(
        paste("Problem in inla",
              "Giving up and returning last successfully obtained result",
              "inla-program exited with an error",
              "maximum number of tries has been reached",
              "Newton-Raphson optimizer did not converge",
              sep = "|"), error_i)
      
      data.frame(truth_form = truth_i, 
                 fitted_form = fitted_i, 
                 rep = rep_i, 
                 catalogue_seed = catalogue_seed_i,
                 fit_seed = fit_seed_i,
                 n_fit = nrow(catalogue_i), 
                 runtime_minutes = NA_real_,
                 fit_status = ifelse(inla_failure_i, "inla_failure", "error"),
                 converged = FALSE, hit_max = NA, inla_failure = inla_failure_i,
                 nan_inf_logl = NA, vb_aborted = NA, n_iter = NA_integer_, 
                 min_sd = NA_real_, degenerate = NA, usable = FALSE,
                 status = "error", error = error_i,
                 file = basename(outfile_i))
    })
  }
  
  do.call(rbind, fit_rows)
}

#-------------------------------------------------------------------------------
# Parallel execution
#-------------------------------------------------------------------------------

# Six parallel workers, with one numerical thread per worker. Adjust according to
# available CPU cores and memory.

n_workers <- min(6, future::availableCores())

message("Running ", n_fits, " fits using ", n_workers, " parallel workers.")

future::plan(future::multisession, workers = n_workers)

fit_rows <- future.apply::future_lapply(
  execution_order, fit_catalogue,
  future.seed = 12345,
  future.packages = c("ETAS.inlabru", "inlabru", "INLA"),
  future.scheduling = 25)

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Save fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- do.call(rbind, fit_rows)
rownames(fit_manifest) <- NULL

fit_manifest <- fit_manifest[
  order(match(fit_manifest$truth_form, candidate_forms),
        fit_manifest$rep,
        match(fit_manifest$fitted_form, candidate_forms)), ]

write.csv(fit_manifest, file.path(fit_dir, "fit_manifest.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Problem fits
#-------------------------------------------------------------------------------

problem_fits <- fit_manifest[fit_manifest$status == "error" |
                             is.na(fit_manifest$usable) |
                             !fit_manifest$usable, ]

write.csv(problem_fits, file.path(fit_dir, "problem_fits.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Final completeness checks
#-------------------------------------------------------------------------------

stopifnot(nrow(fit_manifest) == n_fits, 
          !anyDuplicated(fit_manifest[, c("truth_form", "rep", "fitted_form")]))

message("Usable fits: ", sum(fit_manifest$usable %in% TRUE), "/", n_fits)
message("Problem fits: ", nrow(problem_fits), "/", n_fits)

message("Finished main simulation fitting: ",
        sum(fit_manifest$status != "error"), "/", n_fits,
        " fits available.")