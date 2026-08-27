#===============================================================================
# Rescue non-converged simulation fits
#===============================================================================

library(ETAS.inlabru)
library(future)
library(future.apply)
library(here)
library(inlabru)

source(here::here("analyses", "simulation", "00_design.R"))

main_fit_dir <- file.path(fit_dir, "main_simulation")
rescue_dir <- file.path(main_fit_dir, "rescue_maxiter200")
dir.create(rescue_dir, recursive = TRUE, showWarnings = FALSE)

rescue_max_iter <- 200

# Original production failures
problem_fits <- read.csv(file.path(main_fit_dir, "problem_fits.csv"),
                         stringsAsFactors = FALSE)

stopifnot(
  nrow(problem_fits) == 114,
  all(problem_fits$status != "error"),
  all(!problem_fits$converged),
  all(problem_fits$hit_max),
  all(problem_fits$n_iter == 100)
)

# Catalogue manifest
simulation_manifest <- read.csv(file.path(catalogue_dir, "simulation_manifest.csv"),
                                stringsAsFactors = FALSE)

simulation_manifest$catalogue_file <- file.path(catalogue_dir,
                                                simulation_manifest$file)

for (kernel in c("ou", "mse", "rate_state")) {
  dir.create(file.path(rescue_dir, paste0("truth_", kernel)),
             recursive = TRUE, showWarnings = FALSE)
}

#-------------------------------------------------------------------------------
# Fit one failed job
#-------------------------------------------------------------------------------

rescue_fit <- function(i) {
  
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  INLA::inla.setOption(num.threads = "1:1")
  
  job <- problem_fits[i, ]
  
  truth_i <- job$truth_kernel
  fitted_i <- job$fitted_kernel
  rep_i <- job$rep
  catalogue_seed_i <- job$catalogue_seed
  fit_seed_i <- job$fit_seed
  
  cat_row <- simulation_manifest[
    simulation_manifest$kernel == truth_i &
      simulation_manifest$rep == rep_i, ]
  
  stopifnot(nrow(cat_row) == 1)
  
  obj_cat <- readRDS(cat_row$catalogue_file)
  
  stopifnot(
    obj_cat$truth_kernel == truth_i,
    obj_cat$rep == rep_i,
    obj_cat$seed == catalogue_seed_i
  )
  
  catalogue_i <- obj_cat$catalogue
  catalogue_i <- catalogue_i[
    catalogue_i$ts >= T_fit_start & catalogue_i$ts <= T_fit_end, ]
  catalogue_i <- catalogue_i[order(catalogue_i$ts), ]
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  outfile_i <- file.path(
    rescue_dir, paste0("truth_", truth_i),
    sprintf("rep_%04d_fit_%s.rds", rep_i, fitted_i)
  )
  
  # Resume safely if this rescue fit already exists
  if (file.exists(outfile_i)) {
    
    existing_i <- readRDS(outfile_i)
    
    stopifnot(
      existing_i$truth_kernel == truth_i,
      existing_i$fitted_kernel == fitted_i,
      existing_i$rep == rep_i,
      existing_i$fit_seed == fit_seed_i,
      existing_i$max_iter == rescue_max_iter,
      isTRUE(all.equal(existing_i$rel_tol, fit_control$rel_tol)),
      isTRUE(all.equal(existing_i$binning_parameters, temporal_binning))
    )
    
    return(data.frame(
      truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
      catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
      n_fit = existing_i$n_fit, runtime_minutes = existing_i$runtime_minutes,
      converged = existing_i$converged, hit_max = existing_i$hit_max,
      n_iter = existing_i$n_iter, n_cpo_fail = existing_i$n_cpo_fail,
      status = "existing", error = NA_character_
    ))
  }
  
  tryCatch({
    
    set.seed(fit_seed_i)
    
    link_i <- make_links_P0(fitted_i)
    bru_i <- make_bru_options_P0(
      fitted_i, rel_tol = fit_control$rel_tol, max_iter = rescue_max_iter
    )
    
    # Keep the same INLA outputs as the production run
    bru_i$control.compute <- list(
      config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE
    )
    
    start_i <- Sys.time()
    
    fit_i <- ETAS.inlabru::Temporal.ETAS(
      total.data = catalogue_i, M0 = M0, T1 = T_fit_start, T2 = T_fit_end,
      link.functions = link_i, coef.t. = temporal_binning$coef.t,
      delta.t. = temporal_binning$delta.t, N.max. = temporal_binning$N.max,
      bru.opt = bru_i, kernel = fitted_i
    )
    
    runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
    
    log_i <- as.character(inlabru::bru_log(fit_i))
    converged_i <- any(grepl("Convergence criterion met", log_i, fixed = TRUE))
    hit_max_i <- any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
    n_iter_i <- max(fit_i$bru_iinla$track$iteration, na.rm = TRUE)
    n_cpo_fail_i <- if (!is.null(fit_i$cpo$failure)) {
      sum(fit_i$cpo$failure != 0, na.rm = TRUE)
    } else {
      NA_integer_
    }
    
    output_i <- list(
      fit = fit_i,
      link.functions = link_i,
      truth_kernel = truth_i,
      truth_parameters = obj_cat$truth_parameters,
      fitted_kernel = fitted_i,
      rep = rep_i,
      catalogue_seed = catalogue_seed_i,
      fit_seed = fit_seed_i,
      n_fit = nrow(catalogue_i),
      binning_parameters = temporal_binning,
      rel_tol = fit_control$rel_tol,
      max_iter = rescue_max_iter,
      runtime_minutes = runtime_i,
      converged = converged_i,
      hit_max = hit_max_i,
      n_iter = n_iter_i,
      n_cpo_fail = n_cpo_fail_i
    )
    
    # Atomic save, as in 06
    tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i),
                      fileext = ".rds")
    saveRDS(output_i, tmp_i)
    
    if (!file.rename(tmp_i, outfile_i)) {
      stop("Could not move temporary fit file to final output.")
    }
    
    data.frame(
      truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
      catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
      n_fit = nrow(catalogue_i), runtime_minutes = runtime_i,
      converged = converged_i, hit_max = hit_max_i, n_iter = n_iter_i,
      n_cpo_fail = n_cpo_fail_i,
      status = "fitted", error = NA_character_
    )
    
  }, error = function(e) {
    
    data.frame(
      truth_kernel = truth_i, fitted_kernel = fitted_i, rep = rep_i,
      catalogue_seed = catalogue_seed_i, fit_seed = fit_seed_i,
      n_fit = nrow(catalogue_i), runtime_minutes = NA_real_,
      converged = FALSE, hit_max = NA, n_iter = NA_integer_,
      n_cpo_fail = NA_integer_, status = "error", error = conditionMessage(e)
    )
  })
}

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------

n_workers <- min(4, future::availableCores())
future::plan(future::multisession, workers = n_workers)

rescue_rows <- future.apply::future_lapply(
  seq_len(nrow(problem_fits)), rescue_fit,
  future.seed = 12345,
  future.packages = c("ETAS.inlabru", "inlabru", "INLA"),
  future.scheduling = 25
)

future::plan(future::sequential)

rescue_manifest <- do.call(rbind, rescue_rows)
write.csv(rescue_manifest,
          file.path(rescue_dir, "rescue_manifest.csv"),
          row.names = FALSE)

print(table(rescue_manifest$converged, useNA = "ifany"))
print(table(rescue_manifest$truth_kernel,
            rescue_manifest$fitted_kernel,
            rescue_manifest$converged))