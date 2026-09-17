#===============================================================================
# Shared fit diagnostics
#===============================================================================

# Extract common numerical diagnostics from an ETAS.inlabru fit.
# Used by pilot, simulation and Ridgecrest fitting workflows.

fit_diagnostics <- function(fit, sd_tol = 1e-6) {
  
  log_i <- as.character(inlabru::bru_log(fit))
  
  # Detect known INLA failure messages
  inla_failure <- any(grepl(
    paste("Problem in inla",
          "Giving up and returning last successfully obtained result",
          "inla-program exited with an error",
          "maximum number of tries has been reached",
          "Newton-Raphson optimizer did not converge",
          sep = "|"), log_i))
  
  # Detect additional numerical warnings
  nan_inf_logl <- any(grepl("NAN/INF values in logl", log_i, fixed = TRUE))
  vb_aborted <- any(grepl("max_correction|vb.correction.*aborted", log_i))
  
  # Classify outer inlabru convergence status
  converged <- !inla_failure &&
               any(grepl("Convergence criterion met", log_i, fixed = TRUE))
  hit_max <- !inla_failure &&
             any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
  
  fit_status <- dplyr::case_when(inla_failure ~ "inla_failure",
                                 converged ~ "converged",
                                 hit_max ~ "hit_max",
                                 TRUE ~ "unknown")
  
  n_iter <- max(fit$bru_iinla$track$iteration, na.rm = TRUE)
  
  # Flag near-degenerate posterior marginals
  min_sd <- min(fit$summary.fixed$sd, na.rm = TRUE)
  degenerate <- isTRUE(min_sd < sd_tol)
  
  list(fit_status = fit_status,
       converged = converged,
       hit_max = hit_max,
       inla_failure = inla_failure,
       nan_inf_logl = nan_inf_logl,
       vb_aborted = vb_aborted,
       n_iter = n_iter,
       min_sd = min_sd,
       degenerate = degenerate)
}

# Determine whether a fitted model passes the common numerical-quality checks.
fit_is_usable <- function(diag) {
  
  isTRUE(diag$converged) &&
    !isTRUE(diag$hit_max) &&
    !isTRUE(diag$degenerate) &&
    !isTRUE(diag$vb_aborted) &&
    !isTRUE(diag$inla_failure) &&
    !isTRUE(diag$nan_inf_logl)
}