#===============================================================================
# Shared temporal ETAS fitting helpers
#===============================================================================

# Prepare an event catalogue for temporal ETAS fitting by restricting it to
# the fitting window, ordering events by time and assigning parent indices.
prepare_temporal_catalogue <- function(catalogue, T1, T2) {
  
  catalogue <- catalogue[catalogue$ts >= T1 & catalogue$ts <= T2, ]
  catalogue <- catalogue[order(catalogue$ts), ]
  catalogue$idx.p <- seq_len(nrow(catalogue))
  
  catalogue
}

# Fit one temporal ETAS model using the common prior transformations, numerical
# controls and diagnostics. Optionally request model comparison criteria
fit_temporal_etas <- function(catalogue, fitted_form, binning, fit_control,
                              M0, T1, T2, compute_model_criteria = FALSE,
                              seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Construct common parameter transformations and inlabru options
  link_i <- make_links_P0(fitted_form)
  
  bru_i <- make_bru_options_P0(fitted_form, 
                               rel_tol = fit_control$rel_tol,
                               max_iter = fit_control$max_iter)
  
  # Retain INLA configuration for posterior sampling
  bru_i$control.compute <- list(config = TRUE)
  
  # Request model-comparison criteria only when required
  if (compute_model_criteria) {
    bru_i$control.compute <- list(config = TRUE, dic = TRUE, mlik = TRUE)
  }
  
  start_i <- Sys.time()
  
  # Fit temporal ETAS model
  fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue,
                                       M0 = M0,
                                       T1 = T1,
                                       T2 = T2,
                                       link.functions = link_i,
                                       coef.t. = binning$coef.t,
                                       delta.t. = binning$delta.t,
                                       N.max. = binning$N.max,
                                       bru.opt = bru_i,
                                       form = fitted_form)
  
  runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
  
  # Apply common numerical-quality diagnostics
  diag_i <- fit_diagnostics(fit_i)
  
  list(fit = fit_i,
       link.functions = link_i,
       runtime_minutes = runtime_i,
       diagnostics = diag_i,
       usable = fit_is_usable(diag_i))
}