#===============================================================================
# Ridgecrest magnitude-threshold sensitivity
#===============================================================================

library(dplyr)
library(ETAS.inlabru)
library(here)
library(inlabru)

source(here("analyses", "simulation", "00_design.R"))
source(here("src", "eda-helpers", "mc_helpers.R"))

#-------------------------------------------------------------------------------
# Load and prepare pre-threshold Ridgecrest catalogue
#-------------------------------------------------------------------------------

cat_rc <- readRDS(here("data", "derived", "ridgecrest_temporal_catalogue.rds"))

# Mainshock defines t = 0
mainshock_i <- which.max(cat_rc$mag)
mainshock_datetime <- cat_rc$datetime[mainshock_i]

# Convert to form expected by Temporal.ETAS()
cat_rc <- cat_rc %>%
  mutate(ts = as.numeric(difftime(datetime, mainshock_datetime, units = "days")),
         magnitudes = mag) %>%
  arrange(ts)

cat_rc$idx.p <- seq_len(nrow(cat_rc))

#-------------------------------------------------------------------------------
# Ridgecrest observation window
#-------------------------------------------------------------------------------

window_start <- as.POSIXct("2016-01-01 00:00:00", tz = "UTC")
window_end <- as.POSIXct("2025-12-31 23:59:59", tz = "UTC")

T_fit_start <- as.numeric(
  difftime(window_start, mainshock_datetime, units = "days"))

T_fit_end <- as.numeric(
  difftime(window_end, mainshock_datetime, units = "days"))

#-------------------------------------------------------------------------------
# Fixed fitting design
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

# Baseline M0 = 2.5 already fitted in 01_fit_baseline.R
M0_extra <- c(3.0, 3.3)
M0_all <- c(2.5, 3.0, 3.3)

stopifnot(temporal_binning$N.max == 14)
stopifnot(fit_control$rel_tol == 0.1)

sensitivity_dir <- here(
  "outputs", "ridgecrest", "sensitivity", "magnitude_threshold")

ridgecrest_fit_dir <- file.path(sensitivity_dir, "fits")
ridgecrest_summary_dir <- file.path(sensitivity_dir, "summaries")

dir.create(ridgecrest_fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ridgecrest_summary_dir, recursive = TRUE, showWarnings = FALSE)

# Separate seeds by threshold and fitted model
fit_seeds <- matrix(
  c(410001, 410002, 410003, 420001, 420002, 420003),
  nrow = length(M0_extra),
  byrow = TRUE,
  dimnames = list(as.character(M0_extra),candidate_kernels))

#-------------------------------------------------------------------------------
# Fit diagnostics
#-------------------------------------------------------------------------------

fit_diagnostics <- function(fit) {
  log_i <- as.character(inlabru::bru_log(fit))
  
  inla_failure <- any(grepl(
    "Problem in inla|Giving up and returning last successfully obtained result|inla-program exited with an error|maximum number of tries has been reached|Newton-Raphson optimizer did not converge",
    log_i))
  
  nan_inf_logl <- any(grepl("NAN/INF values in logl", log_i, fixed = TRUE))
  vb_aborted <- any(grepl("max_correction|vb.correction.*aborted", log_i))
  
  converged <- !inla_failure &&
    any(grepl("Convergence criterion met", log_i, fixed = TRUE))
  
  hit_max <- !inla_failure &&
    any(grepl("Maximum iterations reached", log_i, fixed = TRUE))
  
  fit_status <- case_when(inla_failure ~ "inla_failure",
                          converged ~ "converged",
                          hit_max ~ "hit_max",
                          TRUE ~ "unknown")
  
  n_iter <- max(fit$bru_iinla$track$iteration, na.rm = TRUE)
  
  list(fit_status = fit_status,
       converged = converged,
       hit_max = hit_max,
       inla_failure = inla_failure,
       nan_inf_logl = nan_inf_logl,
       vb_aborted = vb_aborted,
       n_iter = n_iter)
}

#-------------------------------------------------------------------------------
# Posterior diagnostics
#-------------------------------------------------------------------------------

posterior_diagnostics <- function(fit) {
  x <- fit$summary.fixed
  
  data.frame(latent_parameter = rownames(x),
             mean = x$mean,
             sd = x$sd,
             q025 = x$`0.025quant`,
             median = x$`0.5quant`,
             q975 = x$`0.975quant`,
             degenerate = !is.finite(x$sd) | x$sd < 1e-6,
             row.names = NULL)
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
# Catalogue summaries at each magnitude threshold
#-------------------------------------------------------------------------------

catalogue_summary <- vector("list", length(M0_all))

for (i in seq_along(M0_all)) {
  
  M0_i <- M0_all[i]
  
  catalogue_i <- cat_rc %>%
    filter(magnitudes >= M0_i) %>%
    arrange(ts)
  
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  stopifnot(is.data.frame(catalogue_i))
  stopifnot(all(c("ts", "magnitudes", "idx.p") %in% names(catalogue_i)))
  stopifnot(all(is.finite(catalogue_i$ts)))
  stopifnot(all(is.finite(catalogue_i$magnitudes)))
  stopifnot(all(catalogue_i$magnitudes >= M0_i))
  stopifnot(all(diff(catalogue_i$ts) >= 0))
  stopifnot(any(catalogue_i$ts == 0))
  stopifnot(min(catalogue_i$ts) >= T_fit_start)
  stopifnot(max(catalogue_i$ts) <= T_fit_end)
  
  b_i <- b_aki(catalogue_i$magnitudes, mc = M0_i, dm = 0.1)
  b_hat_i <- b_i[["b"]]
  b_err_i <- b_i[["b_err"]]
  beta_i <- b_hat_i * log(10)
  
  catalogue_summary[[i]] <- data.frame(M0 = M0_i,
                                       n_fit = nrow(catalogue_i),
                                       b_hat = b_hat_i,
                                       b_err = b_err_i,
                                       beta = beta_i)
}

catalogue_summary <- bind_rows(catalogue_summary)

print(catalogue_summary)

write.csv(catalogue_summary,
          file.path(ridgecrest_summary_dir, "catalogue_summary.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Fit OU, MSE and rate-state at alternative thresholds
#-------------------------------------------------------------------------------

fit_manifest <- vector("list", length(M0_extra) * length(candidate_kernels))

row_id <- 1

for (M0_i in M0_extra) {
  
  catalogue_i <- cat_rc %>%
    filter(magnitudes >= M0_i) %>%
    arrange(ts)
  
  catalogue_i$idx.p <- seq_len(nrow(catalogue_i))
  
  summary_i <- catalogue_summary %>%
    filter(M0 == M0_i)
  
  b_hat_i <- summary_i$b_hat
  beta_i <- summary_i$beta
  
  cat("============================================================\n")
  cat("M0 =", M0_i,
      "| n =", nrow(catalogue_i),
      "| b =", round(b_hat_i, 4),
      "| beta =", round(beta_i, 4), "\n")
  cat("============================================================\n")
  
  for (k in seq_along(candidate_kernels)) {
    
    fitted_i <- candidate_kernels[k]
    
    M0_tag <- gsub("\\.", "p", as.character(M0_i))
    
    outfile_i <- file.path(ridgecrest_fit_dir,
                           paste0("fit_M0_", M0_tag, "_", fitted_i, ".rds"))
    
    message("============================================================")
    message("Ridgecrest | M0 = ", M0_i, " | Fit: ", fitted_i)
    message("============================================================")
    
    if (file.exists(outfile_i)) {
      
      message("Already exists: ", basename(outfile_i))
      output_i <- readRDS(outfile_i)
      
    } else {
      
      set.seed(fit_seeds[as.character(M0_i), fitted_i])
      
      link_i <- make_links_P0(fitted_i)
      
      bru_i <- make_bru_options_P0(fitted_i,
                                   rel_tol = fit_control$rel_tol,
                                   max_iter = 200)
      
      bru_i$control.compute <- list(config = TRUE, mlik = TRUE)
      
      start_i <- Sys.time()
      
      fit_i <- ETAS.inlabru::Temporal.ETAS(total.data = catalogue_i,
                                           M0 = M0_i,
                                           T1 = T_fit_start,
                                           T2 = T_fit_end,
                                           link.functions = link_i,
                                           coef.t. = temporal_binning$coef.t,
                                           delta.t. = temporal_binning$delta.t,
                                           N.max. = temporal_binning$N.max,
                                           bru.opt = bru_i,
                                           kernel = fitted_i)
      
      runtime_i <- as.numeric(difftime(Sys.time(), start_i, units = "mins"))
      
      diag_i <- fit_diagnostics(fit_i)
      posterior_diag_i <- posterior_diagnostics(fit_i)
      lml_i <- extract_lml(fit_i)
      
      output_i <- list(fit = fit_i,
                       link.functions = link_i,
                       fitted_kernel = fitted_i,
                       M0 = M0_i,
                       T1 = T_fit_start,
                       T2 = T_fit_end,
                       n_fit = nrow(catalogue_i),
                       binning_parameters = temporal_binning,
                       rel_tol = fit_control$rel_tol,
                       max_iter = 200,
                       fit_seed = fit_seeds[as.character(M0_i), fitted_i],
                       runtime_minutes = runtime_i,
                       fit_status = diag_i$fit_status,
                       converged = diag_i$converged,
                       hit_max = diag_i$hit_max,
                       inla_failure = diag_i$inla_failure,
                       nan_inf_logl = diag_i$nan_inf_logl,
                       vb_aborted = diag_i$vb_aborted,
                       n_iter = diag_i$n_iter,
                       posterior_diagnostics = posterior_diag_i,
                       lml = lml_i)
      
      tmp_i <- tempfile(pattern = "fit_", tmpdir = dirname(outfile_i),
                        fileext = ".rds")
      
      saveRDS(output_i, tmp_i)
      
      if (!file.rename(tmp_i, outfile_i)) {
        stop("Could not move temporary fit file to final output.")
      }
      
      message("Saved ", basename(outfile_i),
              " | runtime = ", round(runtime_i, 2), " min",
              " | LML = ", round(lml_i, 2))
    }
    
    posterior_diag_i <- posterior_diagnostics(output_i$fit)
    
    n_degenerate_i <- sum(posterior_diag_i$degenerate)
    posterior_valid_i <- n_degenerate_i == 0
    
    if (!posterior_valid_i) {
      warning(fitted_i, " fit contains ",
              n_degenerate_i,
              " near-degenerate latent posterior marginal(s).")
    }
    
    fit_manifest[[row_id]] <- data.frame(M0 = M0_i,
                                         b_hat = b_hat_i,
                                         beta = beta_i,
                                         fitted_kernel = fitted_i,
                                         n_fit = output_i$n_fit,
                                         lml = output_i$lml,
                                         runtime_minutes = output_i$runtime_minutes,
                                         fit_status = output_i$fit_status,
                                         converged = output_i$converged,
                                         posterior_valid = posterior_valid_i,
                                         hit_max = output_i$hit_max,
                                         inla_failure = output_i$inla_failure,
                                         nan_inf_logl = output_i$nan_inf_logl,
                                         vb_aborted = output_i$vb_aborted,
                                         n_iter = output_i$n_iter,
                                         n_degenerate = n_degenerate_i,
                                         file = basename(outfile_i))
    
    row_id <- row_id + 1
  }
}

#-------------------------------------------------------------------------------
# Save alternative-threshold fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- bind_rows(fit_manifest) %>%
  arrange(M0, desc(lml))

print(fit_manifest)

write.csv(fit_manifest, file.path(ridgecrest_fit_dir, "fit_manifest.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Add baseline M0 = 2.5 fits
#-------------------------------------------------------------------------------

baseline_fit_dir <- here("outputs", "ridgecrest", "baseline")

baseline_summary <- catalogue_summary %>%
  filter(M0 == 2.5)

baseline_rows <- vector("list", length(candidate_kernels))

for (k in seq_along(candidate_kernels)) {
  
  fitted_i <- candidate_kernels[k]
  
  output_i <- readRDS(
    file.path(baseline_fit_dir, paste0("fit_", fitted_i, ".rds")))
  
  baseline_rows[[k]] <- data.frame(M0 = 2.5,
                                   b_hat = baseline_summary$b_hat,
                                   beta = baseline_summary$beta,
                                   fitted_kernel = fitted_i,
                                   n_fit = output_i$n_fit,
                                   lml = output_i$lml)
}

baseline_rows <- bind_rows(baseline_rows)

#-------------------------------------------------------------------------------
# Combine marginal likelihoods across thresholds
#-------------------------------------------------------------------------------

model_scores <- bind_rows(fit_manifest %>%
                            select(M0, b_hat, beta, fitted_kernel, n_fit, lml),
                          baseline_rows) %>%
  arrange(M0, fitted_kernel)

print(model_scores)

write.csv(model_scores, file.path(ridgecrest_summary_dir, "model_scores.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Pairwise marginal-likelihood differences within each threshold
#-------------------------------------------------------------------------------

comparison_pairs <- list(c("ou", "mse"), c("ou", "rate_state"), 
                         c("mse", "rate_state"))

comparison_rows <- list()
row_id <- 1

for (M0_i in M0_all) {
  
  scores_i <- model_scores %>%
    filter(M0 == M0_i)
  
  stopifnot(nrow(scores_i) == 3)
  
  for (pair_i in comparison_pairs) {
    
    model_a <- pair_i[1]
    model_b <- pair_i[2]
    
    lml_a <- scores_i$lml[scores_i$fitted_kernel == model_a]
    lml_b <- scores_i$lml[scores_i$fitted_kernel == model_b]
    
    stopifnot(length(lml_a) == 1)
    stopifnot(length(lml_b) == 1)
    
    delta_i <- lml_a - lml_b
    
    comparison_rows[[row_id]] <- data.frame(M0 = M0_i,
                                            n_fit = unique(scores_i$n_fit),
                                            b_hat = unique(scores_i$b_hat),
                                            beta = unique(scores_i$beta),
                                            model_a = model_a,
                                            model_b = model_b,
                                            lml_a = lml_a,
                                            lml_b = lml_b,
                                            delta_lml = delta_i,
                                            preferred_model = 
                                              ifelse(delta_i > 0, model_a, model_b))
    
    row_id <- row_id + 1
  }
}

model_comparison <- bind_rows(comparison_rows)

print(model_comparison)

write.csv(model_comparison, 
          file.path(ridgecrest_summary_dir, "model_comparison.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Compact sensitivity summary
#-------------------------------------------------------------------------------

ou_mse <- model_comparison %>%
  filter(model_a == "ou", model_b == "mse") %>%
  select(M0, n_fit, b_hat, beta, delta_ou_mse = delta_lml)

ou_rs <- model_comparison %>%
  filter(model_a == "ou", model_b == "rate_state") %>%
  select(M0, delta_ou_rs = delta_lml)

mse_rs <- model_comparison %>%
  filter(model_a == "mse", model_b == "rate_state") %>%
  select(M0, delta_mse_rs = delta_lml)

sensitivity_summary <- ou_mse %>%
  left_join(ou_rs, by = "M0") %>%
  left_join(mse_rs, by = "M0") %>%
  arrange(M0)

print(sensitivity_summary)

write.csv(sensitivity_summary,
          file.path(ridgecrest_summary_dir, "magnitude_threshold_summary.csv"),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(sensitivity_dir, "sessionInfo.txt"))

message("Finished Ridgecrest magnitude-threshold sensitivity.")