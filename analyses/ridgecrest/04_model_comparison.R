#===============================================================================
# Validate and compare Ridgecrest baseline fits
#===============================================================================

library(dplyr)
library(here)

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

ridgecrest_fit_dir <- here("outputs", "ridgecrest", "baseline")

ridgecrest_summary_dir <- here(ridgecrest_fit_dir, "summaries")

dir.create(ridgecrest_summary_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load baseline fits
#-------------------------------------------------------------------------------

fit_objects <- setNames(
  lapply(candidate_kernels, function(k) {
    readRDS(file.path(ridgecrest_fit_dir, paste0("fit_", k, ".rds")))
  }),
  candidate_kernels)

#-------------------------------------------------------------------------------
# Validate fits
#-------------------------------------------------------------------------------

validation_rows <- vector("list", length(candidate_kernels))

for (i in seq_along(candidate_kernels)) {
  
  k <- candidate_kernels[i]
  obj_i <- fit_objects[[k]]
  fixed_i <- obj_i$fit$summary.fixed
  
  n_degenerate_i <- sum(!is.finite(fixed_i$sd) | fixed_i$sd < 1e-6)
  
  lml_i <- obj_i$lml
  
  validation_rows[[i]] <- data.frame(
    fitted_kernel = k,
    converged = isTRUE(obj_i$converged),
    hit_max = isTRUE(obj_i$hit_max),
    inla_failure = isTRUE(obj_i$inla_failure),
    nan_inf_logl = isTRUE(obj_i$nan_inf_logl),
    vb_aborted = isTRUE(obj_i$vb_aborted),
    n_iter = obj_i$n_iter,
    n_degenerate = n_degenerate_i,
    posterior_valid = n_degenerate_i == 0,
    lml = lml_i,
    lml_valid = length(lml_i) == 1 && is.finite(lml_i))
}

fit_validation <- bind_rows(validation_rows) %>%
  mutate(usable = converged &
                  !hit_max &
                  !inla_failure &
                  !nan_inf_logl &
                  !vb_aborted &
                  posterior_valid &
                  lml_valid)

print(fit_validation)

write.csv(fit_validation,
          file.path(ridgecrest_summary_dir, "baseline_fit_validation.csv"),
          row.names = FALSE)

# Do not compare models if any baseline fit is unusable
stopifnot(all(fit_validation$usable))

#-------------------------------------------------------------------------------
# Marginal-likelihood comparison
#-------------------------------------------------------------------------------

lml_summary <- fit_validation %>%
  select(fitted_kernel, lml) %>%
  arrange(desc(lml))

print(lml_summary)

write.csv(lml_summary,
          file.path(ridgecrest_summary_dir, "baseline_lml.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Pairwise Bayes factors
#-------------------------------------------------------------------------------

pairs <- combn(candidate_kernels, 2, simplify = FALSE)

comparison_rows <- lapply(pairs, function(pair_i) {
  
  model_a <- pair_i[1]
  model_b <- pair_i[2]
  
  lml_a <- lml_summary$lml[lml_summary$fitted_kernel == model_a]
  
  lml_b <- lml_summary$lml[lml_summary$fitted_kernel == model_b]
  
  delta_lml <- lml_a - lml_b
  
  data.frame(model_a = model_a, model_b = model_b, 
             lml_a = lml_a, lml_b = lml_b, 
             delta_lml = delta_lml,
             bf_a_vs_b = exp(delta_lml),
             preferred_model = ifelse(delta_lml > 0, model_a, model_b))
})

pairwise_comparison <- bind_rows(comparison_rows)

print(pairwise_comparison)

write.csv(pairwise_comparison, file.path(ridgecrest_summary_dir, 
                                         "baseline_pairwise_lml.csv"),
          row.names = FALSE)