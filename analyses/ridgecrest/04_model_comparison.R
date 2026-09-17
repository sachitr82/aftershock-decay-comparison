#===============================================================================
# Validate and compare Ridgecrest baseline fits
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(here)

source(here("analyses", "ridgecrest", "00_design.R"))
candidate_forms <- c("ou", "mse", "rate_state")

#-------------------------------------------------------------------------------
# Load baseline fits
#-------------------------------------------------------------------------------

fit_objects <- setNames(
  lapply(candidate_forms, function(k) {
    readRDS(file.path(baseline_fit_dir, paste0("fit_", k, ".rds")))
  }),
  candidate_forms)

#-------------------------------------------------------------------------------
# Validate fits
#-------------------------------------------------------------------------------

validation_rows <- vector("list", length(candidate_forms))

for (i in seq_along(candidate_forms)) {
  
  k <- candidate_forms[i]
  obj_i <- fit_objects[[k]]
  diag_i <- obj_i$diagnostics
  lml_i <- obj_i$lml
  
  validation_rows[[i]] <- data.frame(
    fitted_form = k,
    converged = diag_i$converged,
    hit_max = diag_i$hit_max,
    inla_failure = diag_i$inla_failure,
    nan_inf_logl = diag_i$nan_inf_logl,
    vb_aborted = diag_i$vb_aborted,
    n_iter = diag_i$n_iter,
    min_sd = diag_i$min_sd,
    degenerate = diag_i$degenerate,
    usable = obj_i$usable,
    lml = lml_i,
    lml_valid = length(lml_i) == 1 && is.finite(lml_i))
}

fit_validation <- bind_rows(validation_rows) %>%
  mutate(usable = usable & lml_valid)

write.csv(fit_validation,
          file.path(model_comparison_table_dir, "baseline_fit_validation.csv"),
          row.names = FALSE)

# Do not compare models if any baseline fit is unusable
stopifnot(all(fit_validation$usable))

#-------------------------------------------------------------------------------
# Marginal-likelihood comparison
#-------------------------------------------------------------------------------

lml_summary <- fit_validation %>%
  select(fitted_form, lml) %>%
  arrange(desc(lml))

write.csv(lml_summary,
          file.path(model_comparison_table_dir, "baseline_lml.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Pairwise Bayes factors
#-------------------------------------------------------------------------------

# Bayes factor as exponential LML difference
pairs <- combn(candidate_forms, 2, simplify = FALSE)

comparison_rows <- lapply(pairs, function(pair_i) {
  
  model_a <- pair_i[1]
  model_b <- pair_i[2]
  
  lml_a <- lml_summary$lml[lml_summary$fitted_form == model_a]
  
  lml_b <- lml_summary$lml[lml_summary$fitted_form == model_b]
  
  delta_lml <- lml_a - lml_b
  
  data.frame(model_a = model_a,
             model_b = model_b, 
             lml_a = lml_a,
             lml_b = lml_b, 
             delta_lml = delta_lml,
             bf_a_vs_b = exp(delta_lml),
             preferred_model = ifelse(delta_lml > 0, model_a, model_b))
})

pairwise_comparison <- bind_rows(comparison_rows)

write.csv(pairwise_comparison, file.path(model_comparison_table_dir, 
                                         "baseline_pairwise_lml.csv"),
          row.names = FALSE)

message("Saved baseline model-comparison tables to: ",
        model_comparison_table_dir)
message("Finished Ridgecrest baseline fit comparison.")