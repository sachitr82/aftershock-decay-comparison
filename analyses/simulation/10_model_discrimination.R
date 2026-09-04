#===============================================================================
# Model discrimination under known truth
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(patchwork)

source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------------------

main_fit_dir <- file.path(fit_dir, "main_simulation")
comparison_dir <- file.path(main_fit_dir, "model_discrimination")
figure_dir <- file.path(comparison_dir, "figures")

dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

candidate_kernels <- c("ou", "mse", "rate_state")

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(main_fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)

fit_manifest$fit_file <- file.path(
  main_fit_dir, paste0("truth_", fit_manifest$truth_kernel), fit_manifest$file)

stopifnot(nrow(fit_manifest) == 900)

fit_manifest <- fit_manifest %>%
  mutate(usable = converged & !hit_max & !degenerate & !vb_aborted & 
           !inla_failure & !nan_inf_logl)

message("Usable fits: ", sum(fit_manifest$usable, na.rm = TRUE), 
        "/", nrow(fit_manifest))

#-------------------------------------------------------------------------------
# Extract INLA model-comparison quantities (full DIC breakdown)
#-------------------------------------------------------------------------------

extract_mlik <- function(fit, type = c("integration", "Gaussian")) {
  
  type <- match.arg(type)
  
  if (is.null(fit$mlik)) return(NA_real_)
  
  rn <- rownames(fit$mlik)
  if (is.null(rn)) return(NA_real_)
  
  idx <- grep(type, rn, ignore.case = TRUE)
  
  if (length(idx) != 1)  return(NA_real_)
  
  as.numeric(fit$mlik[idx, 1])
}


usable_manifest <- fit_manifest %>% filter(usable)

score_rows <- vector("list", nrow(usable_manifest))

for (i in seq_len(nrow(usable_manifest))) {
  
  row_i <- usable_manifest[i, ]
  
  message(i, "/", nrow(usable_manifest), " | truth ", row_i$truth_kernel, 
          " | fit ", row_i$fitted_kernel, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  fit_i <- obj_i$fit
  
  score_rows[[i]] <- data.frame(truth_kernel = row_i$truth_kernel,
                                fitted_kernel = row_i$fitted_kernel,
                                rep = row_i$rep,
    
    dic = if (!is.null(fit_i$dic$dic)) {
      as.numeric(fit_i$dic$dic)
    } else {
      NA_real_
    },
    mean_deviance = if (!is.null(fit_i$dic$mean.deviance)) {
      as.numeric(fit_i$dic$mean.deviance)
    } else {
      NA_real_
    },

    p_eff = if (!is.null(fit_i$dic$p.eff)) {
      as.numeric(fit_i$dic$p.eff)
    } else {
      NA_real_
    },
    
    lml_integration = extract_mlik(fit_i, "integration"),
    lml_gaussian = extract_mlik(fit_i, "Gaussian")
  )
}

model_scores_full <- bind_rows(score_rows)

#-------------------------------------------------------------------------------
# Retain criteria used for model discrimination (both LMLs agree, retain one)
#-------------------------------------------------------------------------------

model_scores <- model_scores_full %>%
  select(truth_kernel, fitted_kernel, rep, dic, lml = lml_integration)

#-------------------------------------------------------------------------------
# Make all criteria share orientation (higher = better, DIC becomes -DIC)
#-------------------------------------------------------------------------------

score_long <- model_scores %>%
  pivot_longer(cols = c(dic, lml),
               names_to = "criterion", values_to = "score") %>%
  mutate(oriented_score = if_else(criterion == "dic", -score, score))

#-------------------------------------------------------------------------------
# Catalogue-level selection function
#-------------------------------------------------------------------------------

make_decisions <- function(dat, candidates, truths = candidates) {
  
  # Retain only candidate models
  comparison_data <- dat %>%
    filter(truth_kernel %in% truths, fitted_kernel %in% candidates)
  
  # Remove catalogues with infinite score and which don't have all required fits
  # to compare
  eligible <- comparison_data %>%
    group_by(truth_kernel, rep, criterion) %>%
    summarise(n_models = n_distinct(fitted_kernel),
      all_finite = all(is.finite(oriented_score)),
      .groups = "drop") %>%
    filter(n_models == length(candidates), all_finite)
  
  comparison_data <- comparison_data %>%
    inner_join(eligible %>% select(truth_kernel, rep, criterion),
               by = c("truth_kernel", "rep", "criterion"))
  
  keys <- comparison_data %>% distinct(truth_kernel, rep, criterion)
  
  decision_rows <- vector("list", nrow(keys))
  
  for (j in seq_len(nrow(keys))) {
    
    key_j <- keys[j, ]
    
    d <- comparison_data %>%
      filter(truth_kernel == key_j$truth_kernel,
             rep == key_j$rep, criterion == key_j$criterion)
    
    best_score <- max(d$oriented_score)
    
    stopifnot(sum(d$oriented_score == best_score) == 1)
    
    best_idx <- which.max(d$oriented_score)
    
    selected_kernel <- d$fitted_kernel[best_idx]
    
    truth_score <- d$oriented_score[d$fitted_kernel == key_j$truth_kernel]
    
    wrong_scores <- d$oriented_score[d$fitted_kernel != key_j$truth_kernel]
    
    truth_margin <- truth_score - max(wrong_scores)
    
    decision_rows[[j]] <- data.frame(
      truth_kernel = key_j$truth_kernel,
      rep = key_j$rep,
      criterion = key_j$criterion,
      selected_kernel = selected_kernel,
      truth_margin = truth_margin,
      correct = selected_kernel == key_j$truth_kernel
    )
  }
  
  bind_rows(decision_rows)
}

#-------------------------------------------------------------------------------
# Three-way discrimination: OU vs MSE vs RS
#-------------------------------------------------------------------------------

three_way <- make_decisions(score_long, candidates = candidate_kernels, 
                            truths = candidate_kernels)

# Check eligible catalogue counts
three_way_n <- three_way %>%
  count(criterion, truth_kernel, name = "n_eligible")

print(three_way_n)

#-------------------------------------------------------------------------------
# Three-way selection proportions
#-------------------------------------------------------------------------------

three_way_selection <- three_way %>%
  count(criterion, truth_kernel, selected_kernel, name = "n_selected") %>%
  group_by(criterion, truth_kernel) %>%
  complete(selected_kernel = candidate_kernels,
           fill = list(n_selected = 0)) %>%
  mutate(n_eligible = sum(n_selected),
         proportion = n_selected / n_eligible) %>%
  ungroup()

three_way_selection_compact <- three_way_selection %>%
  filter(n_selected > 0) %>%
  select(criterion, truth_kernel, selected_kernel,
         n_selected, n_eligible, proportion)

write.csv(three_way_selection_compact,
          file.path(comparison_dir, "three_way_selection.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Pair-specific criterion differences - how strongly each is favoured
#-------------------------------------------------------------------------------

# Extract criterion values from the correctly specified fit for each catalogue
truth_scores <- model_scores %>%
  filter(fitted_kernel == truth_kernel) %>%
  select(truth_kernel, rep, dic, lml) %>%
  rename(dic_truth = dic, lml_truth = lml)

# Match each misspecified fit to the correct fit for the same catalogue
pairwise_differences <- model_scores %>%
  filter(fitted_kernel != truth_kernel) %>%
  select(truth_kernel, fitted_kernel, rep, dic, lml) %>%
  rename(competitor_kernel = fitted_kernel,
         dic_competitor = dic,
         lml_competitor = lml) %>%
  inner_join(truth_scores, by = c("truth_kernel", "rep")) %>%
  mutate(delta_dic = dic_competitor - dic_truth,
         delta_lml = lml_truth - lml_competitor)

#-------------------------------------------------------------------------------
# Long format
#-------------------------------------------------------------------------------

pairwise_differences_long <- pairwise_differences %>%
  select(truth_kernel, competitor_kernel, rep, delta_dic, delta_lml) %>%
  pivot_longer(cols = c(delta_dic, delta_lml),
               names_to = "criterion", values_to = "delta") %>%
  mutate(criterion = recode(criterion, delta_dic = "DIC", delta_lml = "LML"),
         truth_label = kernel_labels[truth_kernel],
         competitor_label = kernel_labels[competitor_kernel],
         comparison = paste0(truth_label, " truth vs ", competitor_label)) %>%
  filter(is.finite(delta))
#-------------------------------------------------------------------------------
# Pair-specific summaries
#-------------------------------------------------------------------------------

pairwise_difference_summary <- pairwise_differences_long %>%
  group_by(criterion, truth_kernel, competitor_kernel, comparison) %>%
  summarise(n = n(),
            n_truth_favoured = sum(delta > 0),
            prop_truth_favoured = mean(delta > 0),
            median_delta = median(delta),
            q10_delta = quantile(delta, 0.10),
            q90_delta = quantile(delta, 0.90),
            .groups = "drop")

print(pairwise_difference_summary)

write.csv(pairwise_difference_summary,
          file.path(comparison_dir, "pairwise_criterion_difference_summary.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Order comparisons consistently in figures
#-------------------------------------------------------------------------------

comparison_levels <- c("OU truth vs MSE",
                       "OU truth vs RS",
                       "MSE truth vs OU",
                       "MSE truth vs RS",
                       "RS truth vs OU",
                       "RS truth vs MSE")

#-------------------------------------------------------------------------------
# DIC difference plot
#-------------------------------------------------------------------------------

dic_plot_data <- pairwise_differences_long %>%
  filter(criterion == "DIC") %>%
  group_by(comparison) %>%
  mutate(n_pair = n(),
         comparison_n = paste0(comparison, " (n = ", n_pair, ")")) %>%
  ungroup()

dic_levels <- dic_plot_data %>%
  distinct(comparison, comparison_n) %>%
  mutate(comparison = factor(comparison, levels = comparison_levels)) %>%
  arrange(comparison) %>%
  pull(comparison_n)

dic_plot_data <- dic_plot_data %>%
  mutate(comparison_n = factor(comparison_n, levels = dic_levels))

p_dic <- ggplot(dic_plot_data, aes(x = delta, y = comparison_n)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_jitter(position = position_jitter(height = 0.15, width = 0, seed = 123),
              alpha = 0.35, size = 1.3) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.25) +
  labs(x = expression(Delta*"DIC = DIC"[alt] - "DIC"[true]),
       y = NULL) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.margin = margin(t = 0, r = 30, b = 0, l = 0))

#-------------------------------------------------------------------------------
# Log marginal-likelihood difference plot
#-------------------------------------------------------------------------------

lml_plot_data <- pairwise_differences_long %>%
  filter(criterion == "LML") %>%
  group_by(comparison) %>%
  mutate(n_pair = n(),
         comparison_n = paste0(comparison, " (n = ", n_pair, ")")) %>%
  ungroup()

lml_levels <- lml_plot_data %>%
  distinct(comparison, comparison_n) %>%
  mutate(comparison = factor(comparison, levels = comparison_levels)) %>%
  arrange(comparison) %>%
  pull(comparison_n)

lml_plot_data <- lml_plot_data %>%
  mutate(comparison_n = factor(comparison_n, levels = lml_levels))

p_lml <- ggplot(lml_plot_data, aes(x = delta, y = comparison_n)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_jitter(position = position_jitter(height = 0.15, width = 0, seed = 123),
              alpha = 0.35, size = 1.3) +
  geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.25) +
  labs(x = expression(Delta*"LML = LML"[true] - "LML"[alt]),
       y = NULL) +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 30))

p_comparison <- (p_dic | p_lml) +
  plot_layout(axes = "collect_y") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

print(p_comparison)

ggsave(file.path(figure_dir, "pairwise_model_comparison.pdf"), p_comparison,
       width = 12, height = 4.5)

#-------------------------------------------------------------------------------
# Identify marginal-likelihood misclassifications
#-------------------------------------------------------------------------------

lml_errors <- pairwise_differences %>%
  filter(is.finite(delta_lml), delta_lml < 0) %>%
  mutate(bf_truth_vs_alt = exp(delta_lml), bf_alt_vs_truth = exp(-delta_lml)) %>%
  select(truth_kernel, competitor_kernel, rep, delta_lml,
         bf_truth_vs_alt, bf_alt_vs_truth)

print(lml_errors)
#-------------------------------------------------------------------------------
# DIC diagnostics - unclear
#-------------------------------------------------------------------------------

dic_component_summary <- model_scores_full %>%
  filter(fitted_kernel != truth_kernel) %>%
  select(truth_kernel, fitted_kernel, rep, mean_deviance, p_eff) %>%
  rename(competitor_kernel = fitted_kernel,
         mean_deviance_competitor = mean_deviance,
         p_eff_competitor = p_eff) %>%
  inner_join(truth_dic, by = c("truth_kernel", "rep")) %>%
  mutate(delta_mean_deviance = mean_deviance_competitor - mean_deviance_truth,
         delta_p_eff = p_eff_competitor - p_eff_truth) %>%
  group_by(truth_kernel, competitor_kernel) %>%
  summarise(
    median_delta_mean_deviance = median(delta_mean_deviance),
    median_delta_p_eff = median(delta_p_eff),
    prop_mean_deviance_favours_truth = mean(delta_mean_deviance > 0),
    prop_p_eff_favours_truth = mean(delta_p_eff > 0),
    .groups = "drop")

print(dic_component_summary)
