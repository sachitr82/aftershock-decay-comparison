#===============================================================================
# Posterior dependence among parameters for correctly specified fits
#===============================================================================

#-------------------------------------------------------------------------------
# Packages and design
#-------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(ETAS.inlabru)

source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Load fit manifest
#-------------------------------------------------------------------------------

fit_manifest <- read.csv(file.path(fit_dir, "fit_manifest.csv"),
                         stringsAsFactors = FALSE)

fit_manifest$fit_file <- file.path(fit_dir,
                                   paste0("truth_", fit_manifest$truth_form),
                                   fit_manifest$file)

correct_manifest <- fit_manifest %>%
  filter(usable, truth_form == fitted_form)

stopifnot(nrow(correct_manifest) == 300)

#-------------------------------------------------------------------------------
# Posterior parameter dependence
#-------------------------------------------------------------------------------

# Number of joint posterior draws
n_corr_samp <- 1000

correlation_pars <- list(
  ou = c("mu", "K", "alpha", "c", "p"),
  mse = c("mu", "K", "alpha", "d", "rho", "gamma"),
  rate_state = c("mu", "K", "alpha", "B", "ta"))

posterior_corr_rows <- vector("list", nrow(correct_manifest))

for (i in seq_len(nrow(correct_manifest))) {
  
  row_i <- correct_manifest[i, ]
  form_i <- row_i$truth_form
  
  message(i, "/", nrow(correct_manifest), " | posterior correlation | ",
          form_i, " | rep ", row_i$rep)
  
  obj_i <- readRDS(row_i$fit_file)
  
  set.seed(5000000 + i)
  
  # n_corr_samp joint posterior draws
  post_i <- ETAS.inlabru::post_sampling(
    input.list = list(model.fit = obj_i$fit,
                      link.functions = obj_i$link.functions,
                      form = form_i),
    n.samp = n_corr_samp,
    max.batch = n_corr_samp)
  
  pars_i <- post_i %>%
    select(all_of(correlation_pars[[form_i]]))
  
  if (form_i == "rate_state") {
    pars_i <- pars_i %>%
      mutate(`1-B` = 1 - B) %>%
      select(mu, K, alpha, `1-B`, ta)
  }
  
  # Compute within-catalogue Spearman correlations across joint posterior draws
  corr_i <- cor(pars_i, method = "spearman", use = "pairwise.complete.obs")
  
  pair_i <- as.data.frame(as.table(corr_i)) %>%
    rename(parameter_1 = Var1, parameter_2 = Var2, correlation = Freq) %>%
    mutate(parameter_1 = as.character(parameter_1),
           parameter_2 = as.character(parameter_2)) %>%
    filter(parameter_1 < parameter_2) %>%
    select(parameter_1, parameter_2, correlation) %>%
    mutate(truth_form = form_i, rep = row_i$rep) %>%
    select(truth_form, rep, parameter_1, parameter_2, correlation)
  
  posterior_corr_rows[[i]] <- pair_i
  
  rm(obj_i, post_i, pars_i, corr_i, pair_i)
  
  if (i %% 25 == 0) gc()
}

posterior_parameter_correlations <- bind_rows(posterior_corr_rows)

stopifnot(
  sum(posterior_parameter_correlations$truth_form == "ou") == 1000,
  sum(posterior_parameter_correlations$truth_form == "mse") == 1500,
  sum(posterior_parameter_correlations$truth_form == "rate_state") == 1000)

#-------------------------------------------------------------------------------
# Median within-catalogue correlation across the 100 repeated catalogues
#-------------------------------------------------------------------------------

posterior_corr_summary <- posterior_parameter_correlations %>%
  group_by(truth_form, parameter_1, parameter_2) %>%
  summarise(n = n(),
            median_correlation = median(correlation),
            .groups = "drop")

#-------------------------------------------------------------------------------
# Temporal-parameter subset
#-------------------------------------------------------------------------------

temporal_pars <- list(ou = c("c", "p"),
                      mse = c("d", "rho", "gamma"),
                      rate_state = c("1-B", "ta"))

posterior_temporal_corr_summary <- posterior_corr_summary %>%
  rowwise() %>%
  filter(parameter_1 %in% temporal_pars[[truth_form]],
         parameter_2 %in% temporal_pars[[truth_form]]) %>%
  ungroup()

#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

posterior_corr_plot_data <- posterior_parameter_correlations %>%
  mutate(form = factor(truth_form, levels = c("ou", "mse", "rate_state"),
                  labels = c("OU", "MSE", "RS")),
         pair = paste(parameter_1, parameter_2, sep = " – "))

p_corr <- ggplot(posterior_corr_plot_data, aes(x = pair, y = correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  facet_wrap(~ form, scales = "free_x") +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = NULL, y = "Posterior Spearman correlation") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))

#-------------------------------------------------------------------------------
# Save
#-------------------------------------------------------------------------------
write.csv(posterior_temporal_corr_summary,
          file.path(recovery_table_dir,
                    "posterior_temporal_correlations_summary.csv"),
          row.names = FALSE)

ggsave(file.path(recovery_figure_dir,
                 "posterior_parameter_correlations.pdf"),
       p_corr, width = 10, height = 4.5)

message("Saved posterior dependence tables to: ", recovery_table_dir)
message("Saved posterior dependence figure to: ", recovery_figure_dir)
message("Finished posterior dependence analysis.")