#===============================================================================
# Ridgecrest posterior parameter summaries
#===============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(ETAS.inlabru)

#-------------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------------

candidate_kernels <- c("ou", "mse", "rate_state")

ridgecrest_fit_dir <- here("outputs", "ridgecrest", "fits", "baseline")

ridgecrest_summary_dir <- here("outputs", "ridgecrest", "summaries")

ridgecrest_figure_dir <- here("outputs", "ridgecrest", "figures")

dir.create(ridgecrest_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ridgecrest_figure_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Plot labels
#-------------------------------------------------------------------------------

kernel_labels <- c(ou = "OU", mse = "MSE", rate_state = "RS")

parameter_labels <- c(mu = "mu", K = "K", alpha = "alpha", c = "c",
                      p = "p", d = "d", rho = "rho", gamma = "gamma", 
                      `1-B` = "1-B", ta = "t[a]")

#-------------------------------------------------------------------------------
# Load baseline fits
#-------------------------------------------------------------------------------

fit_objects <- setNames(
  lapply(candidate_kernels, function(k) {
    readRDS(file.path(ridgecrest_fit_dir, paste0("fit_", k, ".rds")))
  }),
  candidate_kernels)

#-------------------------------------------------------------------------------
# Extract marginal posterior distributions and summaries
#-------------------------------------------------------------------------------

posterior_rows <- vector("list", length(candidate_kernels))
marginal_rows <- vector("list", length(candidate_kernels))

for (i in seq_along(candidate_kernels)) {
  
  kernel_i <- candidate_kernels[i]
  obj_i <- fit_objects[[kernel_i]]
  
  message("Posterior marginals | ", kernel_i)
  
  post_out_i <- ETAS.inlabru::get_posterior_param(
    list(model.fit = obj_i$fit, link.functions = obj_i$link.functions,
         kernel = kernel_i))
  
  # Numerical summaries
  post_i <- post_out_i$post.summary %>%
    rename(parameter = param, q025 = q0.025, q975 = q0.975) %>%
    mutate(model = kernel_i)
  
  # Marginal density curves
  marg_i <- post_out_i$post.df %>%
    rename(parameter = param) %>%
    mutate(model = kernel_i)
  
  posterior_rows[[i]] <- post_i %>%
    select(model, parameter, mean, q025, median, q975)
  
  marginal_rows[[i]] <- marg_i %>%
    select(model, parameter, x, y)
}

posterior_summary <- bind_rows(posterior_rows)
posterior_marginals <- bind_rows(marginal_rows)

#-------------------------------------------------------------------------------
# Express RS parameter as 1-B
#-------------------------------------------------------------------------------

posterior_summary <- posterior_summary %>%
  mutate(mean_old = mean,
         q025_old = q025,
         median_old = median,
         q975_old = q975,
         parameter = ifelse(model == "rate_state" & parameter == "B",
                            "1-B", parameter),
         mean = ifelse(parameter == "1-B", 1 - mean_old, mean_old),
         q025 = ifelse(parameter == "1-B", 1 - q975_old, q025_old),
         median = ifelse(parameter == "1-B", 1 - median_old, median_old),
         q975 = ifelse(parameter == "1-B", 1 - q025_old, q975_old)) %>%
  select(model, parameter, mean, q025, median, q975)

posterior_marginals <- posterior_marginals %>%
  mutate(is_rs_B = model == "rate_state" & parameter == "B",
         x = ifelse(is_rs_B, 1 - x, x),
         parameter = ifelse(is_rs_B, "1-B", parameter)) %>%
  select(-is_rs_B) %>%
  arrange(model, parameter, x)

#-------------------------------------------------------------------------------
# Labels and ordering
#-------------------------------------------------------------------------------

posterior_summary <- posterior_summary %>%
  mutate(model = factor(model, 
                        levels = candidate_kernels, 
                        labels = c("OU", "MSE", "RS")))

posterior_marginals <- posterior_marginals %>%
  mutate(model = factor(model,
                        levels = candidate_kernels,
                        labels = c("OU", "MSE", "RS")))

#-------------------------------------------------------------------------------
# Numerical posterior summaries
#-------------------------------------------------------------------------------

print(as_tibble(posterior_summary), n = Inf)

write.csv(posterior_summary, file.path(ridgecrest_summary_dir,
                                       "baseline_posterior_parameter_summary.csv"),
  row.names = FALSE)

shared_summary <- posterior_summary %>%
  filter(parameter %in% c("mu", "K", "alpha"))

print(shared_summary)

##-------------------------------------------------------------------------------
# Shared ETAS parameters
#-------------------------------------------------------------------------------

shared_marginals <- posterior_marginals %>%
  filter(parameter %in% c("mu", "K", "alpha"))

p_shared <- shared_marginals %>%
  ggplot(aes(x = x, y = y, colour = model, linetype = model)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ parameter, scales = "free",
             labeller = as_labeller(parameter_labels, label_parsed)) +
  labs(x = "Parameter value", y = "Posterior density",
       colour = "Temporal form", linetype = "Temporal form") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

print(p_shared)

ggsave(file.path(ridgecrest_figure_dir,
                 "posterior_shared_parameters.pdf"),
       p_shared, width = 9, height = 3.5)

#-------------------------------------------------------------------------------
# Model-specific temporal parameters
#-------------------------------------------------------------------------------

temporal_marginals <- posterior_marginals %>%
  filter(!parameter %in% c("mu", "K", "alpha"))

p_temporal <- temporal_marginals %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ model + parameter, scales = "free",
             labeller = labeller(
               parameter = as_labeller(parameter_labels, label_parsed))) +
  labs(x = "Parameter value", y = "Posterior density") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

print(p_temporal)

ggsave(file.path(ridgecrest_figure_dir,
                 "posterior_temporal_parameters.pdf"),
       p_temporal, width = 10, height = 5.5)