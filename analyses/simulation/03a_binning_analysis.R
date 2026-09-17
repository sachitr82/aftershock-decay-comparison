#===============================================================================
# Analyse temporal binning check
#===============================================================================

library(ggplot2)
library(dplyr)
library(here)

source(here("analyses", "simulation", "00_design.R"))
#-------------------------------------------------------------------------------
# Load numerical results from 03
#-------------------------------------------------------------------------------

results <- read.csv(file.path(binning_fit_dir, "binning_manifest.csv"))
results$binning <- factor(results$binning, 
                          levels = c("N8", "N10", "N12", "N14", "N16"))
stopifnot(nrow(results) == 15)

#-------------------------------------------------------------------------------
# Changes relative to N8 (iteration and run time)
#-------------------------------------------------------------------------------

results <- results %>%
  group_by(truth_form, fitted_form) %>%
  mutate(iter_diff = n_iter-n_iter[binning == "N8"],
         runtime_ratio = runtime_minutes / runtime_minutes[binning == "N8"]) %>%
  ungroup()

#-------------------------------------------------------------------------------
# Numerical summary
#-------------------------------------------------------------------------------

summary_tab <- results %>%
  group_by(binning) %>%
  summarise(n_fits = n(),
            n_usable = sum(usable),
            n_converged = sum(converged),
            n_hit_max = sum(hit_max),
            n_degenerate = sum(degenerate),
            n_inla_failure = sum(inla_failure),
            n_vb_aborted = sum(vb_aborted),
            n_nan_inf_logl = sum(nan_inf_logl),
            median_iter = median(n_iter),
            max_iter = max(n_iter),
            median_iter_diff = median(iter_diff),
            median_runtime = median(runtime_minutes),
            max_runtime = max(runtime_minutes),
            median_runtime_ratio = median(runtime_ratio),
            .groups = "drop")

write.csv(summary_tab, file.path(binning_table_dir, "binning_summary.csv"),
          row.names=FALSE)

#-------------------------------------------------------------------------------
# Identify fits failing common numerical-quality checks
#-------------------------------------------------------------------------------

problem_fits <- results[!results$usable, ]
write.csv(problem_fits,
          file.path(binning_table_dir, "binning_problem_fits.csv"),
          row.names = FALSE)

if (nrow(problem_fits) > 0) {
  stop("One or more binning-check fits failed the common numerical-quality checks.")
}

#-------------------------------------------------------------------------------
# Extract physical-scale posterior marginals
#-------------------------------------------------------------------------------

temporal_parameters <- list(ou = c("K", "c", "p"), 
                            mse = c("K", "d", "rho", "gamma"),
                            rate_state = c("K", "B", "ta"))

posterior_curves <- list()
j <- 1

for (form in names(temporal_parameters)) {
  
  for (binning in c("N8", "N10", "N12", "N14", "N16")) {
    
    file_i <- file.path(binning_fit_dir,
                        paste0("truth_", form, "_fit_", form, "_", binning, ".rds"))
    obj_i <- readRDS(file_i)
    
    post_i <- ETAS.inlabru::get_posterior_param(
      list(model.fit = obj_i$fit,
           link.functions = obj_i$link.functions,
           form = form))$post.df
    
    post_i <- post_i[post_i$param %in% temporal_parameters[[form]], ]
    
    posterior_curves[[j]] <- data.frame(form = form, 
                                        binning = binning, 
                                        parameter = post_i$param,
                                        value = post_i$x,
                                        density = post_i$y)
    j <- j + 1
  }
}

posterior_density <- do.call(rbind, posterior_curves)

posterior_density$binning <- factor(posterior_density$binning,
                                    levels = c("N8", "N10", "N12", "N14", "N16"))

#-------------------------------------------------------------------------------
# Posterior marginal overlays
#-------------------------------------------------------------------------------

bin_cols <- c(N8  = "firebrick", N10 = "darkorange", N12 = "forestgreen", 
              N14 = "steelblue", N16 = "mediumpurple4")
bin_lty <- c(N8  = "dashed", N10 = "solid", N12 = "dotdash", N14 = "longdash",
             N16 = "dotted")

p_post <- ggplot(
  posterior_density,
  aes(x = value, y = density, colour = binning, linetype = binning)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ form + parameter, scales = "free") +
  scale_colour_manual(values = bin_cols) +
  scale_linetype_manual(values = bin_lty) +
  labs(x = "Parameter value", y = "Posterior density", colour = "N.max", 
       linetype = "N.max") +
  theme_bw()

ggsave(file.path(binning_figure_dir, "binning_posterior_marginals.png"),
       p_post, width = 10, height = 7, dpi = 300)

#-------------------------------------------------------------------------------
# Posterior overlap under sequential refinement (total variation distance)
#-------------------------------------------------------------------------------

calc_overlap <- function(x1, y1, x2, y2, n = 2000) {
  
  grid <- seq(min(c(x1, x2)), max(c(x1, x2)), length.out = n)
  dx <- diff(grid)
  
  d1 <- approx(x1, y1, xout = grid, yleft = 0, yright = 0)$y
  d2 <- approx(x2, y2, xout = grid, yleft = 0, yright = 0)$y
  
  int1 <- sum((d1[-1] + d1[-length(d1)]) * dx / 2)
  int2 <- sum((d2[-1] + d2[-length(d2)]) * dx / 2)
  
  d1 <- d1 / int1
  d2 <- d2 / int2
  
  TV <- 0.5 * sum((abs(d1 - d2)[-1] + abs(d1 - d2)[-length(d1)]) * dx / 2)
  
  TV <- min(max(TV, 0), 1)
  
  c(overlap = 1 - TV, TV = TV)
}

refinements <- list("N8 -> N10" = c("N8", "N10"), 
                    "N10 -> N12" = c("N10", "N12"),
                    "N12 -> N14" = c("N12", "N14"), 
                    "N14 -> N16" = c("N14", "N16"))

overlap_rows <- list()
j <- 1

for (form in names(temporal_parameters)) {
  
  for (parameter in temporal_parameters[[form]]) {
    
    for (comparison in names(refinements)) {
      
      b1 <- refinements[[comparison]][1]
      b2 <- refinements[[comparison]][2]
      
      d1 <- posterior_density[posterior_density$form == form &
                              posterior_density$parameter == parameter &
                              posterior_density$binning == b1, ]
      
      d2 <- posterior_density[posterior_density$form == form &
                              posterior_density$parameter == parameter &
                              posterior_density$binning == b2, ]
      
      d <- calc_overlap(d1$value, d1$density, d2$value, d2$density)
      
      overlap_rows[[j]] <- data.frame(form = form,
                                      parameter = parameter,
                                      comparison = comparison,
                                      overlap = d["overlap"],
                                      TV = d["TV"])
      j <- j + 1
    }
  }
}

posterior_overlap <- do.call(rbind, overlap_rows)
rownames(posterior_overlap) <- NULL

write.csv(posterior_overlap, 
          file.path(binning_table_dir, "posterior_overlap.csv"),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# Functional stability under sequential N.max refinement
#-------------------------------------------------------------------------------

n_samp <- 1000
t_grid <- 10^seq(-3, log10(T_fit_end), length.out = 300)

functional_rows <- list()
j <- 1

set.seed(123)

for (form in c("ou", "mse", "rate_state")) {
  
  for (binning in c("N8", "N10", "N12", "N14", "N16")) {
    
    obj_i <- readRDS(file.path(
      binning_fit_dir,
      paste0("truth_", form, "_fit_", form, "_", binning, ".rds")))
    
    temp_i <- ETAS.inlabru::posterior_temporal_summary(
      list(model.fit = obj_i$fit, link.functions = obj_i$link.functions, form = form),
      t.eval = t_grid, n.samp = n_samp)
    
    functional_rows[[j]] <- temp_i$summary %>%
                            mutate(form = form, binning = binning) %>%
                            select(form, binning, quantity, time, median)
    
    j <- j + 1
  }
}

functional_post <- do.call(rbind, functional_rows)

functional_post$binning <- factor(functional_post$binning,
                                  levels = c("N8", "N10", "N12", "N14", "N16"))

# Functional posterior overlays
p_function <- ggplot(
  functional_post,
  aes(x = time, y = median, colour = binning, linetype = binning)) +
  geom_line(linewidth = 0.8) +
  facet_grid(form ~ quantity, scales = "free_y") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = bin_cols) +
  scale_linetype_manual(values = bin_lty) +
  labs(x = "Time since parent event (days)", y = "Posterior decay function",
       colour = "N.max", linetype = "N.max") +
  theme_bw()

ggsave(file.path(binning_figure_dir, "binning_functional_stability.png"),
       p_function, width = 10, height = 8, dpi = 300)
message("Saved binning-check analysis outputs to ", binning_analysis_dir)
message("Finished temporal binning analysis.")
