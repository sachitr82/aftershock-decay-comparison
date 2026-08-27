#===============================================================================
# Analyse rel_tol check
#===============================================================================

library(ggplot2)
library(dplyr)
library(here)

source(here::here("analyses", "simulation", "00_design.R"))

rel_tol_dir <- file.path(fit_dir, "rel_tol_check")
analysis_dir <- file.path(rel_tol_dir, "analysis")
dir.create(analysis_dir,recursive=TRUE, showWarnings=FALSE)

#-------------------------------------------------------------------------------
# Load numerical results from 05
#-------------------------------------------------------------------------------

results <- read.csv(file.path(rel_tol_dir, "rel_tol_manifest.csv"))
results$tolerance <- factor(results$tolerance, levels = 
                              c("RT01", "RT005", "RT001"))

stopifnot(nrow(results) == 9)

#-------------------------------------------------------------------------------
# Paired changes relative to rel_tol = 0.1
#-------------------------------------------------------------------------------

results <- results %>%
  group_by(truth_kernel, fitted_kernel) %>%
  mutate(iter_diff = n_iter-n_iter[tolerance == "RT01"],
         runtime_ratio = runtime_minutes / runtime_minutes[tolerance == "RT01"]) %>%
  ungroup()

#-------------------------------------------------------------------------------
# Numerical summary
#-------------------------------------------------------------------------------

summary_tab <- results %>%
  group_by(tolerance) %>%
  summarise(
    n_fits=n(),
    n_converged = sum(converged),
    n_hit_max = sum(hit_max),
    median_iter = median(n_iter),
    max_iter = max(n_iter),
    median_iter_diff = median(iter_diff),
    median_runtime = median(runtime_minutes),
    max_runtime = max(runtime_minutes),
    median_runtime_ratio=median(runtime_ratio),
    .groups="drop"
  )

print(summary_tab)
write.csv(summary_tab, file.path(analysis_dir, "rel_tol_summary.csv"),
          row.names=FALSE)

#-------------------------------------------------------------------------------
# Identify problematic fits (didn't converge or hit max iter)
#-------------------------------------------------------------------------------

problem_fits <- results[!results$converged | results$hit_max,]
print(problem_fits)
write.csv(problem_fits, file.path(analysis_dir, "rel_tol_problem_fits.csv"), 
          row.names=FALSE)

#-------------------------------------------------------------------------------
# Extract physical-scale posterior marginals
#-------------------------------------------------------------------------------

temporal_parameters <- list(
  ou = c("K", "c", "p"),
  mse = c("K", "d", "rho", "gamma"),
  rate_state = c("K", "B", "ta"))

posterior_curves <- list()
j <- 1

for (kernel in names(temporal_parameters)) {
  
  for (tolerance in c("RT01", "RT005", "RT001")) {
    
    file_i <- file.path(
      rel_tol_dir,
      paste0("truth_", kernel, "_fit_", kernel, "_", tolerance, ".rds"))
    obj_i <- readRDS(file_i)
    
    post_i <- ETAS.inlabru::get_posterior_param(list(
      model.fit = obj_i$fit,
      link.functions = obj_i$link.functions,
      kernel = kernel))$post.df
    
    post_i <- post_i[
      post_i$param %in% temporal_parameters[[kernel]], ]
    
    posterior_curves[[j]] <- data.frame(
      kernel = kernel,
      tolerance = tolerance,
      parameter = post_i$param,
      value = post_i$x,
      density = post_i$y)
    
    j <- j + 1
  }
}

posterior_density <- do.call(rbind, posterior_curves)

posterior_density$tolerance <- factor(
  posterior_density$tolerance,
  levels = c("RT01", "RT005", "RT001"))

#-------------------------------------------------------------------------------
# Posterior marginal overlays
#-------------------------------------------------------------------------------

tol_cols <- c(RT01  = "firebrick", RT005 = "forestgreen", RT001 = "steelblue")
tol_lty <- c(RT01 = "dashed", RT005 = "solid", RT001 = "dotted")
tol_labels <- c(RT01 = "0.1", RT005 = "0.05", RT001 = "0.01")

p_post <- ggplot(
  posterior_density,
  aes(x = value, y = density, colour = tolerance, linetype = tolerance)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ kernel + parameter, scales = "free") +
  scale_colour_manual(values = tol_cols, labels = tol_labels) +
  scale_linetype_manual(values = tol_lty, labels = tol_labels) +
  labs(x = "Parameter value", y = "Posterior density", colour = "rel_tol", 
       linetype = "rel_tol") +
  theme_bw()

print(p_post)

ggsave(file.path(analysis_dir, "rel_tol_posterior_marginals.png"),
       p_post, width = 10, height = 7, dpi = 300)

#-------------------------------------------------------------------------------
# Posterior overlap under sequential refinement
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

refinements <- list("0.1 -> 0.05" = c("RT01", "RT005"), 
                    "0.05 -> 0.01" = c("RT005", "RT001"))

overlap_rows <- list()
j <- 1

for (kernel in names(temporal_parameters)) {
  
  for (parameter in temporal_parameters[[kernel]]) {
    
    for (comparison in names(refinements)) {
      
      t1 <- refinements[[comparison]][1]
      t2 <- refinements[[comparison]][2]
      
      d1 <- posterior_density[
        posterior_density$kernel == kernel &
          posterior_density$parameter == parameter &
          posterior_density$tolerance == t1, ]
      
      d2 <- posterior_density[
        posterior_density$kernel == kernel &
          posterior_density$parameter == parameter &
          posterior_density$tolerance == t2, ]
      
      d <- calc_overlap(d1$value, d1$density, d2$value, d2$density)
      
      overlap_rows[[j]] <- data.frame(
        kernel = kernel,
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

print(posterior_overlap)

write.csv(
  posterior_overlap,
  file.path(analysis_dir, "posterior_overlap.csv"),
  row.names = FALSE)

#-------------------------------------------------------------------------------
# # Functional stability under sequential rel_tol refinement
#-------------------------------------------------------------------------------

n_samp <- 1000
t_grid <- 10^seq(-3, log10(T_fit_end), length.out = 300)

functional_rows <- list()
j <- 1

set.seed(123)

for (kernel in c("ou", "mse", "rate_state")) {
  
  for (tolerance in c("RT01", "RT005", "RT001")) {
    
    obj_i <- readRDS(file.path(
      rel_tol_dir,
      paste0("truth_", kernel, "_fit_", kernel, "_", tolerance, ".rds")))
    
    samp_i <- ETAS.inlabru::post_sampling(
      list(
        model.fit = obj_i$fit,
        link.functions = obj_i$link.functions,
        kernel = kernel),
      n.samp = n_samp)
    
    g_mat <- matrix(NA_real_, nrow = n_samp, ncol = length(t_grid))
    
    for (s in seq_len(n_samp)) {
      
      g_mat[s, ] <- ETAS.inlabru::temporal_kernel(
        dt = t_grid, theta = as.list(samp_i[s, ]),kernel = kernel)
    }
    
    Kg_mat <- sweep(g_mat, 1, samp_i$K, "*")
    
    functional_rows[[j]] <- rbind(
      data.frame(
        kernel = kernel,
        tolerance = tolerance,
        quantity = "g(t)",
        time = t_grid,
        median = apply(g_mat, 2, median)),
      
      data.frame(
        kernel = kernel,
        tolerance = tolerance,
        quantity = "K*g(t)",
        time = t_grid,
        median = apply(Kg_mat, 2, median))
    )
    
    j <- j + 1
  }
}

functional_post <- do.call(rbind, functional_rows)

functional_post$tolerance <- factor(
  functional_post$tolerance,
  levels = c("RT01", "RT005", "RT001"))

#-------------------------------------------------------------------------------
# Functional posterior overlays
#-------------------------------------------------------------------------------

p_function <- ggplot(
  functional_post,
  aes(x = time, y = median, colour = tolerance, linetype = tolerance)) +
  geom_line(linewidth = 0.8) +
  facet_grid(kernel ~ quantity, scales = "free_y") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = tol_cols, labels = tol_labels) +
  scale_linetype_manual(values = tol_lty, labels = tol_labels) +
  labs(x = "Time since parent event (days)", y = "Posterior triggering function",
       colour = "rel_tol", linetype = "rel_tol") +
  theme_bw()

print(p_function)

ggsave(file.path(analysis_dir, "rel_tol_functional_stability.png"),
       p_function, width = 10, height = 8, dpi = 300)
