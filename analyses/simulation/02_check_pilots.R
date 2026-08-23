#===============================================================================
# Check pilot catalogues and validate candidate synthetic truths
#===============================================================================
library(ETAS.inlabru)
library(ggplot2)
library(here)

source(here::here("analyses", "simulation", "00_design.R"))

future::plan(future::sequential)

# -------------------------------------------------------------------------
# Catalogue summary helper
# -------------------------------------------------------------------------

summarise_catalogue <- function(cat, kernel) {
  
  data.frame(
    kernel = kernel,
    n_fit = sum(cat$ts > T_fit_start &cat$ts < T_fit_end),
    n_pre_mainshock =  sum(cat$ts > T_fit_start & cat$ts < 0),
    n_post_mainshock = sum(cat$ts >= 0 & cat$ts < T_fit_end),
    n_forecast = sum(cat$ts >= T_fit_end & cat$ts < T_forecast_end),
    n_total =  sum(cat$ts > T_fit_start & cat$ts < T_forecast_end),
    max_magnitude = max(cat$magnitudes),
    max_generation = max(cat$gen)
  )
}

#===============================================================================
# Summarise the three fixed pilots generated in 01
#===============================================================================

pilot_files <- c(
  ou = here::here("results", "simulation", "pilots", "ou_pilot.rds"),
  mse = here::here("results", "simulation", "pilots", "mse_pilot.rds"),
  rate_state = here::here(
    "results", "simulation", "pilots", "rate_state_pilot.rds"))

pilots <- lapply(pilot_files, readRDS)

pilot_summary <- do.call(rbind, lapply(names(pilots), function(kernel) {
  summarise_catalogue(pilots[[kernel]]$catalogue, kernel)
}))

rownames(pilot_summary) <- NULL
print(pilot_summary)

#-------------------------------------------------------------------------------
# Largest randomly generated events; understanding why pilot catalogue counts
# differ
#-------------------------------------------------------------------------------

largest_events <- do.call(rbind, lapply(names(pilots), function(kernel) {
  
  x <- pilots[[kernel]]$catalogue
  x <- x[x$gen != -1, ]
  x <- x[order(x$magnitudes, decreasing = TRUE), ]
  x <- head(x, 10)
  x$kernel <- kernel
  x[, c("kernel", "ts", "magnitudes", "gen")]
}))

rownames(largest_events) <- NULL
print(largest_events)

#===============================================================================
# Temporary generative validation batch
#===============================================================================

n_calibration <- 100

calibration_index <- expand.grid(
  kernel = names(truths), rep = 1:n_calibration,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

calibration_index$seed <- 910000 + seq_len(nrow(calibration_index))

calibration_results <- vector("list", nrow(calibration_index))

for (i in seq_len(nrow(calibration_index))) {
  
  kernel_i <- calibration_index$kernel[i]
  rep_i <- calibration_index$rep[i]
  seed_i <- calibration_index$seed[i]
  
  set.seed(seed_i)
  
  catalogue_i <- generate_temporal_ETAS_synthetic(
    theta = truths[[kernel_i]], beta.p = beta_true, M0 = M0,
    T1 = T_fit_start, T2 = T_forecast_end, Ht = mainshock_event,
    format = "df", kernel = kernel_i, Mmax = Mmax)
  
  out_i <- summarise_catalogue(catalogue_i, kernel_i)
  out_i$rep <- rep_i
  out_i$seed <- seed_i
  
  calibration_results[[i]] <- out_i
  
  message(kernel_i, " rep ", rep_i, ": ", nrow(catalogue_i), " events")
}

calibration_results <- do.call(rbind, calibration_results)
rownames(calibration_results) <- NULL

#===============================================================================
# Summarise calibration distributions
#===============================================================================

calibration_summary <- do.call(rbind, lapply(names(truths), function(kernel){

  x <- calibration_results[calibration_results$kernel == kernel, ]
  
  data.frame(
    kernel = kernel,
    
    pre_q10 = quantile(x$n_pre_mainshock, 0.10),
    pre_median = median(x$n_pre_mainshock),
    pre_q90 = quantile(x$n_pre_mainshock, 0.90),
    
    post_q10 = quantile(x$n_post_mainshock, 0.10),
    post_median = median(x$n_post_mainshock),
    post_q90 = quantile(x$n_post_mainshock, 0.90),
    
    fit_q10 = quantile(x$n_fit, 0.10),
    fit_median = median(x$n_fit),
    fit_q90 = quantile(x$n_fit, 0.90),
    
    total_q10 = quantile(x$n_total, 0.10),
    total_median = median(x$n_total),
    total_q90 = quantile(x$n_total, 0.90),
    
    max_generation_median = median(x$max_generation))
}))

rownames(calibration_summary) <- NULL
print(calibration_summary)

#===============================================================================
# Pilot plots
#===============================================================================

pilot_catalogues <- do.call(rbind, lapply(names(pilots), function(kernel) {
  x <- pilots[[kernel]]$catalogue
  x$kernel <- kernel
  x
}))

pilot_catalogues$kernel <- factor(
  pilot_catalogues$kernel,
  levels = c("ou", "mse", "rate_state"),
  labels = c("OU", "MSE", "Rate-state"))

p_mag <- ggplot(pilot_catalogues, aes(ts, magnitudes)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = T_fit_end, linetype = "dotted") +
  facet_wrap(~kernel, ncol = 1) +
  labs(x = "Time relative to mainshock (days)", y = "Magnitude") +
  theme_bw()

ggsave(file.path(pilot_dir, "pilot_magnitude_time.pdf"), p_mag, width = 8, height = 7)

#===============================================================================
# E. Save summaries
#===============================================================================

write.csv(pilot_summary, file.path(pilot_dir, "pilot_summary.csv"), 
          row.names = FALSE)
write.csv(largest_events, file.path(pilot_dir, "pilot_largest_events.csv"), 
          row.names = FALSE)
write.csv(calibration_results, file.path(pilot_dir, "calibration_batch_counts.csv"), 
          row.names = FALSE)
write.csv(calibration_summary, file.path(pilot_dir, "calibration_batch_summary.csv"), 
          row.names = FALSE)