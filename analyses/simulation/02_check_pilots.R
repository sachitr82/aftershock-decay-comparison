#===============================================================================
# Check pilot catalogues and validate candidate synthetic truths
#===============================================================================

library(ETAS.inlabru)
library(ggplot2)
library(here)
library(future)

source(here("analyses", "simulation", "00_design.R"))

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Catalogue summary helper
#-------------------------------------------------------------------------------

summarise_catalogue <- function(cat, form) {
  
  data.frame(form = form,
             n_total = nrow(cat),
             n_pre_mainshock = sum(cat$ts < 0),
             n_post_mainshock = sum(cat$ts >= 0),
             max_magnitude = max(cat$magnitudes),
             max_generation = max(cat$gen))
}

#-------------------------------------------------------------------------------
# Summarise the three fixed pilots generated in 01
#-------------------------------------------------------------------------------

pilot_files <- c(ou = file.path(pilot_catalogue_dir, "ou_pilot.rds"),
                 mse = file.path(pilot_catalogue_dir, "mse_pilot.rds"),
                 rate_state = file.path(pilot_catalogue_dir, "rate_state_pilot.rds"))

pilots <- lapply(pilot_files, readRDS)

pilot_summary <- do.call(rbind, lapply(names(pilots), function(form) {
  summarise_catalogue(pilots[[form]]$catalogue, form)
}))

rownames(pilot_summary) <- NULL

#-------------------------------------------------------------------------------
# Generative validation batch - 100 catalogues under each generating form
#-------------------------------------------------------------------------------

n_validation <- 100

validation_index <- expand.grid(form = names(truths), 
                                 rep = 1:n_validation,
                                 KEEP.OUT.ATTRS = FALSE, 
                                 stringsAsFactors = FALSE)

validation_index$seed <- 910000 + seq_len(nrow(validation_index))

validation_results <- vector("list", nrow(validation_index))

for (i in seq_len(nrow(validation_index))) {
  
  form_i <- validation_index$form[i]
  rep_i <- validation_index$rep[i]
  seed_i <- validation_index$seed[i]
  
  set.seed(seed_i)
  
  catalogue_i <- generate_temporal_ETAS_synthetic(theta = truths[[form_i]],
                                                  beta.p = beta_true,
                                                  M0 = M0,
                                                  T1 = T_fit_start, 
                                                  T2 = T_fit_end, 
                                                  Ht = mainshock_event,
                                                  format = "df",
                                                  form = form_i,
                                                  Mmax = Mmax)
  
  out_i <- summarise_catalogue(catalogue_i, form_i)
  out_i$rep <- rep_i
  out_i$seed <- seed_i
  
  validation_results[[i]] <- out_i
  
  message(form_i, " rep ", rep_i, ": ", nrow(catalogue_i), " events")
}

validation_results <- do.call(rbind, validation_results)
rownames(validation_results) <- NULL

#-------------------------------------------------------------------------------
# Summarise validation distributions
#-------------------------------------------------------------------------------

validation_summary <- do.call(rbind, lapply(names(truths), function(form){

  x <- validation_results[validation_results$form == form, ]
  
  data.frame(form = form,
            
             pre_q10 = quantile(x$n_pre_mainshock, 0.10),
             pre_median = median(x$n_pre_mainshock),
             pre_q90 = quantile(x$n_pre_mainshock, 0.90),
             
             post_q10 = quantile(x$n_post_mainshock, 0.10),
             post_median = median(x$n_post_mainshock),
             post_q90 = quantile(x$n_post_mainshock, 0.90),
             
             total_q10 = quantile(x$n_total, 0.10),
             total_median = median(x$n_total),
             total_q90 = quantile(x$n_total, 0.90),
            
             max_generation_median = median(x$max_generation))
}))

rownames(validation_summary) <- NULL

#-------------------------------------------------------------------------------
# Pilot plots
#-------------------------------------------------------------------------------

pilot_catalogues <- do.call(rbind, lapply(names(pilots), function(form) {
  x <- pilots[[form]]$catalogue
  x$form <- form
  x
}))

pilot_catalogues$form <- factor(pilot_catalogues$form,
                                  levels = c("ou", "mse", "rate_state"),
                                  labels = c("OU", "MSE", "Rate-state"))

p_mag <- ggplot(pilot_catalogues, aes(ts, magnitudes)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~form, ncol = 1) +
  labs(x = "Time relative to mainshock (days)", y = "Magnitude") +
  theme_bw()

ggsave(file.path(pilot_validation_dir, "pilot_magnitude_time.pdf"), p_mag,
       width = 8, height = 7)

#-------------------------------------------------------------------------------
# Save summaries
#-------------------------------------------------------------------------------

write.csv(pilot_summary, file.path(pilot_validation_dir, "pilot_summary.csv"), 
          row.names = FALSE)
write.csv(validation_results, file.path(pilot_validation_dir, "validation_batch_counts.csv"), 
          row.names = FALSE)
write.csv(validation_summary, file.path(pilot_validation_dir, "validation_batch_summary.csv"), 
          row.names = FALSE)
message("Saved pilot-validation outputs to ", pilot_validation_dir)
message("Finished pilot validation.")