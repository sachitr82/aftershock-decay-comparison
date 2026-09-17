#===============================================================================
# Generate final synthetic catalogue ensemble
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(here)
library(future)

source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Sequential execution
#-------------------------------------------------------------------------------

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Check final simulation index
#-------------------------------------------------------------------------------

stopifnot(nrow(sim_index) == 3 * n_rep,
          !anyDuplicated(sim_index$seed),
          setequal(sim_index$form, names(truths)),
          all(table(sim_index$form) == n_rep))

#-------------------------------------------------------------------------------
# Storage for simulation manifest
#-------------------------------------------------------------------------------

simulation_manifest <- vector("list", nrow(sim_index))

#-------------------------------------------------------------------------------
# Generate final catalogues
#-------------------------------------------------------------------------------

for (i in seq_len(nrow(sim_index))) {
  
  form_i <- sim_index$form[i]
  rep_i <- sim_index$rep[i]
  seed_i <- sim_index$seed[i]
  
  set.seed(seed_i)
  
  message("Generating ", form_i, " catalogue ", rep_i, "/", n_rep,
          " (seed = ", seed_i, ")...")
  
  catalogue_i <- generate_temporal_ETAS_synthetic(theta = truths[[form_i]],
                                                  beta.p = beta_true,
                                                  M0 = M0,
                                                  T1 = T_fit_start,
                                                  T2 = T_fit_end,
                                                  Ht = mainshock_event,
                                                  format = "df",
                                                  form = form_i,
                                                  Mmax = Mmax)
  
  #-----------------------------------------------------------------------------
  # Structural checks
  #-----------------------------------------------------------------------------
  
  stopifnot(
    is.data.frame(catalogue_i),
    all(c("ts", "magnitudes", "gen") %in% names(catalogue_i)),
    all(catalogue_i$ts >= T_fit_start), all(catalogue_i$ts <= T_fit_end),
    all(catalogue_i$magnitudes >= M0), all(catalogue_i$magnitudes <= Mmax),
    all(diff(catalogue_i$ts) >= 0)
  )
  
  imposed_i <- catalogue_i[catalogue_i$gen == -1, ]
  stopifnot(nrow(imposed_i) == 1, imposed_i$ts == 0, imposed_i$magnitudes == 7.1)
  
  #-----------------------------------------------------------------------------
  # Save catalogue and generating information
  #-----------------------------------------------------------------------------
  
  outfile_i <- file.path(catalogue_dir, sprintf("%s_%03d.rds", form_i, rep_i))
  
  output_i <- list(catalogue = catalogue_i,
                   truth_form = form_i,
                   truth_parameters = truths[[form_i]],
                   rep = rep_i,
                   seed = seed_i,
                   design = list(M0 = M0,
                                 Mmax = Mmax,
                                 b = b_true,
                                 beta = beta_true,
                                 fit_start_date = fit_start_date,
                                 fit_end_date = fit_end_date,
                                 T_fit_start = T_fit_start,
                                 T_fit_end = T_fit_end,
                                 mainshock = mainshock_event,
                                 n_ou_anchor = n_ou_anchor,
                                 G_mainshock_obs_target = G_mainshock_obs_target,
                                 G_mainshock_obs = G_mainshock_obs[[form_i]],
                                 G_inf = G_inf[[form_i]],
                                 branching_ratio = branching_ratio[[form_i]],
                                 mainshock_direct_obs = mainshock_direct_obs[[form_i]]))
  
  saveRDS(output_i, outfile_i)
  
  simulation_manifest[[i]] <- data.frame(form = form_i, 
                                         rep = rep_i,
                                         seed = seed_i,
                                         n_events = nrow(catalogue_i),
                                         file = basename(outfile_i))
  
  message("Saved ", basename(outfile_i), " | ", nrow(catalogue_i), " events")
}

#-------------------------------------------------------------------------------
# Save simulation manifest
#-------------------------------------------------------------------------------

simulation_manifest <- do.call(rbind, simulation_manifest)
rownames(simulation_manifest) <- NULL

manifest_file <- file.path(catalogue_dir, "simulation_manifest.csv")

write.csv(simulation_manifest, manifest_file, row.names = FALSE)

message("Saved ", manifest_file)

#-------------------------------------------------------------------------------
# Final completeness checks
#-------------------------------------------------------------------------------

stopifnot(nrow(simulation_manifest) == 3 * n_rep,
          all(table(simulation_manifest$form) == n_rep),
          all(file.exists(file.path(catalogue_dir, simulation_manifest$file))))

message("Finished generating all ",
        nrow(simulation_manifest), " final synthetic catalogues.")