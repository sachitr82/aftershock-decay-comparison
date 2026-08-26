#===============================================================================
# Generate final synthetic catalogue ensemble
#===============================================================================

#-------------------------------------------------------------------------------
# Load packages and design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(here)

source(here::here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Sequential execution
#-------------------------------------------------------------------------------

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Check final simulation index
#-------------------------------------------------------------------------------

stopifnot(
  nrow(sim_index) == 3 * n_rep,
  !anyDuplicated(sim_index$seed),
  setequal(sim_index$kernel, names(truths)),
  all(table(sim_index$kernel) == n_rep)
)

#-------------------------------------------------------------------------------
# Storage for simulation manifest
#-------------------------------------------------------------------------------

simulation_manifest <- vector("list", nrow(sim_index))

#===============================================================================
# Generate final catalogues
#===============================================================================

for (i in seq_len(nrow(sim_index))) {
  
  kernel_i <- sim_index$kernel[i]
  rep_i <- sim_index$rep[i]
  seed_i <- sim_index$seed[i]
  
  set.seed(seed_i)
  
  message("Generating ", kernel_i, " catalogue ", rep_i, "/", n_rep,
          " (seed = ", seed_i, ")...")
  
  catalogue_i <- generate_temporal_ETAS_synthetic(
    theta = truths[[kernel_i]], beta.p = beta_true, M0 = M0,
    T1 = T_fit_start, T2 = T_fit_end, Ht = mainshock_event,
    format = "df", kernel = kernel_i, Mmax = Mmax)
  
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
  
  outfile_i <- file.path(catalogue_dir, sprintf("%s_%03d.rds", kernel_i, rep_i))
  
  output_i <- list(
    catalogue = catalogue_i,
    truth_kernel = kernel_i,
    truth_parameters = truths[[kernel_i]],
    rep = rep_i,
    seed = seed_i,
    design = list(
      M0 = M0, Mmax = Mmax, b = b_true, beta = beta_true,
      fit_start_date = fit_start_date, fit_end_date = fit_end_date,
      T_fit_start = T_fit_start, T_fit_end = T_fit_end,
      mainshock = mainshock_event,
      n_ou_anchor = n_ou_anchor,
      G_mainshock_obs_target = G_mainshock_obs_target,
      G_mainshock_obs = G_mainshock_obs[[kernel_i]],
      G_inf = G_inf[[kernel_i]],
      branching_ratio = branching_ratio[[kernel_i]],
      mainshock_direct_obs = mainshock_direct_obs[[kernel_i]]
    )
  )
  
  saveRDS(output_i, outfile_i)
  
  simulation_manifest[[i]] <- data.frame(
    kernel = kernel_i, rep = rep_i, seed = seed_i,
    n_events = nrow(catalogue_i),
    file = basename(outfile_i)
  )
  
  message("Saved ", basename(outfile_i), " | ", nrow(catalogue_i), " events")
}

#===============================================================================
# Save simulation manifest
#===============================================================================

simulation_manifest <- do.call(rbind, simulation_manifest)
rownames(simulation_manifest) <- NULL

write.csv(simulation_manifest, 
          file.path(catalogue_dir, "simulation_manifest.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Final completeness checks
#-------------------------------------------------------------------------------

stopifnot(
  nrow(simulation_manifest) == 3 * n_rep,
  all(table(simulation_manifest$kernel) == n_rep),
  all(file.exists(file.path(catalogue_dir, simulation_manifest$file)))
)

message("Finished generating all ",
        nrow(simulation_manifest), " final synthetic catalogues.")