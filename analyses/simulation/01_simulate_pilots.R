#===============================================================================
# Generate one pilot catalogue from each kernel
#===============================================================================

#-------------------------------------------------------------------------------
# Load package and fixed simulation design
#-------------------------------------------------------------------------------

library(ETAS.inlabru)
library(here)

source(here::here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Sequential execution for pilot simulations
#-------------------------------------------------------------------------------

future::plan(future::sequential)

#-------------------------------------------------------------------------------
# Function to generate one pilot catalogue for a specified temporal kernel
# Uses the fixed truth, seed and observation window defined in 00_design.R
#-------------------------------------------------------------------------------

simulate_pilot <- function(kernel) {
  
  stopifnot(kernel %in% names(truths))
  
  # Fixed kernel-specific pilot seed
  seed <- pilot_seeds[[kernel]]
  
  set.seed(seed)
  
  message("Generating ", kernel," pilot catalogue (seed = ", seed, ")...")
  
  #-----------------------------------------------------------------------------
  # Generate complete synthetic catalogue
  #-----------------------------------------------------------------------------
  
  catalogue <- generate_temporal_ETAS_synthetic(
                theta = truths[[kernel]], beta.p = beta_true, M0 = M0,
                T1 = T_fit_start, T2 = T_forecast_end, Ht = mainshock_event,
                format = "df", kernel = kernel, Mmax = Mmax)
  
  #-----------------------------------------------------------------------------
  # Check valid catalogues returned
  #-----------------------------------------------------------------------------
  
  stopifnot(
    is.data.frame(catalogue),
    all(c("ts", "magnitudes", "gen") %in% names(catalogue)),
    all(catalogue$ts >= T_fit_start), all(catalogue$ts <= T_forecast_end),
    all(catalogue$magnitudes >= M0), all(catalogue$magnitudes <= Mmax),
    all(diff(catalogue$ts) >= 0))
  
  
  # Check that the imposed mainshock is present only once
  imposed <- catalogue[catalogue$gen == -1, ]
  stopifnot(nrow(imposed) == 1, imposed$ts == 0, imposed$magnitudes == 7.1)
  
  #-----------------------------------------------------------------------------
  # Save catalogue and generating information
  #-----------------------------------------------------------------------------
  
  output <- list(
    catalogue = catalogue, truth_kernel = kernel, 
    truth_parameters = truths[[kernel]], seed = seed,
    design = list(
      M0 = M0, Mmax = Mmax, b = b_true, beta = beta_true,
      fit_start_date = fit_start_date, fit_end_date = fit_end_date,
      forecast_end_date = forecast_end_date, T_fit_start = T_fit_start,
      T_fit_end = T_fit_end, T_forecast_end = T_forecast_end,
      mainshock = mainshock_event, G_mainshock_obs_target = G_mainshock_obs_target,
      G_mainshock_obs = G_mainshock_obs[[kernel]], G_inf = G_inf[[kernel]],
      branching_ratio = branching_ratio[[kernel]], 
      mainshock_direct_obs = mainshock_direct_obs[[kernel]]))
  
  #-----------------------------------------------------------------------------
  # Save
  #-----------------------------------------------------------------------------
  
  outfile <- file.path(pilot_dir, paste0(kernel, "_pilot.rds"))
  
  saveRDS(output, outfile)
  
  message("Saved ", outfile, " | ", nrow(catalogue), " events")
  
  invisible(output)
}

#------------------------------------------------------------------------------
# Generate one pilot from each true kernel
#------------------------------------------------------------------------------

for (kernel in names(truths)) simulate_pilot(kernel)

message("Finished generating all three pilot catalogues.")