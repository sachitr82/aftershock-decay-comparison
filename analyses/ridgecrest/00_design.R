#===============================================================================
# Ridgecrest analysis design
#===============================================================================

library(here)

# Shared model, prior and numerical settings
source(here("analyses", "simulation", "00_design.R"))

#-------------------------------------------------------------------------------
# Output directories
#-------------------------------------------------------------------------------

ridgecrest_output_dir <- here("outputs", "ridgecrest")

baseline_dir <- file.path(ridgecrest_output_dir, "baseline")

baseline_fit_dir <- file.path(baseline_dir, "fits")
baseline_fit_summary_dir <- file.path(baseline_dir, "fit_summary")

parameter_dir <- file.path(baseline_dir, "parameters")
parameter_table_dir <- file.path(parameter_dir, "tables")
parameter_figure_dir <- file.path(parameter_dir, "figures")

temporal_function_dir <- file.path(baseline_dir, "temporal_functions")
temporal_function_table_dir <- file.path(temporal_function_dir, "tables")
temporal_function_figure_dir <- file.path(temporal_function_dir, "figures")

model_comparison_dir <- file.path(baseline_dir, "model_comparison")
model_comparison_table_dir <- file.path(model_comparison_dir, "tables")

sensitivity_dir <- file.path(ridgecrest_output_dir, "sensitivity")

magnitude_dir <- file.path(sensitivity_dir, "magnitude")
magnitude_fit_dir <- file.path(magnitude_dir, "fits")
magnitude_table_dir <- file.path(magnitude_dir, "tables")

prior_dir <- file.path(sensitivity_dir, "prior")
prior_fit_dir <- file.path(prior_dir, "fits")
prior_table_dir <- file.path(prior_dir, "tables")

#-------------------------------------------------------------------------------
# Create leaf directories
#-------------------------------------------------------------------------------

dir.create(baseline_fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(baseline_fit_summary_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(parameter_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(parameter_figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(temporal_function_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temporal_function_figure_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(model_comparison_table_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(magnitude_fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(magnitude_table_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(prior_fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prior_table_dir, recursive = TRUE, showWarnings = FALSE)