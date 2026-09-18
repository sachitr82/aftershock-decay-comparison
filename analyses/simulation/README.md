# Simulation study

This directory contains the simulation study comparing the Omori–Utsu (OU), modified stretched exponential (MSE), and rate-state (RS) temporal ETAS forms.

The study generates 100 catalogues from each temporal form (300 total) and fits all three candidate models to every catalogue, giving 900 production fits.

## Workflow

Scripts are numbered in their intended order:

| Script | Purpose | Compute note |
|------------------------|------------------------|------------------------|
| `00_design.R` | Defines simulation truths, priors, seeds, fitting controls, and output paths | Sourced automatically by other scripts |
| `00a_prior_calibration.R` | Checks induced temporal priors | Lightweight |
| `00b_true_plots.R` | Plots the generating temporal decay functions | Lightweight |
| `01_simulate_pilots.R` | Generates one pilot catalogue per form | Lightweight |
| `02_check_pilots.R` | Checks pilots and performs generative validation | — |
| `03_binning_check.R` | Tests five temporal-binning settings | \~15 min, \~150 MB of fits |
| `03a_binning_analysis.R` | Analyses binning sensitivity | Requires fits from `03` |
| `04_rel_tol_check.R` | Tests three convergence tolerances | \~15 min, \~90 MB of fits |
| `04a_rel_tol_analysis.R` | Analyses tolerance sensitivity | Requires fits from `04` |
| `05_simulate_all.R` | Generates the final 300 synthetic catalogues | Catalogues are retained in Git |
| `06_fit_simulation.R` | Fits OU, MSE and RS to all 300 catalogues | **\~20 h, \~10–12 GB, 6 workers on Apple M4** |
| `07_convergence_analysis.R` | Summarises convergence and runtimes | Uses the retained fit manifest |
| `07a_misspecification_diagnostics.R` | Targeted investigation of difficult RS/MSE misspecification | **\~2 h, \~200 MB** |
| `08_posterior_parameter_recovery.R` | Parameter and functional recovery under correct specification | **\~30–60 min** |
| `08a_posterior_temporal_dependence.R` | Posterior parameter-dependence analysis | — |
| `09_temporal_form_misspecification.R` | Functional recovery under misspecification | **\~30 min; requires output from `08`** |
| `10_model_discrimination.R` | DIC and log-marginal-likelihood model comparison | — |

Runtime estimates are shown only for the more computationally demanding stages. They were recorded on the Apple M4 MacBook Air used for the dissertation analysis and will vary with hardware and system configuration.

## Fitted-model objects

Large fitted-model `.rds` files beneath `fits/` directories are intentionally excluded from Git because they are reproducible intermediate files and would add substantial storage overhead. Individual fitted-model objects are typically around 10 MB, with the 900 production fits requiring approximately 10–12 GB in total.

The repository instead retains the synthetic catalogues, fitting manifests, summaries, tables and figures needed to document the original analysis. Scripts that require fitted-model objects, including `07a` and `08`–`10`, therefore require the corresponding fits to be regenerated first.

`07_convergence_analysis.R` is an exception because it operates from the retained `fit_manifest.csv`.

## Reproducing the main simulation analysis

The final 300 synthetic catalogues are included in the repository, so reproducing the main simulation analysis can begin with `06_fit_simulation.R`.

After the production fits have been regenerated, scripts `07`–`10` reproduce the downstream convergence, recovery, misspecification and model-discrimination analyses. `07a_misspecification_diagnostics.R` reproduces the additional convergence investigation and can be run separately.

To reproduce the simulation design and numerical checks as well, run scripts `00a`–`05` in filename order before `06`.

Software versions used for the dissertation analysis are recorded in `outputs/session_info.txt`.

The shared fitting helpers are in `src/fit_helpers/`, and the modified temporal ETAS implementation is provided in the `software/ETAS.inlabru` submodule.
