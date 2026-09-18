# Report materials

This directory contains the final MSc Statistics dissertation **Comparing Aftershock Decay Functions for the 2019 Ridgecrest Earthquake Sequence,** and the source materials used to compile it.

The dissertation was prepared in Overleaf. A read-only version of the Overleaf project is available at:

<https://www.overleaf.com/read/qjtvgstygjry#bc5cb0>

## Contents

- `source/` — a copy of the final Overleaf project source, including the LaTeX files, bibliography, figures, class files, and other files required to compile the report
- the compiled dissertation PDF, retained in this directory

The statistical analysis underlying the dissertation is contained elsewhere in the repository. Analysis scripts are under `../analyses/`, shared functions are under `../src/`, data are under `../data/`, and retained results and figures are under `../outputs/`.

Figures in `source/images/` are copies used for report compilation. Their original generated versions are retained under `../outputs/`. The table below records the provenance of each report figure.

## Figure provenance

| Report figure | Original generated file | Generating script |
|-----------------|----------------------------------|---------------------|
| `decay-comparison.pdf` | `outputs/comparison-figure/decay-comparison.pdf` | `analyses/comparison-figure/01_plot-decay-comparison.R` |
| `ridgecrest_completeness.pdf` | `outputs/eda/figures/ridgecrest_completeness.pdf` | `analyses/eda/eda_ridgecrest.qmd` |
| `ridgecrest_window.pdf` | `outputs/eda/figures/ridgecrest_window.pdf` | `analyses/eda/eda_ridgecrest.qmd` |
| `synthetic-decay-design.pdf` | `outputs/simulation/design_figures/synthetic-decay-design.pdf` | `analyses/simulation/00b_true_plots.R` |
| `convergence_3x3.pdf` | `outputs/simulation/analysis/convergence/figures/convergence_3x3.pdf` | `analyses/simulation/07_convergence_analysis.R` |
| `parameter_posterior_overlays.pdf` | `outputs/simulation/analysis/parameter_recovery/figures/parameter_posterior_overlays.pdf` | `analyses/simulation/08_posterior_parameter_recovery.R` |
| `functional_recovery.pdf` | `outputs/simulation/analysis/parameter_recovery/figures/functional_recovery.pdf` | `analyses/simulation/08_posterior_parameter_recovery.R` |
| `posterior_parameter_correlations.pdf` | `outputs/simulation/analysis/parameter_recovery/figures/posterior_parameter_correlations.pdf` | `analyses/simulation/08a_posterior_temporal_dependence.R` |
| `posterior_functional_tv.pdf` | `outputs/simulation/analysis/misspecification/figures/posterior_functional_tv.pdf` | `analyses/simulation/09_temporal_form_misspecification.R` |
| `pairwise_model_comparison.pdf` | `outputs/simulation/analysis/model_discrimination/figures/pairwise_model_comparison.pdf` | `analyses/simulation/10_model_discrimination.R` |
| `ridgecrest-functional-posteriors.pdf` | `outputs/ridgecrest/baseline/temporal_functions/figures/ridgecrest-functional-posteriors.pdf` | `analyses/ridgecrest/03_temporal_functions.R` |

To regenerate a report figure, run the corresponding analysis script and use the resulting file under `../outputs/`.

See the root `README.md` for instructions for reproducing the analysis.
