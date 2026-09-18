# Comparing Aftershock Decay Functions for the 2019 Ridgecrest Earthquake Sequence

Code and reproducibility materials for the MSc Statistics dissertation by Sachit Ramakrishnan, **Comparing Aftershock Decay Functions for the 2019 Ridgecrest Earthquake Sequence**, at Imperial College London.

This project compares three temporal decay specifications within a Bayesian temporal epidemic-type aftershock sequence (ETAS) model:

- **Omori–Utsu (OU)**
- **modified stretched exponential (MSE)**
- **rate-state (RS)**

The project extends `ETAS.inlabru` beyond its original OU implementation, evaluates the three forms in a known-truth simulation study, and applies the extended framework to the 2019 Ridgecrest earthquake sequence.

## Repository structure

| Path | Contents |
|------------------------------------|------------------------------------|
| `analyses/eda/` | Ridgecrest catalogue exploration and preparation |
| `analyses/simulation/` | Simulation design, fitting and downstream analyses |
| `analyses/ridgecrest/` | Empirical Ridgecrest fitting and sensitivity analyses |
| `analyses/comparison-figure/` | Comparison of the three temporal decay functions |
| `data/raw/` | Original SCEDC/SCSN catalogue and provenance |
| `data/derived/` | Analysis-ready Ridgecrest catalogue |
| `src/` | Shared EDA and fitting helper functions |
| `software/ETAS.inlabru/` | Modified `ETAS.inlabru` package used in this project |
| `outputs/` | Retained catalogues, manifests, tables, diagnostics and figures |
| `reports/` | Report-related materials |

Detailed workflows are documented in:

- `analyses/simulation/README.md`
- `analyses/ridgecrest/README.md`
- `outputs/README.md`
- `software/ETAS.inlabru/README.md`

## Software setup

### 1. Clone the repository

Clone the dissertation repository and run analysis scripts from the repository root:

``` bash
git clone https://github.com/sachitr82/aftershock-decay-comparison.git
cd aftershock-decay-comparison
```

The modified `ETAS.inlabru` package is installed separately in Step 3. Its source is also linked under `software/ETAS.inlabru/` as a Git submodule for reference.

### 2. Install R dependencies

The analysis requires `inlabru`, R-INLA and the dissertation version of `ETAS.inlabru`, together with standard analysis packages.

Install `inlabru`:

``` r
install.packages("inlabru")
```

Install R-INLA:

``` r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"),
  dep = TRUE)
```

Install the remaining analysis packages:

``` r
install.packages(c("remotes",
                    "here",
                    "tidyverse",
                    "future",
                    "future.apply",
                    "patchwork",
                    "zoo",
                    "sf",
                    "maps",
                    "lubridate",
                    "ggrepel",
                    "geosphere",
                    "scales",
                    "knitr"))
```

The Ridgecrest EDA is a Quarto document, so Quarto is also required to render `analyses/eda/eda_ridgecrest.qmd`.

### 3. Install the modified `ETAS.inlabru`

Install the dissertation fork from GitHub:

``` r
remotes::install_github("sachitr82/ETAS.inlabru")
```

The package modifications and upstream `ETAS.inlabru` documentation are provided in `software/ETAS.inlabru/README.md`.

The recorded analysis environment used R 4.6.1, `ETAS.inlabru` 1.1.1.9003, `inlabru` 2.14.1 and INLA 26.06.08. Full session information is stored in `outputs/session_info.txt`.

## Reproducing the analysis

The repository already contains:

- the raw SCEDC/SCSN earthquake catalogue,
- the derived Ridgecrest catalogue,
- and the final 300 synthetic catalogues.

This means the full workflow can be reproduced from the retained raw inputs, while the main fitting analyses can begin from the retained derived inputs.

### 1. Ridgecrest data preparation

The raw catalogue is stored at:

`data/raw/2026-07-23_scedc-scsn_socal_2016-2025_raw.txt`

Source and acquisition information is recorded in `data/raw/metadata.txt`.

Run:

`analyses/eda/eda_ridgecrest.qmd`

The EDA performs catalogue checks and filtering, defines the study window and magnitude threshold, and produces the analysis-ready catalogue:

`data/derived/ridgecrest_temporal_catalogue.rds`

The derived catalogue is already retained in Git, so this step can be skipped when reproducing Ridgecrest model-fitting analyses.

### 2. Simulation study

See `analyses/simulation/README.md` for the full dependency structure.

For complete reproduction of the simulation design, validation and catalogue generation, run the numbered scripts through `05_simulate_all.R` in order. `00_design.R` contains shared settings and is sourced automatically.

Because the final 300 synthetic catalogues are already retained, the main simulation analysis can instead begin with:

`analyses/simulation/06_fit_simulation.R`

This fits OU, MSE and RS to each of the 300 catalogues, producing **900 fitted models**.

After the production fits have been regenerated, scripts `07`–`10` reproduce the downstream convergence, recovery, misspecification and model-comparison analyses. `07a_misspecification_diagnostics.R` reproduces the additional convergence investigation and can be run separately.

### 3. Ridgecrest analysis

See `analyses/ridgecrest/README.md` for the full workflow.

Because the derived Ridgecrest catalogue is already retained, the empirical analysis can begin with:

`analyses/ridgecrest/01_fit_baseline.R`

This generates the baseline OU, MSE and RS fits.

The remaining numbered scripts reproduce the posterior parameter, temporal-function, model-comparison and sensitivity analyses.

The Ridgecrest analysis uses the shared priors and numerical settings defined by the simulation design, but reproducing the 900 simulation fits is **not** a prerequisite for fitting the Ridgecrest models.

### 4. Additional outputs

The standalone temporal-decay comparison figure is generated from:

`analyses/comparison-figure/`

The software environment can be recorded again by running:

`analyses/99_session_info.R`

## Computational requirements

Most preprocessing, plotting and summary stages are relatively lightweight. The main computational cost comes from fitting the temporal ETAS models.

| Stage                                   | Approximate requirement           |
|---------------------------------------|---------------------------------|
| Simulation binning check (`03`)         | \~15 min, \~150 MB                |
| Simulation tolerance check (`04`)       | \~15 min, \~90 MB                 |
| Production simulation fits (`06`)       | **\~20 h, \~10–12 GB, 6 workers** |
| Misspecification diagnostics (`07a`)    | **\~2 h, \~200 MB**               |
| Posterior recovery (`08`)               | **\~30–60 min**                   |
| Temporal misspecification (`09`)        | **\~30 min**                      |
| Ridgecrest baseline fits (`01`)         | 3 fits, \~30 MB                   |
| Ridgecrest magnitude sensitivity (`05`) | **\~20–30 min, \~60 MB**          |
| Ridgecrest prior sensitivity (`06`)     | 2 fits, \~20 MB                   |

Runtime estimates were recorded on the Apple M4 MacBook Air used for the dissertation analysis and will vary with hardware and system configuration.

## Important notes

### Fitted-model objects are not stored in Git

Individual fitted-model `.rds` objects are typically around 10 MB. The 900 production simulation fits alone require approximately 10–12 GB.

Large fit objects are therefore intentionally excluded from Git. The repository instead retains the input catalogues and compact manifests, summaries, diagnostics, tables and figures.

A fresh clone must regenerate the relevant fits before running downstream analyses that require the full fitted models. See `outputs/README.md` for details.

### Take care before running `06_fit_simulation.R`

`06_fit_simulation.R` is the most computationally demanding stage of the project. Ensure sufficient disk space and compute resources are available before running it.

Existing completed fits are reused, so an interrupted fitting run can be continued without intentionally refitting completed catalogue/model combinations.

### Run dependent scripts in order

Later scripts may require fitted objects or outputs created by earlier scripts. The exact dependencies are documented in the simulation and Ridgecrest READMEs.

### Numerical diagnostics

Some misspecification and sensitivity analyses investigate fits that were numerically difficult or did not satisfy the convergence criterion. These diagnostics form part of the original analysis workflow and should not be interpreted as indicating that every attempted fit converged.

### Software versions

The original software environment is recorded in `outputs/session_info.txt`.

## Modified `ETAS.inlabru`

The dissertation-specific package is included at:

`software/ETAS.inlabru/`

It extends the original OU temporal implementation to support MSE and RS decay while retaining OU as the default.

The project builds on the `ETAS.inlabru` package developed by Francesco Serafini, Mark Naylor, Finn Lindgren and Kirsty Bayliss, using `inlabru` and R-INLA for Bayesian inference.

See `software/ETAS.inlabru/README.md` for details of the package modifications, upstream authorship and the original project.
