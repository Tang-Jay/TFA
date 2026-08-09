# EL-Based Tolerance Testing for Algorithm Fairness

This repository contains the R code and data used to study empirical-likelihood (EL) based tolerance tests for algorithmic fairness. It includes theoretical rejection-probability calculations, simulation studies under several data-generating processes, and a COMPAS case study comparing EL, constrained EL (CEL), split EL (SEL), and adjusted split EL (ASEL) confidence intervals.

## Repository structure

```text
.
├── Fig1.R                 # Theoretical local rejection-probability curves
├── Fig2-5.R               # Power plots from saved simulation results
├── Fig6-Data.R            # COMPAS audit: African-American subgroups
├── Fig7-Data.R            # COMPAS audit: sex-by-age subgroups
├── Fig8-Data.R            # COMPAS audit: female race-by-age subgroups
├── Fig6-8.R               # Confidence-interval plots for the COMPAS audits
├── Tab1.R                 # Table 1 from the saved audit results
├── SEL-{E,N,T,C}-data.R   # Simulation scripts for four distributions
├── Functions/             # EL, CEL, SEL, and ASEL confidence-interval routines
└── data/                  # COMPAS data and precomputed simulation/audit results
```

The simulation suffixes identify the data-generating process:

- `E`: exponential
- `N`: normal
- `T`: Student's *t* with 3 degrees of freedom
- `C`: chi-square configuration (with a normal fallback for `test2`/`test3`)

## Requirements

- R
- R packages: `emplik`, `ggplot2`, `dplyr`, `readr`, `Cairo`, `scales`, and `tibble`

Install the required packages from CRAN:

```r
install.packages(c(
  "emplik", "ggplot2", "dplyr", "readr",
  "Cairo", "scales", "tibble"
))
```

## Quick start

Run all commands from the repository root because the scripts use relative paths.

The repository already includes precomputed `.RData` and `.csv` files, so the main results can be inspected without rerunning the simulations:

```sh
Rscript Fig1.R
Rscript Fig2-5.R
Rscript Fig6-8.R
Rscript Tab1.R
```

By default, most plotting scripts display a plot without writing a file. Set `save_plot <- TRUE` near the top of the relevant script to save the figure. Similarly, set `save_data <- TRUE` in a data-generation script when regenerated results should be written to `data/`.

## Reproducing the simulations

Each simulation script exposes its main settings near the bottom of the file, including the test, sample sizes, number of Monte Carlo repetitions, and output switches. After choosing the settings, run the required scripts:

```sh
Rscript SEL-E-data.R
Rscript SEL-N-data.R
Rscript SEL-T-data.R
Rscript SEL-C-data.R
Rscript Fig2-5.R
```

The supported test identifiers are:

- `test1` / `greater`: upper one-sided alternative
- `test2` / `less`: lower one-sided alternative
- `test3` / `or`: two-sided alternative outside a tolerance interval

Simulation can be computationally expensive. The checked-in `.RData` files can be used directly by `Fig2-5.R`; make sure its `test`, `nsims`, `ns`, and `a2` settings match the filenames available in `data/`.

## Reproducing the COMPAS analysis

The COMPAS workflow has two stages:

1. Run `Fig6-Data.R`, `Fig7-Data.R`, or `Fig8-Data.R` to compute audit confidence intervals. Enable `save_data` to write the results.
2. Set the desired `alpha` and `test` in `Fig6-8.R`, then run it to draw the corresponding confidence-interval figures.

For example:

```sh
Rscript Fig7-Data.R
Rscript Fig6-8.R
Rscript Tab1.R
```

The audit scripts compare four methods:

- `ELCI`: empirical-likelihood confidence interval
- `CELCI`: constrained empirical-likelihood confidence interval
- `SELCI`: split empirical-likelihood confidence interval
- `ASELCI`: adjusted split empirical-likelihood confidence interval

The input dataset is `data/compas_data.csv`. Generated audit summaries follow the naming convention
`<group>-Fairness-Audit-<confidence>-<test>-<method>.csv`; the `ELCI` filenames omit the test because that interval is two-sided.

## Notes

- Open `SEL.Rproj` in RStudio for a project configured at the repository root.
- Generated figures are EPS by default; some simulation scripts can also save PNG files.
- Existing files in `data/` are research artifacts and may reflect parameter settings different from a script's current defaults. Check the filename and the parameter block before comparing or overwriting results.
