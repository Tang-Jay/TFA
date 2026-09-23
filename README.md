# Tolerance-Based Fairness Auditing

This repository contains the R code and data used to study empirical-likelihood (EL) based tolerance tests for algorithmic fairness. It includes theoretical rejection-probability calculations, single-test power simulations, multiple-testing simulations with Benjamini–Hochberg (BH) adjustment, and a COMPAS case study comparing EL, constrained EL (CEL), split EL (SEL), and adjusted split EL (ASEL) confidence intervals.

## Repository structure

```text
.
├── Fig1.R                 # Theoretical local rejection-probability curves
├── Fig2-4.R               # Power plots from saved single-test simulation results
├── Fig5-Data.R            # Multiple-testing simulations with BH adjustment
├── Fig5.R                 # FFR–Power plots from multiple-testing CSV results
├── Fig6-Data.R            # COMPAS audit: African-American subgroups
├── Fig7-Data.R            # COMPAS audit: sex-by-age subgroups
├── Fig8-Data.R            # COMPAS audit: female race-by-age subgroups
├── Figs6-8.R              # Confidence-interval plots for the COMPAS audits
├── Tabs1-2.R              # FFR and power tables from multiple-testing results
├── Tab3.R                 # African-American All-group confidence bounds
├── SEL-{E,N,T,C}-data.R    # Single-test simulation scripts for four distributions
├── Functions/             # EL, CEL, SEL, and ASEL confidence-interval routines
├── TFA.Rproj              # RStudio project
└── data/                  # COMPAS data and precomputed simulation/audit results
```

The suffixes in `SEL-{E,N,T,C}-data.R` identify the data-generating process:

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

## Figure and table map

| Output | Data source or generation script | Plotting or table script |
| --- | --- | --- |
| Figure 1 | Theoretical calculations | `Fig1.R` |
| Figures 2–4 | `SEL-{E,N,T,C}-data.R` or included `.RData` files | `Fig2-4.R` |
| Figure 5 | `Fig5-Data.R` | `Fig5.R` |
| Tables 1–2 | `Fig5-Data.R` | `Tabs1-2.R` |
| Figures 6–8 | `Fig6-Data.R`, `Fig7-Data.R`, `Fig8-Data.R`, or included audit CSV files | `Figs6-8.R` |
| Table 3 | African-American audit CSV files | `Tab3.R` |

## Reproducing the single-test simulations: Figures 2–4

Each `SEL-*-data.R` script exposes its settings near the bottom, including sample sizes, Monte Carlo repetitions, test selection, and output switches. Before regenerating results for `Fig2-4.R`:

1. Set `save_data <- TRUE` and `power_dir <- "data"` in the required generation scripts. Their current output directory is `"."`, whereas `Fig2-4.R` reads from `data/`.
2. Choose `ns`, `nsim`, and `a2`, and match them to `ns`, `nsims`, and `a2` in `Fig2-4.R`. The current generator defaults differ across distributions and do not regenerate every included file.
3. Select the same test in the generator and plotter. Passing a test as the first command-line argument selects that test only, as in this example:

```sh
Rscript SEL-E-data.R test1
Rscript SEL-N-data.R test1
Rscript SEL-T-data.R test1
Rscript SEL-C-data.R test1
Rscript Fig2-4.R
```

The supported test identifiers are:

- `test1` / `greater`: upper one-sided alternative
- `test2` / `less`: lower one-sided alternative
- `test3` / `or`: two-sided alternative outside a tolerance interval

Without command-line arguments, each generator uses its own `test` and `run_all_tests` settings. When `run_all_tests` is `TRUE`, it runs all three tests.

Simulation can be computationally expensive. To use the included `.RData` files directly, match `sim_prefixes`, `test`, `nsims`, `ns`, and `a2` in `Fig2-4.R` to the available filenames. The plotter warns and skips missing files. Results use the naming convention:

```text
SEL-<distribution>-<test>-nsim-<repetitions>-a2-<four-decimal-value>-n-<sample-size>.RData
```

## Reproducing the multiple-testing simulations: Figure 5 and Tables 1–2

`Fig5-Data.R` simulates multiple group tests and, by default, compares CEL-BH, SEL-BH, ASEL-BH, and t-BH. `Fig5.R` plots FFR and power; `Tabs1-2.R` prints the global-null FFR table and selected all-alternative power results. The CSV column `fdr` stores the mean of replicate-level false discovery proportions and is labeled FFR in these plotting and table scripts.

First align the data paths:

1. In `Fig5-Data.R`, set `plot_dir <- "data"` and keep `save_data <- TRUE`. The current `plot_dir <- "."` writes CSV files to the working directory.
2. In both `Fig5.R` and `Tabs1-2.R`, replace the personal `project_dir <- path.expand(...)` setting with `project_dir <- "."`. Their existing `data_dir <- file.path(project_dir, "data")` will then read from the repository's `data/` directory. `Fig5.R` also saves its EPS files there by default.

Check the parameter sections before running. The current generation defaults are `n = 1000`, `ks = c(26)`, Normal and Exponential distributions, `nsim = 100`, `qs = 0.05`, `a2 = 0.08`, `epsilon0 = 0.25`, and `seed = 2`. Choose the repetition count needed for your analysis; match the sample size, group count, methods, and effect/null grids in the readers to the generated data.

After applying those settings, run:

```sh
Rscript Fig5-Data.R
Rscript Fig5.R
Rscript Tabs1-2.R
```

The generator writes one CSV per method, distribution, sample size, group count, and alternative proportion, combining the effect or null grid within each file:

```text
<method>-<distribution>-n<sample-size>-K<group-count>-pi<alternative-proportion>.csv
# Example: CEL-N-n1000-K26-pi0.csv
```

For the default readers, distribution codes are `N` (Normal) and `E` (Exponential). Figure 5 requires proportions `0`, `0.5`, and `1`; Tables 1–2 require `0` and `1`. The readers check required columns, scenario settings, and batch consistency. Use files from a consistent simulation configuration. Same-name files are overwritten by the generator.

## Reproducing the COMPAS analysis

The input dataset is `data/compas_data.csv`. To regenerate the audit results:

1. Set `alpha`, `test`, and `save_data <- TRUE` in the required audit script, then run it to write CSV results to `data/`.
2. Set the matching `alpha` and `test` in `Figs6-8.R`, then run it to draw the confidence intervals. Enable `save_plot` for EPS output.

| Audit | Generation script | Plotter selection |
| --- | --- | --- |
| African-American subgroups | `Fig6-Data.R` | `test <- "test1"` |
| Sex-by-age subgroups | `Fig7-Data.R` | `test <- "test1"` |
| Female race-by-age subgroups | `Fig8-Data.R` | `test <- "test2"` |

Here `alpha` denotes the **confidence level**, either `0.90` or `0.95`. In contrast, `alpha` in `Fig1.R` denotes a significance level, currently `0.05`. Match confidence levels explicitly: `Fig6-Data.R` and `Fig7-Data.R` currently default to `0.95`, while `Fig8-Data.R` and `Figs6-8.R` default to `0.90`.

For example, after selecting the same confidence level in `Fig7-Data.R` and `Figs6-8.R` and enabling data saving:

```sh
Rscript Fig7-Data.R
Rscript Figs6-8.R
```

The audit scripts compare four methods:

- `ELCI`: empirical-likelihood confidence interval
- `CELCI`: constrained empirical-likelihood confidence interval
- `SELCI`: split empirical-likelihood confidence interval
- `ASELCI`: adjusted split empirical-likelihood confidence interval

Generated audit summaries follow the naming convention
`<group>-Fairness-Audit-<confidence>-<test>-<method>.csv`; the `ELCI` filenames omit the test because that interval is two-sided.

`Tab3.R` reads the African-American `All`-group lower bounds at both 90% and 95% confidence levels and prints Table 3. These files are already included. To regenerate that table's inputs, run `Fig6-Data.R` separately at each confidence level with `test <- "test1"` and `save_data <- TRUE`, then run:

```sh
Rscript Tab3.R
```
