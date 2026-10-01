# Gaussian REPAIR pathways

Executed brms, rstanarm, Mplus Bayes random-slope, and Mplus ML two-factor examples feed the **existing**, unchanged `R/stage2_estimators.R` Fuller implementation.

Start with `report.html` (self-contained rendered report) or `report.qmd` (source).

From the parent repository:

```sh
Rscript gaussian_repair_pathways/run_models.R all
Rscript gaussian_repair_pathways/analyse.R 250
Rscript gaussian_repair_pathways/tests/check_outputs.R
quarto render gaussian_repair_pathways/report.qmd
```

Required: R with brms, cmdstanr + CmdStan, rstanarm, rstan, posterior, MplusAutomation, dplyr, purrr, tibble, geigen, ggplot2, tidyr, knitr, rmarkdown; working licensed Mplus; Quarto.

`run_models.R data` prepares deterministic demo datasets. Subsequent `brms`, `rstanarm`, and `mplus` phases can run separately. The render step only reads completed results. The brms fit is cached; rstanarm and Mplus phases refit when invoked.

- `R/helpers.R`: direct Gaussian and conditional-prior-removal adapters, simulation, and calls to existing Fuller code.
- `run_models.R`: actual software fits and extraction.
- `analyse.R`: pathway identities, exported Mplus score check, nondiagonal-R check, zero-ME limit, conditional Monte Carlo study.
- `mplus/`: inputs, measurements, and original outputs; the report uses real Mplus fits.
- `results/`: portable plug-in estimates and outputs, all MC replicates, source hashes, diagnostics, and sessions.
- `cache/`: large posterior fit objects, ignored by Git.

The repeated-cohort study is conditional: it does not refit the three Bayesian packages each replication. It re-estimates pooled residual variance for stress, holds the other shared parameters at truth, and holds the factor measurement parameters at truth. It validates the ME/Fuller interface, not unconditional first-stage uncertainty propagation. Actual fitted-model examples establish the software pathways separately.

Mplus printed estimates are rounded. The report distinguishes machine-precision equivalence within a consistent parameter set from the small rounding discrepancies when unshrinking actual full-precision saved factor scores using reconstructed covariances.
