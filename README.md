# REPAIR: Random-effect prediction adjustment for inference in regression

This release contains the sanitized code and documentation for the REPAIR
simulation study by Zach Vig and Michael Hallquist (2026). It contains the
study-specific simulation drivers under `vig_hallquist_2026/` and the shared R
helpers required to run them under `R/`.

## Contents

- `vig_hallquist_2026/`: REPAIR study designs, simulation runner, analyses, and
  batch scripts.
- `R/`: shared helper modules sourced by the REPAIR runner.
- `data/README.md`: access instructions and checksums for the separately
  archived simulation data.

The large Parquet simulation datasets are available from the data record at
[10.5281/zenodo.23965547](https://doi.org/10.5281/zenodo.23965547).
The software release is archived at
[10.5281/zenodo.23982336](https://doi.org/10.5281/zenodo.23982336); this
identifier refers to the earlier full-repository archive. The DOI for this
sanitized release will be added after Zenodo processes the release.

See `vig_hallquist_2026/README.Rmd` for study details and execution guidance.
