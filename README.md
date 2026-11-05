# REPAIR: Random-effect prediction adjustment for inference in regression

This release contains the amended **v2 simulation pipeline for Studies 1--4** used in the REPAIR manuscript (Vig Hallquist, 2026). It is a sanitized release tree: legacy v1 designs, Study 5, ICC bridge code, variance-ratio sensitivity code, and manuscript analysis scripts are not included.

The simulation data are archived separately in Zenodo record [10.5281/zenodo.23965547](https://doi.org/10.5281/zenodo.23965547). See [`data/README.md`](data/README.md) for file checksums and loading instructions.

## Running the simulations

From this directory, install the R package dependencies used by the shared helpers, then run:

```sh
Rscript vig_hallquist_2026/random_effects_structural_simulation.R
```

The default selection is `allv2`. Individual amended studies can be selected with `1v2`, `2v2`, `3v2`, or `4v2`.
