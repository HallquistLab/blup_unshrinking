# REPAIR: Random-effect prediction adjustment for inference in regression

This release contains the **simulation pipeline for Studies 1--4** used in "REPAIR: Random-effect prediction adjustment for inference in regression" by Zachary A. Vig and Michael N. Hallquist.

The simulation data are archived separately in Zenodo record [10.5281/zenodo.23965547](https://doi.org/10.5281/zenodo.23965547). See [`data/README.md`](data/README.md) for file checksums and loading instructions.

## Running the simulations

From this directory, install the R package dependencies used by the shared helpers, then run:

```sh
Rscript vig_hallquist_2026/random_effects_structural_simulation.R
```

The default selection is `allv2`. Individual studies can be selected with `1v2`, `2v2`, `3v2`, or `4v2`.
