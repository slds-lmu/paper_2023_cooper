# Simulation Code for Cooperative Penalized Regression (CooPeR) Project

Simulation code for [`cooper`](https://github.com/jemus42/cooper) in competing
risk settings, based on the [original fwelnet implementation](https://github.com/kjytay/fwelnet/)
([preprint](https://arxiv.org/pdf/2006.01395.pdf)).

## Structure

- `R/`: Helper functions for simulation and plotting. 
   - Loaded automatically via `.Rprofile`.
- `n-<name>/`: Code for simulation experiments and application example.
   - `1-proof-of-concept`: Supports section 3.1. Proof of Concept.
   - `2-variable-selection`: Supports section 3.2. High-Dimensional Data.
   - `3-example-bladder`: Supprots section 4. Application Example.
- `data-raw/`: Preprocessing script and intermediate / raw files for the bladder cancer example application.
- `data/`: Processed data in different forms depending on use-case.
- `run-all.R`: Wraps relevant scripts needed to reproduce results, with additional instructions for simulations. See below.
- `registries/`: Holds `batchtools` registries as created by `n-run-batchtools.R` scripts.
- `renv/`, `renv.lock`, `.renvignore` ensure fixed R package versions via `renv`.

Note the project uses `renv` to ensure consistent R dependencies.  
Upon first load it's going to bug you to run `renv::restore()` to locally 
install dependencies as specified in `renv.lock`.

## Reproducing results 

1. Run `renv::restore()` to install all dependencies.
   Please note that you may want to install the appropriate R version (see e.g. [rig](https://github.com/r-lib/rig))
2. Refer to `run-all.R` to reproduce figures and tables from simulations results and data example.
   2.1 To re-run simulations, the `n-run-batchtools.R` scripts will need to be adapted to the local environment as they expect an HPC setting
      Results are stored in `results/` as .rds files and are used to create plots and tables.
   2.2 The remaining scripts can be run locally, yet scripts in `3-example-bladder/` will take some time to fit models.