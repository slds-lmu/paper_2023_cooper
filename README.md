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
   - `3-example-bladder`: Supports section 4. Application Example.
- `data-raw/`: Preprocessing script and intermediate / raw files for the bladder cancer example application.
- `data/`: Processed data in different forms depending on use-case.
- `run-all.R`: Wraps relevant scripts needed to reproduce results, with additional instructions for simulations. See below.
- `registries/`: Holds `batchtools` registries as created by `n-run-batchtools.R` scripts.
- `renv/`, `renv.lock`, `.renvignore` ensure fixed R package versions via `renv`.

Note the project uses `renv` to ensure consistent R dependencies.  
Upon first load it's going to bug you to run `renv::restore()` to locally 
install dependencies as specified in `renv.lock`.

## Reproducing results 

- Run `renv::restore()` to install all dependencies.
   Please note that you may want to install the appropriate R version (see e.g. [rig](https://github.com/r-lib/rig))
- Refer to `run-all.R` to reproduce figures and tables from simulations results and data example.

### Section 2: Proof of Concept

- To re-run simulations, the `1-run-batchtools.R` scripts will need to be adapted to the local
environment as they expect an HPC setting.
- Lines 129f will need to be adapted to a local environment.
- The number of replications can be changed using the `repls` variable in line 13.
- Results are stored in `results/` as .rds files `1-results-poc.rds` and `1-results-long-poc.rds`.
- Figure 1 of the manuscript is then produced by `1-figure-1-proof-of-concept-bias.R`.
- Figures 2 through 4 are produced by `2-figure-2-3-4-high-dim-sim-variable-selection.R`.
- Figure 5 is produced by `2-figure-5-high-dim-sim-performance.R`.
- Tables 1 through 4 are produced by `2-tables-1-2-3-4-high-dim-sim-variable-selection.R`.

### Section 3: High-Dimensional Data

- A batchtools script analogous to `1-run-batchtools.R`
- Lines 106f will need to be adapted to a local environment.
- The number of replications can be changed using the `repls` variable in line 15.
- Results are stored in `results/` as .rds files
   - `2-results-varsel-csc-varsel.rds` for the variable selection performance scores.
   - `2-results-varsel-csc-perf.rds` for the prediction performance scores.


### Section 4. Bladder Cancer Data

- The data is described by Dyrskjøt et al. (2007).
- Visit https://aacrjournals.org/clincancerres/article/13/12/3545/13137/Gene-Expression-Signatures-Predict-Outcome-in-Non
   - Download supplemental tables 1 and to to `data-raw`.
   - The script `data-raw/preprocess-bladder-data.R` downloads the remaining data automatically given an internet connection is available

The expected folder structure is then

```
data-raw
├── 10780432ccr062940-sup-supplemental_file_1.xls
├── 10780432ccr062940-sup-supplemental_file_2.xls
├── GSE5479_Final_processed_data_1.txt
├── GSE5479_Final_processed_data_2.txt
└── preprocess-bladder-data.R
```

- The remaining scripts in `3-example-bladder/` will take some time to fit models but produce results prefixed with `3-` in `results/`
- Figure 6 is produced by `3-figure-6-bladder-performance.R`
