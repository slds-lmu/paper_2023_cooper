# Recreate all results
renv_state <- renv::status()

if (isFALSE(renv_state$synchronized)) {
  warning("Run renv::restore() to ensure fully reproducible results!")
}

# Run batchtools -- will need manual setup / an HPC environment or similar!
source(here::here("1-proof-of-concept/1-run-batchtools.R"))
source(here::here("2-variable-selection-sim/2-run-batchtools.R"))

# Produce plots and tables from simulation results
source(here::here("1-proof-of-concept/1-analysis.R"))
source(here::here("2-variable-selection-sim/2-analysis.R"))

# Application example: Fit models for prediction
source(here::here("3-example-bladder/3-example-bladder-pred.R"))
# Produce plots
source(here::here("3-example-bladder/3-example-bladder-pred-analysis.R"))
# Fit models for variable selection and print selected variables etc.
source(here::here("3-example-bladder/3-example-bladder-varsel.R"))
