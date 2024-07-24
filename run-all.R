# Recreate all results
renv_state <- renv::status()

# Number of requested cores for later model fits. Does not affect simulation
n_cores <- 20

if (isFALSE(renv_state$synchronized)) {
  warning("Run renv::restore() to ensure fully reproducible results!")
}

# Run batchtools -- will need manual setup / an HPC environment or similar!
if (FALSE) {
  source(here::here("1-proof-of-concept/1-run-batchtools.R"))
  source(here::here("2-variable-selection-sim/2-run-batchtools.R"))
}


# Producing figures and tables ---------------------------------------------------------------

# Proof of Concept: Figure 1
cli::cli_alert_info("Recreating results for 1-proof-of-concept")
source(here::here("1-figure-1-proof-of-concept-bias.R"))

# High-Dimensional data simulation
cli::cli_alert_info("Recreating results for 2-variable-selection-sim")
source(here::here("2-figure-2-3-4-high-dim-sim-variable-selection.R"))
source(here::here("2-tables-1-2-3-4-high-dim-sim-variable-selection.R"))
source(here::here("2-figure-5-high-dim-sim-performance.R"))


# Application example: Bladder data -----------------------------------------------------------

# Fit models for prediction
# If more than n-1 available cores are requested we tone it down a notch.
n_cores <- min(n_cores, parallel::detectCores() - 1)

cli::cli_alert_info("Trying to run subsequent model fits using {n_cores} parallel threads (rfsrc, CoxBoost)")
options(rf.cores = n_cores, mc.cores = n_cores)

cli::cli_alert_info("Running bladder cancer prediction models")
source(here::here("3-example-bladder/3-example-bladder-pred.R"))

# Fit models for variable selection
cli::cli_alert_info("Running bladder cancer variable selection models")
source(here::here("3-example-bladder/3-example-bladder-varsel.R"))

# Plot performance for prediction
cli::cli_alert_info("Producing bladder cancer prediction performance plots")
source(here::here("3-figure-6-bladder-performance.R"))
