# Simulation Code for Cooperative Penalized Regression (CooPeR) Project

Simulation code for [`cooper`](https://github.com/jemus42/cooper) in competing
risk settings, based on the [original fwelnet implementation](https://github.com/kjytay/fwelnet/)
([preprint](https://arxiv.org/pdf/2006.01395.pdf)).

## Structure

- `R/`: Helper functions for simulation and plotting. 
   - Loaded automatically via `.Rprofile`.
- `n-<name>/` folders with code for specific settings

Note the project uses `renv` to ensure consistent R dependencies.  
Upon first load it's going to bug you to run `renv::restore()` to locally 
install dependencies as specified in `renv.lock`.

Run `renv::update()` / `renv::snapshot()` to update & lock dependencies.
Run `renv::restore()` to install all dependencies.
Please note that you may want to install the appropriate R version (see e.g. [rig](https://github.com/r-lib/rig))
