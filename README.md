# fwelnet_cr_sim

Simulation code for [`fwelnet`](https://github.com/jemus42/fwelnet) in competing
risk settings, based on the [original fwelnet implementation](https://github.com/kjytay/fwelnet/)
([preprint](https://arxiv.org/pdf/2006.01395.pdf)).

## Structure

- `R/`: Helper functions for simulation and plotting. Load via `devtools::load_all` 
    (hence the dummy `DESCRIPTION`) file.
- `*.Rmd`: Quick / interactive experiments to get a feeling for what's worth including in a simulation

Note the project uses `renv` to ensure consistent R dependencies.  
Upon first load it's likely going to bug you to run `renv::restore()` to locally 
install dependencies as specified in `renv.lock`.

Run `renv::update()` / `renv::snapshot()` to update & lock dependencies.
