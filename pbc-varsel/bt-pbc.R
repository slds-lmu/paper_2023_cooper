# source(here::here("2-variable-selection-sim/2-simulation.R"))
# source(here::here("2-variable-selection-sim/2-algorithms.R"))

library(batchtools)
# library(randomForestSRC)
# library(CoxBoost)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  repls = 10
)

set.seed(config$global_seed)
# Set to FALSE to remove + recreate registry
continue_bt <- FALSE

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "pbc_varsel"
reg_dir <- here::here("registries", reg_name)

if (continue_bt) {
  loadRegistry(reg_dir, writeable = TRUE)
} else {
  unlink(reg_dir, recursive = TRUE)
  makeExperimentRegistry(
    file.dir = reg_dir,
    # packages = c("randomForestSRC", "CoxBoost"),
    seed = config$global.seed,
    source = c(here::here("pbc-varsel", "pbc-varsel.R"))
  )
}

# Problems -----------------------------------------------------------
addProblem(name = "get_pbc_data", fun = get_pbc_data, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "pbc_varsel_score", fun = pbc_varsel_score)
# addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_wrapper)
# addAlgorithm(name = "coxboost", fun = coxboost_varselect_wrapper)


# Experiments -----------------------------------------------------------
prob_design <- list(
  get_pbc_data = expand.grid(
    num_noise = sort(c(0, as.vector(sapply(c(1, 2, 5), \(x) x * 10^(0:2)))))
  )
)

algo_design <- list(
  pbc_varsel_score = expand.grid(
    mt_max_iter = 3,
    alpha = c(1),
    t = c(100),
    thresh = c(1e-7),
    conf_level = c(0.9, 0.95, 0.99)
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobTable(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------
if (FALSE) testJob(id = 1)

# Submit -----------------------------------------------------------

ids <- findNotStarted()
submitJobs(ids)
waitForJobs()

res <- ijoin(tidyr::unnest(reduceResultsDataTable(), cols = "result"), jobtbl)
saveRDS(res, here::here("pbc-varsel", "results.rds"))
