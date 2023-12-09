# source(here::here("2-variable-selection-sim/2-simulation.R"))
# source(here::here("2-variable-selection-sim/2-algorithms.R"))

library(batchtools)
options(batchtools.progress = TRUE)


# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  repls = 100
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
addAlgorithm(name = "cooper_varsel_pbc", fun = cooper_varsel_pbc)
addAlgorithm(name = "rfsrc_pbc_varsel", fun = rfsrc_pbc_varsel)
addAlgorithm(name = "coxboost_pbc_varsel", fun = coxboost_pbc_varsel)


# Experiments -----------------------------------------------------------
prob_design <- list(
  get_pbc_data = expand.grid(
    conf_level = c(0.9, 0.95, 0.99),
    num_noise = sort(c(0, as.vector(sapply(c(1, 2, 5), \(x) x * 10^(1:3)))))
  )
)

algo_design <- list(
  cooper_varsel_pbc = expand.grid(
    mt_max_iter = 3,
    alpha = c(1),
    t = c(100),
    thresh = c(1e-7)
  ),
  rfsrc_pbc_varsel = expand.grid(
    nodesize = 15,
    splitrule = "logrank",
    importance = "random",
    cutoff_method = "vita"
  ),
  coxboost_pbc_varsel = expand.grid(
    stepno = 200,
    penalty = 1000
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobTable(), c("algo.pars", "prob.pars"))
jobtbl[, chunk := lpt(log10(num_noise + 1), n.chunks = 1000)]

jobtbl[, (n = .N), by = .(chunk)]




# Test jobs -----------------------------------------------------------
if (FALSE) testJob(id = 250)

if (FALSE) {
  jobtbl[findNotSubmitted()] |>
    submitJobs()

  jobtbl |>
    dplyr::filter(num_noise == 10) |>
    dplyr::group_by(algorithm, conf_level) |>
    dplyr::slice_sample(n = 1) |>
    submitJobs()
}

# Submit -----------------------------------------------------------

#submitJobs(jobtbl)
#waitForJobs()

#message("Done!")

# res_file <- here::here("pbc-varsel", "results.rds")
# if (file.exists(res_file)) file.remove(res_file)
# res <- ijoin(tidyr::unnest(reduceResultsDataTable(), cols = "result"), jobtbl)
# message("Writing ", res_file)
# saveRDS(res, res_file)
