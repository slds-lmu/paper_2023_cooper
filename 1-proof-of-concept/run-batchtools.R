library(batchtools)
source(here::here("1-proof-of-concept", "sim01-problems.R"))
source(here::here("R", "utils.R"))

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_n = 1000,
  sim_cache = TRUE,
  repls = 100
)

set.seed(config$global_seed)

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "poc"
reg_dir <- here::here("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(
  file.dir = reg_dir,
  source = here::here("R", "utils.R"),
  seed = config$global_seed
)

# Problems -----------------------------------------------------------
addProblem(name = "sim_a", fun = sim_a, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_b", fun = sim_b, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_c", fun = sim_c, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_d", fun = sim_d, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------

cooper_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 2,
    t = 1, a = 0.5, thresh = 1e-3
) {

  cooper::cooper(
    instance$data,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh
  )
}


addAlgorithm(name = "cooper", fun = cooper_wrapper)

# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_a = expand.grid(n = config$sim_n),
  sim_b = expand.grid(n = config$sim_n),
  sim_c = expand.grid(n = config$sim_n),
  sim_d = expand.grid(n = config$sim_n)
)

algo_design <- list(
  cooper = expand.grid(
    mt_max_iter = 5,
    alpha = 1,
    t = c(1, 50, 100),
    thresh = c(1e-3, 1e-7, 0)
  )
)

addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
unwrap(getJobPars(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------
if (interactive()) testJob(id = 1600)

# Submit -----------------------------------------------------------
if (grepl("blog\\d{1}", Sys.info()[["nodename"]])) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(
    ids = ids,
    resources = list(memory = 2048, walltime = 6*3600)
  )
} else {
  ids <- findNotStarted()
  submitJobs(ids = ids)
}
waitForJobs()

# Get results -------------------------------------------------------------
res <-  ijoin(reduceResultsDataTable(), flatten(getJobPars()))
res

saveRDS(res, here::here("results", "results-poc.rds"))
