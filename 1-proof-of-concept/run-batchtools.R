library(batchtools)

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
reg_name <- "fwel_simulations"
reg_name <- "DEBUG"
reg_dir <- here::here("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, source = here::here("1-proof-of-concept", "sim01-sim_cr.R"))

# Problems -----------------------------------------------------------
source(here::here("1-proof-of-concept", "sim01-problems.R"))
addProblem(name = "sim_a", fun = sim_a, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_b", fun = sim_b, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_c", fun = sim_c, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_d", fun = sim_d, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
source(here::here("1-proof-of-concept", "sim01-algorithms.R"))
addAlgorithm(name = "fwel_mt", fun = fwel_mt_wrapper)


# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_a = expand.grid(n = config$sim_n),
  sim_b = expand.grid(n = config$sim_n),
  sim_c = expand.grid(n = config$sim_n),
  sim_d = expand.grid(n = config$sim_n)
)

algo_design <- list(
  fwel_mt = expand.grid(
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
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 50)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = reg_name, chunks.as.arrayjobs = TRUE,
                              ncpus = 1, memory = 6000, walltime = 10*24*3600,
                              max.concurrent.jobs = 40))
} else {
  ids <- findNotStarted()
  submitJobs(ids = ids)
}
waitForJobs()


# Monitor jobs ------------------------------------------------------------
if (interactive()) {
  getStatus()

  ijoin(
    unwrap(getJobPars(findErrors()), c("algo.pars", "prob.pars")),
    getErrorMessages(findErrors())
  )
}

# Get results -------------------------------------------------------------
# res <-  ijoin(reduceResultsDataTable(), flatten(getJobPars()))
# res

