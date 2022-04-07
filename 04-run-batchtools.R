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
reg_dir <- here::here("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir)

# Problems -----------------------------------------------------------
source("02-problems.R")
addProblem(name = "sim_a", fun = sim_a, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_b", fun = sim_b, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_c", fun = sim_c, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_d", fun = sim_d, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
source("03-algorithms.R")
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
    z_scale = c(1, 100),
    z_method = c("original"),
    theta = c("original", 1),
    t = c(1, 10, 50, 100),
    thresh = c(1e-3, 1e-7, 0)
  )
)

# theta == 1 only makes sense if we don't z_scale
algo_design$fwel_mt <- dplyr::filter(algo_design$fwel_mt, !(z_scale > 1 & theta == "1"))
# t == also makes more sense if we don't z_scale
algo_design$fwel_mt <- dplyr::filter(algo_design$fwel_mt, !(z_scale > 1 & t > 1))
# and only if we optimize theta
algo_design$fwel_mt <- dplyr::filter(algo_design$fwel_mt, !(theta != "original" & t > 1))
# threshold only for default z_scale and optimized theta
algo_design$fwel_mt <- dplyr::filter(algo_design$fwel_mt, !(theta != "original" & thresh != 1e-3))
algo_design$fwel_mt <- dplyr::filter(algo_design$fwel_mt, !(z_scale > 1 & thresh < 1e-3))

addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
unwrap(getJobPars(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------
if (interactive()) testJob(id = 5600)

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

