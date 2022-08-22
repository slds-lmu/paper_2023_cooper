library(batchtools)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 100
)

set.seed(config$global_seed)

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "fwel_sim_variableselect"
reg_dir <- here::here("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir)

# Problems -----------------------------------------------------------
addProblem(name = "sim_surv_binder", fun = sim_surv_binder, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "fwel_mt", fun = fwel_mt_varselect_wrapper)


# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_surv_binder = expand.grid(n = 400, p = 5000, ce = 0.5, lambda = 0.1, lambda_c = 0.1)
)

algo_design <- list(
  fwel_mt = expand.grid(
    mt_max_iter = 3,
    alpha = c(1),
    t = c(100),
    thresh = c(1e-3, 1e-7)
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
unwrap(getJobPars(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------
if (interactive()) testJob(id = 4)

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

