library(batchtools)
library(data.table)
source(here::here("1-proof-of-concept", "1-problems.R"))
source(here::here("1-proof-of-concept", "1-algorithms.R"))
invisible(lapply(list.files(here::here("R"), pattern = "*.R", full.names = TRUE), source, echo = FALSE))

if (!dir.exists(here::here("results"))) dir.create(here::here("results"))

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_n = 1000,
  sim_cache = TRUE,
  repls = 1000
)

set.seed(config$global_seed)

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "poc"
reg_dir <- here::here("registries", reg_name)
if (file.exists(reg_dir)) {
  message("Deleting old registry")
  unlink(reg_dir, recursive = TRUE)
  stopifnot("Deleting registry failed" = !file.exists(reg_dir))
}
reg <- makeExperimentRegistry(
  file.dir = reg_dir,
  source = c(
    list.files(here::here("R"), pattern = "*.R", full.names = TRUE),
    here::here("1-proof-of-concept", "1-problems.R")
  ),
  seed = config$global_seed
)

# Problems -----------------------------------------------------------
addProblem(name = "sim_a", fun = sim_a, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_b", fun = sim_b, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_c", fun = sim_c, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "sim_d", fun = sim_d, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
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
    t = 100,
    thresh = 1e-7
  )
)

addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobTable())

# Test jobs -----------------------------------------------------------
if (FALSE) testJob(id = 1600)

# Submit -----------------------------------------------------------
if (grepl("blog\\d{1}", Sys.info()[["nodename"]])) {
  ids <- findNotSubmitted()
  ids[, chunk := chunk(job.id, chunk.size = 20)]

  submitJobs(
    ids = ids,
    resources = list(
      partition = "teton", memory = 1024, walltime = 2 * 3600,
      comment = "cooper-poc",
      chunks.as.arrayjobs = FALSE
    )
  )

  waitForJobs()

} else {
  # Otherwise, run a subset of jobs locally, 5 randomly sampled per simulation setting
  sample_ids = jobtbl[, .SD[sample(nrow(.SD), 5)], by = "problem"]

  message("Submitting subset of jobs only!")
  submitJobs(sample_ids)

  # Wait for jobs to complete before proceeding
  waitForJobs()
}

# Checking for errors
if (nrow(findErrors()) > 0) {
  warning("Errors found!")
  getErrorMessages()
} else {
  message("No errors found!")
}

# Get results -------------------------------------------------------------
res <-  reduceResultsDataTable()
pars <- unwrap(getJobPars())

jobs_done <- 100 * nrow(res)/nrow(getJobTable())
jobs_err <- 100 * nrow(findErrors())/nrow(getJobTable())

if (jobs_err == 1) {
  stop("All jobs failed!")
}

res_long <- rbindlist(lapply(res$job.id, \(id) {
  result <- res[(job.id == id), result][[1]]
  result
  result[, job.id := ..id]
}))

res_long <- merge(pars, res_long, by = "job.id")

cli::cli_alert_info("Saving {nrow(res_long)} results from {jobs_done}% of jobs with {jobs_err}% errors")

saveRDS(res, here::here("results", "1-results-poc.rds"))
saveRDS(res_long, here::here("results", "1-results-long-poc.rds"))
