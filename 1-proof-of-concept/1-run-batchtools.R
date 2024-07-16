library(batchtools)
library(data.table)
source(here::here("1-proof-of-concept", "1-problems.R"))
source(here::here("R", "utils.R"))
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
makeExperimentRegistry(
  file.dir = reg_dir,
  source = c(
    here::here("R", "utils.R"),
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

cooper_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 5,
    t = 1,
    thresh = 1e-3
) {

  fit <- cooper::cooper(
    instance$data,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    stratify_by_status = TRUE,
    alpha = alpha,
    t = t,
    a = 0.5,
    thresh = thresh,
    include_mt_beta_history = FALSE
  )

  res = data.table::rbindlist(lapply(1:2, \(event) {
    rbind(
      data.table::data.table(
        variable = fit$predictors,
        estimate = coef(fit, event = event, use_initial_fit = TRUE),
        cause = event,
        method = "coxnet"
      ),
      data.table::data.table(
        variable = fit$predictors,
        estimate = coef(fit, event = event),
        cause = event,
        method = "cooper"
      )
    )
  }))

  res = merge(res, instance$effects, by = c("variable", "cause"))
  res[, error := (truth - estimate)]

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

  # Get results -------------------------------------------------------------
  res <-  reduceResultsDataTable()
  pars <- unwrap(getJobPars())

  nrow(res)/nrow(getJobTable())
  nrow(findErrors())/nrow(getJobTable())


  res_long <- data.table::rbindlist(lapply(res$job.id, \(id) {
    result <- res[(job.id == id), result][[1]]
    result
    result[, job.id := ..id]
  }))

  res_long <- merge(pars, res_long, by = "job.id")

  saveRDS(res, here::here("results", "1-results-poc.rds"))
  saveRDS(res_long, here::here("results", "1-results-long-poc.rds"))

} else {
  message("Don't know how to properly submit jobs on this platform!")
}