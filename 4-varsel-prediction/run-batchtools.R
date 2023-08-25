source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))
source(here::here("2-variable-selection-sim/20-varsel-simulation.R"))

library(batchtools)
library(randomForestSRC)
library(CoxBoost)
library(riskRegression)
library(rlang)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 100
)

set.seed(config$global_seed)
# Set to FALSE to remove + recreate registry
continue_bt <- FALSE

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "fwel_sim_varsel_predict_2309"
reg_dir <- here::here("registries", reg_name)

if (continue_bt) {
  message("Continuing with existing registry\n")
  loadRegistry(reg_dir, writeable = TRUE)
} else {
  message("Creating new registry!\n")
  unlink(reg_dir, recursive = TRUE)
  makeExperimentRegistry(file.dir = reg_dir, packages = c("randomForestSRC", "CoxBoost", "rlang", "data.table", "riskRegression"),
                         seed = config$global.seed,
                         source = c(here::here("4-varsel-prediction/algorithms.R"),
                                    here::here("4-varsel-prediction/get-bladder-data.R"),
                                    here::here("2-variable-selection-sim/20-varsel-simulation.R"))
  )
}

# Problems -----------------------------------------------------------
addProblem(name = "bladder_geno", fun = get_bladder_data, seed = config$sim_seed, cache = config$sim_cache)
#addProblem(name = "binder_bender", fun = sim_surv_binder_resample, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "fwel_mt", fun = fwel_mt_varselect_pred)
addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_pred)
addAlgorithm(name = "coxboost", fun = coxboost_varselect_pred)


# Experiments -----------------------------------------------------------
prob_design <- list(
  bladder_geno = expand.grid(
    split = 2/3, standardize = TRUE,
    type = c("clinical", "both") # c("clinical", "geno", "both")
  )
  # ,
  # binder_bender = expand.grid(
  #   n_train = 400, n_test = 200, p = 5000,
  #   ce = c(0.5),
  #   lambda1 = 0.1, lambda2 = 0.1, lambda_c = 0.1
  # )
)

algo_design <- list(
  fwel_mt = expand.grid(
    mt_max_iter = 3,
    alpha = c(1, 0.75),
    t = 100,
    thresh = 1e-7,
    stratify_by_status = TRUE,
    nfolds = 5
  ),
  rfsrc = expand.grid(
    importance = "random",
    cutoff_method = "vita",
    mtry = c(1000),
    nodesize = 30,
    splitrule = "logrank"
  ),
  coxboost = expand.grid(
    cmprsk = "csh",
    stepno = 117, # c(100, 300), # 117 is what binder et al report for bladder data
    penalty = 3000 # c(1000, 3000)
  )
)


addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
jobtbl <- unwrap(getJobPars(), c("algo.pars", "prob.pars"))

# Test jobs -----------------------------------------------------------

# jobtbl[algorithm == "fwel_mt"]

# if (interactive()) {
#   testJob(228)
#   testJob(1885)
# }


# Submit -----------------------------------------------------------

testJob(3)
testJob(439)
testJob(684)

jobtbl[algorithm == "fwel_mt", .SD[sample(.N, 5)], by = c("type")] |>
  findNotSubmitted() |>
  submitJobs()


jobtbl[, .SD[sample(.N, 10)], by = c("problem", "algorithm")] |>
  findNotSubmitted() |>
  submitJobs()

jobtbl[findDone(), .(count = .N), by = algorithm]

submitJobs(findNotSubmitted(jobtbl[algorithm == "fwel_mt" & problem == "binder_bender"]))
submitJobs(findNotSubmitted())

# submitJobs(jobtbl[algorithm == "coxboost"])
# submitJobs(jobtbl[algorithm == "rfsrc"])



# Monitor jobs ------------------------------------------------------------
if (interactive()) {
  getStatus()

  ijoin(
    unwrap(getJobPars(findErrors()), c("algo.pars", "prob.pars")),
    getErrorMessages(findErrors())
  )
}

# Get results -------------------------------------------------------------
# res <-  ijoin(reduceResultsDataTable(), unwrap(getJobPars(), c("prob.pars", "algo.pars")))


