source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("2-variable-selection-sim/20-varsel-simulation.R"))
source(here::here("5-prediction/5-prediction-wrapper.R"))


library(batchtools)
library(survival)
library(riskRegression)

# Settings ----------------------------------------------------------------
config <- list(
  global_seed = 563,
  sim_seed = 569,
  sim_cache = FALSE,
  repls = 100
)

# set.seed(config$global_seed)

# Registry ----------------------------------------------------------------
if (!file.exists(here::here("registries"))) dir.create(here::here("registries"))
reg_name <- "5-prediction"
reg_dir <- here::here("registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeExperimentRegistry(file.dir = reg_dir, packages = c("randomForestSRC", "CoxBoost", "data.table", "riskRegression"),
                       seed = config$global.seed,
                       source = c(here::here("4-varsel-prediction/algorithms.R"),
                                  here::here("4-varsel-prediction/get-bladder-data.R"),
                                  here::here("2-variable-selection-sim/20-varsel-simulation.R"))
)
# Problems -----------------------------------------------------------
addProblem(name = "sim_surv_binder", fun = sim_surv_binder, seed = config$sim_seed, cache = config$sim_cache)
addProblem(name = "bladder_geno", fun = get_bladder_data, seed = config$sim_seed, cache = config$sim_cache)

# Algorithms -----------------------------------------------------------
addAlgorithm(name = "fwel_mt", fun = fwel_mt_prediction_wrapper)
addAlgorithm(name = "coxboost", fun = coxboost_varselect_pred)
# addAlgorithm(name = "rfsrc", fun = rfsrc_varselect_pred)


# Experiments -----------------------------------------------------------
prob_design <- list(
  sim_surv_binder = expand.grid(
    n_train = 400,
    n_test = 200,
    p = 5000, ce = c(0.25, 0.5),
    # cause 1 either slightly higher prevalence or noticeably lower
    lambda1 = c(0.1, 0.01),
    lambda_c = 0.1)
  ,
  bladder_geno = expand.grid(
    split = 2/3, standardize = TRUE,
    type = c("clinical", "both") # c("clinical", "geno", "both")
  )
)

algo_design <- list(
  fwel_mt = expand.grid(
    mt_max_iter = 3,
    alpha = c(0.75, 1),
    t = c(100),
    thresh = c(1e-3, 1e-5)
  ),
  # rfsrc = expand.grid(
  #   importance = "random",
  #   cutoff_method = "vita",
  #   mtry = c(1000),
  #   nodesize = 30,
  #   splitrule = "logrank"
  # ),
  coxboost = expand.grid(
    cmprsk = "csh",
    stepno = 117, # c(100, 300), # 117 is what binder et al report for bladder data
    penalty = 3000 # c(1000, 3000)
  )
)

addExperiments(prob_design, algo_design, repls = config$repls)
summarizeExperiments()
unwrap(getJobPars(), c("algo.pars", "prob.pars"))#[algorithm == "fwel_mt" & problem == "bladder_geno",]

# Test jobs -----------------------------------------------------------
if (FALSE) {
  testJob(4)
  testJob(2997)
}

# Submit -----------------------------------------------------------

ids <- findNotStarted()
submitJobs(ids = ids)
# waitForJobs()


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

