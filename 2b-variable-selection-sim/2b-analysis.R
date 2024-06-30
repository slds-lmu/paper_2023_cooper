source(here::here("R/utils.R"))
if (!dir.exists(here::here("results"))) dir.create(here::here("results"))

library(batchtools)
library(data.table)
reg_name <- "varsel-sim-pred"
reg_dir <- here::here("registries", reg_name)

loadRegistry(reg_dir, writeable = FALSE)


res = reduceResultsDataTable()
pars = unwrap(getJobPars())
res$result[[1]]$varsel
res$result[[1]]$scores

res_varsel <- rbindlist(lapply(res$job.id, \(id) {
  varsel_dt = as.data.table(res[job.id == id, result][[1]][["varsel"]])
  varsel_dt[, job.id := id]
  varsel_dt
}))

res_perf <- rbindlist(lapply(res$job.id, \(id) {
  perf_dt = as.data.table(res[job.id == id, result][[1]][["scores"]])
  perf_dt[, job.id := id]
  perf_dt
}))

res_varsel <- ljoin(res_varsel, pars)
res_perf <- ljoin(res_perf, pars)


aggr_perf = res_perf[, .(mean_score = mean(score)),
                     by = .(model, metric, cause, time_quant, lambda1, lambda2)]

