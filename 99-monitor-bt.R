library(batchtools)

reg_dir <- here::here("registries", "fwel_sim_variableselect")
loadRegistry(reg_dir, writeable = FALSE)

getStatus()
