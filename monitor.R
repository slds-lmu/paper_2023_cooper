#! /usr/bin/env Rscript
# monitor current status
library(batchtools)
reg_dir <- here::here("registries", "fwel_sim_varsel")
loadRegistry(reg_dir, writeable = FALSE)


# Running -----------------------------------------------------------------
cli::cli_h1("Running")
getStatus()

cat("\n")

tbl_running <- unwrap(getJobTable(findRunning()))
if (nrow(tbl_running) > 0) {
  tbl_running[, c("job.id", "time.running", "algorithm")]
  tbl_running[, .(count = .N), by = algorithm]
}

# Done --------------------------------------------------------------------
cli::cli_h1("Done")
tbl_done <- unwrap(getJobTable(findDone()))
if (nrow(tbl_done) > 0) {
  tbl_done <- tbl_done[, c("job.id", "time.running", "algorithm")]
  tbl_done[, .(count = .N), by = algorithm]
}
cat("\n")

# Expired -----------
cli::cli_h1("Expired")
tbl_expired <- unwrap(getJobTable(findExpired()))
if (nrow(tbl_expired) > 0) {
  tbl_expired <- tbl_expired[, c("job.id", "time.running", "algorithm")]
  tbl_expired[, .(count = .N), by = algorithm]
}
cat("\n")

# Error'd -----------------------------------------------------------------
cli::cli_h1("Errors")
getErrorMessages(findErrors())
