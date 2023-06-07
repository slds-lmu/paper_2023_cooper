#! /usr/bin/env Rscript
# monitor current status
#

# infer registry name from cli args
argv <- commandArgs(TRUE)

library(batchtools)

if (length(argv) == 0) {
  # If no argument, list all registries
  reg_name <- list.files(here::here("registries"))
} else {
  reg_name <- argv[[1]]
  checkmate::assert_directory_exists(here::here("registries", reg_name))
}

for (current_reg in reg_name) {

  reg_dir <- here::here("registries", current_reg)
  reg <- suppressMessages(loadRegistry(reg_dir, writeable = FALSE, make.default = FALSE))

  cli::cli_h1("Current Status: {current_reg}")
  print(getStatus(reg = reg))
  cat("\n")

  # Running -----------------------------------------------------------------
  tbl_running <- unwrap(getJobTable(findRunning(reg = reg), reg = reg))
  if (nrow(tbl_running) > 0) {
    cli::cli_h2("Running: {current_reg}")
    tbl_running[, c("job.id", "time.running", "problem", "algorithm")]
    print(tbl_running[, .(count = .N), by = algorithm])
  }

  # Done --------------------------------------------------------------------
  tbl_done <- unwrap(getJobTable(findDone(reg = reg), reg = reg))
  if (nrow(tbl_done) > 0) {
    cli::cli_h2("Done: {current_reg}")
    tbl_done <- tbl_done[, c("job.id", "time.running", "problem", "algorithm")]
    print(tbl_done[, .(count = .N), by = algorithm])
  }

  cat("\n")

  # Expired -----------
  tbl_expired <- unwrap(getJobTable(findExpired(reg = reg), reg = reg))
  if (nrow(tbl_expired) > 0) {
    cli::cli_h2("Expired: {current_reg}")
    tbl_expired <- tbl_expired[, c("job.id", "time.running", "problem", "algorithm")]
    print(tbl_expired[, .(count = .N), by = algorithm])
  }

  cat("\n")

  # Error'd -----------------------------------------------------------------
  tbl_errors <- unwrap(getJobTable(findErrors(reg = reg), reg = reg))
  if (nrow(tbl_errors) > 0) {
    cli::cli_h2("Errors: {current_reg}")
    tbl_errors <- tbl_errors[, c("job.id", "time.running", "problem", "algorithm")]
    print(tbl_errors[, .(count = .N), by = algorithm])
    print(getErrorMessages(findErrors(reg = reg), reg = reg))
  }

}
