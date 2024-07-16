source("renv/activate.R")

# Ensure no accidental nested parallelization
Sys.setenv(OMP_NUM_THREADS="1")
Sys.setenv(OPENBLAS_NUM_THREADS="1")
Sys.setenv(MKL_NUM_THREADS="1")

options(
  datatable.print.class = TRUE,
  datatable.print.keys = TRUE
)

# Load helper functions in R/
invisible(lapply(list.files("R", pattern = "*.R", full.names = TRUE), source, echo = FALSE))