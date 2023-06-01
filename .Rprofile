source("renv/activate.R")

Sys.setenv(OMP_NUM_THREADS="1")
Sys.setenv(OPENBLAS_NUM_THREADS="1")
Sys.setenv(MKL_NUM_THREADS="1")

options(
  datatable.print.class = TRUE,
  datatable.print.keys = TRUE,
  batchtools.progress = FALSE # for speedup
)

# Load helper functions in R/
#purrr::walk(list.files(here::here("R"), pattern = "*.R", full.names = TRUE), source, echo = FALSE)
