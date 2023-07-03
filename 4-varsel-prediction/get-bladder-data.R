get_bladder_data <- function(job = NULL, data = NULL, split = 2/3, standardize = TRUE, type = c("clinical", "geno", "both")) {

  bladder_file <- switch(
    as.character(type),
    "clinical" = here::here("data/bladder-binder-clinical.rds"),
    "geno" = here::here("data/bladder-binder-geno.rds"),
    "both" = here::here("data/bladder-binder-clinical_geno.rds")

    # old versions
    # "clinical" = here::here("data/bladder_surv_clinical.rds"),
    # "geno" = here::here("data/bladder_surv_geno.rds"),
    # "both" = here::here("data/bladder_surv_clin_geno.rds")
  )

  if (!file.exists(bladder_file)) {
    message("Recreating bladder dataset from bladder-data-prep.R")
    #source(here::here("4-varsel-prediction/bladder-data-prep.R"))
    source(here::here("4-varsel-prediction/preprocess-orig-binder.R"))
  }
  bladder_surv <- readRDS(bladder_file)

  #prop.table(table(bladder_surv_geno$status))
  bladder_surv <- data.table::as.data.table(bladder_surv)

  bladder_surv[, row_id := .I]
  dim(bladder_surv)

  train <- bladder_surv[ , .SD[sample(.N, size = round(split * .N), replace = FALSE)], by = status]
  test <- bladder_surv[setdiff(row_id, train$row_id), ]
  checkmate::assert_set_equal(union(train$row_id, test$row_id), bladder_surv$row_id)
  #checkmate::assert_true(length(setdiff(train$time, test$time)) == length(train$time)) # non-unique time points

  # Remove row_ids so we don't accidentally standardize or train on them,
  # but store them to triple check sampling works as expected
  train_row_ids <- train$row_id
  test_row_ids <- test$row_id
  train[, row_id := NULL]
  test[, row_id := NULL]


  if (standardize) {
    # get column-wise means and sd of genetic data only (1, 2 are time, status), exclude possible factor vars
    x_numerics <- sapply(train, \(x) is.numeric(x) & !setequal(unique(x), 0:1))
    targets <- names(train) %in% c("time", "status")

    x_to_standardize <- x_numerics & !targets
    vars_remaining <- !(x_to_standardize)

    gen_means <- apply(train[, x_to_standardize, with = FALSE], 2, mean)
    gen_sds <- apply(train[, x_to_standardize, with = FALSE], 2, sd)
    # apply to only genetic data, cbind with time, status, etc. There's probably a more elegant solution, sorry.
    train <- cbind(train[, vars_remaining, with = FALSE], scale(train[, x_to_standardize, with = FALSE], center = gen_means, scale = gen_sds))
    test <- cbind(test[, vars_remaining, with = FALSE], scale(test[, x_to_standardize, with = FALSE], center = gen_means, scale = gen_sds))
  }

  checkmate::assert_count(nrow(train), positive = TRUE)
  checkmate::assert_count(nrow(test), positive = TRUE)
  checkmate::assert_true(ncol(train) == ncol(test))
  checkmate::assert_true(nrow(train) + nrow(test) == nrow(bladder_surv))
  # checkmate::assert_set_equal(setdiff(train$time, test$time), train$time)

  status_distr_orig <- prop.table(table(bladder_surv$status))
  status_distr_train <- prop.table(table(train$status))
  status_distr_test <- prop.table(table(test$status))

  checkmate::assert_true(all(round(status_distr_orig - status_distr_test, 1) == 0))
  checkmate::assert_true(all(round(status_distr_train - status_distr_test, 1) == 0))

  list(
    train = train,
    test = test,
    split = split,
    train_row_ids = train_row_ids,
    test_row_ids = test_row_ids,
    train_event_prop = status_distr_train,
    test_event_prop = status_distr_test
  )
}

if (FALSE) {
  get_bladder_data(standardize = FALSE, type = "clinical")
  get_bladder_data(standardize = TRUE, type = "clinical")

  get_bladder_data(standardize = FALSE, type = "geno")
  get_bladder_data(standardize = TRUE, type = "geno")

  get_bladder_data(standardize = FALSE, type = "both")
  get_bladder_data(standardize = TRUE, type = "both")
}
