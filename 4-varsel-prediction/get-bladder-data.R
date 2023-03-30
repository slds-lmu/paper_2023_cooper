get_bladder_data <- function(job = NULL, data = NULL, split = 2/3, standardize = TRUE) {
# browser()
  bladder_file <- here::here("data/bladder_surv_geno.rds")
  if (!file.exists(bladder_file)) {
    message("Recreating bladder dataset/bladder-data-prep.R")
    source(here::here("4-varsel-prediction/bladder-data-prep.R"))
  }
  bladder_surv_geno <- readRDS(bladder_file)

  #prop.table(table(bladder_surv_geno$status))
  bladder_surv_geno <- data.table::as.data.table(bladder_surv_geno)

  bladder_surv_geno[, row_id := .I]
  dim(bladder_surv_geno)

  train <- bladder_surv_geno[ , .SD[sample(.N, size = round(split * .N), replace = FALSE)], by = status]
  test <- bladder_surv_geno[setdiff(row_id, train$row_id), ]
  checkmate::assert_set_equal(union(train$row_id, test$row_id), bladder_surv_geno$row_id)
  checkmate::assert_true(length(setdiff(train$time, test$time)) == length(train$time))

  # Remove row_ids so we don't accidentally standardize or train on them,
  # but store them to triple check sampling works as expected
  train_row_ids <- train$row_id
  test_row_ids <- test$row_id
  train[, row_id := NULL]
  test[, row_id := NULL]


  if (standardize) {
    # get column-wise means and sd of genetic data only (1, 2 are time, status)
    gen_means <- apply(train[, -c(1,2)], 2, mean)
    gen_sds <- apply(train[, -c(1,2)], 2, sd)
    # apply to only genetic data, cbind with time, status. There's probably a more elegant solution, sorry.
    train <- cbind(train[, 1:2], scale(train[, -c(1,2)], center = gen_means, scale = gen_sds))
    test <- cbind(test[, 1:2], scale(test[, -c(1,2)], center = gen_means, scale = gen_sds))
  }

  checkmate::assert_count(nrow(train), positive = TRUE)
  checkmate::assert_count(nrow(test), positive = TRUE)
  checkmate::assert_true(ncol(train) == ncol(test))
  checkmate::assert_true(nrow(train) + nrow(test) == nrow(bladder_surv_geno))
  checkmate::assert_set_equal(setdiff(train$time, test$time), train$time)

  status_distr_orig <- prop.table(table(bladder_surv_geno$status))
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
