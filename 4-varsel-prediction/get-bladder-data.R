get_bladder_data <- function(split = 2/3, standardize = TRUE) {
  bladder_surv_geno <- readRDS(here::here("data/bladder_surv_geno.rds"))

  #prop.table(table(bladder_surv_geno$status))
  bladder_surv_geno <- data.table::as.data.table(bladder_surv_geno)

  bladder_surv_geno[, row_id := .I]
  dim(bladder_surv_geno)

  train <- bladder_surv_geno[,.SD[sample(.N, split * .N)], by = status]
  test <- bladder_surv_geno[setdiff(row_id, train$row_id), ]
  checkmate::assert_set_equal(union(train$row_id, test$row_id), bladder_surv_geno$row_id)
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

  status_distr_orig <- prop.table(table(bladder_surv_geno$status))
  status_distr_train <- prop.table(table(train$status))
  status_distr_test <- prop.table(table(train$status))

  checkmate::assert_true(all(round(status_distr_orig - status_distr_test, 2) == 0))
  checkmate::assert_true(all(round(status_distr_train - status_distr_test, 3) == 0))

  list(
    train = train,
    test = test,
    split = split,
    train_event_prop = status_distr_train,
    test_event_prop = status_distr_test
  )
}
