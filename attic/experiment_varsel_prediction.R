# dry run varsel prediction with simulated data
#source(here::here("4-varsel-prediction/get-bladder-data.R"))
source(here::here("4-varsel-prediction/algorithms.R"))
#source(here::here("2-variable-selection-sim/20-varsel-algorithms.R"))

library(fwelnet)
library(survival)
library(ggplot2)
library(data.table)
library(riskRegression)

# Standardize predictors
standardize <- TRUE
# holdout split ratio
split <- 2/3

# Get some dataset (sim as placeholder)
xdat <- prodlim::SimCompRisk(500)
xdat <- xdat[c("time", "event", "X1", "X2")]
names(xdat) <- c("time", "status", "X1", "X2")
xdat <- as.data.table(xdat)

# Filter missings if necessary, cry if it breaks everything
xdat <- na.omit(xdat)
checkmate::assert_data_table(xdat, any.missing = FALSE, min.rows = 100, min.cols = 3)

# Need to dummy-encode (without intercept, hence -1, but this fucks up dichotomous factors, so 0/1 encode them beforehand!)
xdat <- as.data.table(model.matrix(~ -1 + ., data = xdat))
# strip spaces from varnames just in case
names(xdat) <- gsub(" ", "", names(xdat))


check_status <- function(xs) {
  tab <- as.integer(table(xs))
  mtab <- rbind(events = tab, prop = 100 * round(tab/sum(tab), 2))
  colnames(mtab) <- unique(xs)
  mtab <- cbind(mtab, sum = rowSums(mtab))
  print(mtab)
}

check_status(xdat$status)

# sanity KM? unify both events just to check
xdat_c1 <- data.table::copy(xdat)
xdat_c1[, status := ifelse(status == 2, 1, status)]
survfit(Surv(time, status) ~ 1, data = xdat_c1) |>
  plot()

# train/test split ----
xdat[, row_id := .I]
dim(xdat)

# stratified sampling by status
train <- xdat[ , .SD[sample(.N, size = round(split * .N), replace = FALSE)], by = status]
test <- xdat[setdiff(row_id, train$row_id), ]
checkmate::assert_set_equal(union(train$row_id, test$row_id), xdat$row_id)
checkmate::assert_true(length(setdiff(train$time, test$time)) == length(train$time))
check_status(train$status)
check_status(test$status)

# Remove row_ids so we don't accidentally standardize or train on them,
# but store them to triple check sampling works as expected
train_row_ids <- train$row_id
test_row_ids <- test$row_id
train[, row_id := NULL]
test[, row_id := NULL]

# Standardize ----
if (standardize) {
  # get column-wise means and sd of genetic data only (1, 2 are time, status)
  # exclude possible factor vars and 0/1 vars from standardization
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
checkmate::assert_true(nrow(train) + nrow(test) == nrow(xdat))
checkmate::assert_set_equal(setdiff(train$time, test$time), train$time)

# create instance object for batchtools pipeline
test_instance <- list(
  train = train, test = test
)

# Fit fwelnet (incl glmnet)  ----
# coxph prediction on selected vars is part of the pipeline alread

fwfit <- fwel_mt_varselect_pred(
  data = NULL, job = NULL,
  instance = test_instance,
  alpha = 1,
  mt_max_iter = 2,
  a = 0.5,
  t = 100,
  thresh = 1e-7
)

# Estimated coefs from cause-specific models based on either algorithm
fwfit$coefs$fwelnet$cause1
fwfit$coefs$glmnet$cause1

fwfit$coefs$fwelnet$cause2
fwfit$coefs$glmnet$cause2
# Don't ask me why I'm storing these weirdly pls it's all scuffed

# Compare cause-1 selected vars
selected_c1 <- lapply(fwfit$coefs, \(x) names(x[[1]]))
setdiff(selected_c1$glmnet, selected_c1$fwelnet)


# Plot tAUC / Brier for fwelnet/glmnet (remove IPA for now)
fwfit$scores |>
  dplyr::filter(metric %in% c("AUC", "Brier")) |>
  ggplot(aes(x = times, y = score, color = model, fill = model)) +
  facet_wrap(vars(metric, cause), scales = "free", labeller = labeller(cause = label_both)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = \(x) x * 100) +
  labs(
    title = "More Scores I guess",
    subtitle = "Data is hard",
    x = "Event Times", y = "Score (%)",
    caption = "Grmrmrlmlm"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title.position = "plot")
