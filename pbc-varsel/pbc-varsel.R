library(data.table)
library(fwelnet)

get_pbc_data <- function(job = NULL, data = NULL, na.rm = TRUE, num_noise = 0) {
  checkmate::assert_int(num_noise, lower = 0)

  ee = new.env(parent = emptyenv(), hash = FALSE)
  data(list = "Pbc3", package = "pec", envir = ee)
  xdf = ee[["Pbc3"]]

  if (na.rm) xdf <- na.omit(xdf)

  xdf$time <- xdf$days
  xdf$ptno <- NULL # patient number
  xdf$days <- NULL # recoded time
  xdf$unit <- NULL # Must be factor, but then too few events

  xdf <- data.table::as.data.table(xdf)

  if (num_noise > 0) {
    noise <- matrix(rnorm(nrow(xdf) * num_noise), ncol = num_noise)
    noise <- data.table::data.table(noise)
    names(noise) <- paste0("xnoise", seq_len(num_noise))
    xdf <- cbind(xdf, noise)
  }
  xdf
}

recode_cr_status <- function(xdf, target_event = 1) {
  checkmate::assert_subset(target_event, choices = 1:2)

  if (target_event == 1) {
    xdf$status[xdf$status == 2] <- 1
  } else if (target_event == 2) {
    xdf$status[xdf$status == 1] <- 0
    xdf$status[xdf$status == 2] <- 1
  }

  xdf
}


# pbc_data <- get_pbc_data(na.rm = TRUE, num_noise = 4)
# target_event <- 1
# pbc_varsel_score(
#   instance = get_pbc_data(na.rm = TRUE, num_noise = 4),
#   target_event = 1,
#   conf_level = 0.95,
#   mt_max_iter = 2, t = 100, thresh = 1e-7, alpha = 1
# )


pbc_varsel_score <- function(job = NULL, data = NULL, instance,
                             target_event = 1, conf_level = 0.95,
                             mt_max_iter = 3, t = 100, thresh = 1e-7, alpha = 1) {
  pbc_data <- instance

  # Fitting CSC should be part of data generating function I guess but otherwise inconvenient
  # "True" covariates depend on confidence level, so easier to re-fit and subset dymanically.
  mod_cph <- survival::coxph(
    survival::Surv(time, status) ~ .,
    # Fit only on original data, not noise vars
    data = recode_cr_status(pbc_data[, which(!startsWith(names(pbc_data), "xnoise")), with = FALSE],
                            target_event = target_event)
  )
  cph_coefs <- broom::tidy(mod_cph)
  coefs_true <- cph_coefs$term[cph_coefs$p.value < (1 - conf_level)]


  cofit <- fwelnet::cooper(
    pbc_data,
    mt_max_iter = mt_max_iter,
    stratify_by_status = TRUE,
    t = t,
    thresh = thresh
  )

  pbc_covars = setdiff(names(pbc_data), c("time", "status"))
  noise_vars = pbc_covars[startsWith(pbc_covars, "xnoise")]


  data.table::rbindlist(list(
    score_model(cofit, target_event = target_event, model =  "cooper", coefs_true, noise_vars, pbc_covars),
    score_model(cofit, target_event = target_event, model =  "coxnet", coefs_true, noise_vars, pbc_covars)
  ))
}

selected_vars <- function(copper_fit, event = 1, use_initial_fit = FALSE) {
  res <- coef(copper_fit, event = event, use_initial_fit = use_initial_fit)
  names(res[res != 0])
}

not_selected_vars <- function(copper_fit, event = 1, use_initial_fit = FALSE) {
  res <- coef(copper_fit, event = event, use_initial_fit = use_initial_fit)
  names(res[res == 0])
}

score_model <- function(cooper_fit, target_event = 1, model = "cooper", coefs_true, noise_vars, pbc_covars) {
  checkmate::assert_subset(model, choices = c("cooper", "coxnet"))


  cooper_selected <- selected_vars(cooper_fit, event = target_event, use_initial_fit = model == "coxnet")
  cooper_not_selected <- not_selected_vars(cooper_fit, event = target_event, use_initial_fit = model == "coxnet")

  which_positive <- pbc_covars %in% coefs_true
  # total_neg based only on noise vars, treating non-signif original vars as neither true nor false
  which_negative <- pbc_covars %in% noise_vars

  total_positive <- sum(which_positive)
  total_negative <- sum(which_negative)

  result <- data.table::data.table(
    total = total_positive + total_negative,
    total_pos = sum(which_positive),
    total_neg = sum(which_negative),
    tp = sum(cooper_selected %in% pbc_covars[which_positive]),
    fp = sum(cooper_selected %in% pbc_covars[which_negative]),
    tn = sum(cooper_not_selected %in% pbc_covars[which_negative]),
    fn = sum(cooper_not_selected %in% pbc_covars[which_positive]),
    selected = list(cooper_selected),
    true = list(coefs_true)
  )

  result[, ppv := tp/(tp + fp)]
  result[, fpr := fp/total_neg]
  result$model <- model
  result
}
