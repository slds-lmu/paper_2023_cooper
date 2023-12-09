library(data.table)
library(fwelnet)
library(survival)

if (FALSE) {

  instance <- get_pbc_data(na.rm = TRUE, num_noise = 4, conf_level = 0.95)

  cooper_varsel_pbc(instance = instance)
  rfsrc_pbc_varsel(instance = instance)
  coxboost_pbc_varsel(instance = instance)

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

fit_csc_wrap <- function(xdf, conf_level = 0.95, target_event = 1) {
  # Fitting CSC should be part of data generating function I guess but otherwise inconvenient
  # "True" covariates depend on confidence level, so easier to re-fit and subset dymanically.
  mod_cph <- survival::coxph(
    survival::Surv(time, status) ~ .,
    # Fit only on original data, not noise vars
    data = recode_cr_status(xdf, target_event = target_event)
  )
  cph_coefs <- broom::tidy(mod_cph)
  coefs_true <- cph_coefs$term[cph_coefs$p.value < (1 - conf_level)]

  list(coefs_true = coefs_true, conf_level = conf_level, target_event = target_event)
}


get_pbc_data <- function(job = NULL, data = NULL, na.rm = TRUE, num_noise = 0, conf_level = 0.95) {
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

  result = list(
    effects = list(
      fit_csc_wrap(xdf, conf_level = conf_level, target_event = 1),
      fit_csc_wrap(xdf, conf_level = conf_level, target_event = 2)
  ))

  if (num_noise > 0) {
    noise <- matrix(rnorm(nrow(xdf) * num_noise), ncol = num_noise)
    noise <- data.table::data.table(noise)
    names(noise) <- paste0("xnoise", seq_len(num_noise))
    xdf <- cbind(xdf, noise)

    result$noise_vars <- names(noise)
  }

  result$data <- xdf
  result$covars <- setdiff(names(xdf), c("time", "status"))

  result

}


score_model <- function(model_fit, target_event = 1, model = "cooper", coefs_true, noise_vars, pbc_covars) {
  checkmate::assert_subset(model, choices = c("cooper", "coxnet", "coxboost", "rsfrc"))

  #browser()
  # I need you to know that I know how S3 classes work ok
  model_coefs <- switch(class(model_fit)[[1]],
                "cooper" = coef(model_fit, event = target_event, use_initial_fit = model == "coxnet"),
                "CoxBoost" =  coef(model_fit)[[target_event]],
                "rfsrc" = rsf_vita_selected(model_fit, target_event = target_event)
  )
  # Selected and not selected variable names
  model_selected <- names(model_coefs[model_coefs != 0])
  model_not_selected <-  names(model_coefs[model_coefs == 0])

  which_positive <- pbc_covars %in% coefs_true
  # total_neg based only on noise vars, treating non-signif original vars as neither true nor false
  which_negative <- pbc_covars %in% noise_vars

  total_positive <- sum(which_positive)
  total_negative <- sum(which_negative)

  result <- data.table::data.table(
    total = total_positive + total_negative,
    total_pos = sum(which_positive),
    total_neg = sum(which_negative),
    tp = sum(model_selected %in% pbc_covars[which_positive]),
    fp = sum(model_selected %in% pbc_covars[which_negative]),
    tn = sum(model_not_selected %in% pbc_covars[which_negative]),
    fn = sum(model_not_selected %in% pbc_covars[which_positive]),
    selected = list(model_selected),
    true = list(coefs_true)
  )

  result[, ppv := tp/(tp + fp)]
  result[, fpr := fp/total_neg]
  result$model <- model
  result$target_event <- target_event
  result
}


# Algorithms ---------------------------------------------------------------------------------

cooper_varsel_pbc <- function(job = NULL, data = NULL, instance,
                              conf_level = 0.95,
                              mt_max_iter = 3, t = 100, thresh = 1e-7, alpha = 1) {


  cooperfit <- fwelnet::cooper(
    instance$data,
    mt_max_iter = mt_max_iter,
    stratify_by_status = TRUE,
    t = t,
    thresh = thresh
  )

  data.table::rbindlist(list(
    score_model(cooperfit, target_event = 1, model =  "cooper", instance$effects[[1]]$coefs_true, instance$noise_vars, instance$covars),
    score_model(cooperfit, target_event = 2, model =  "cooper", instance$effects[[2]]$coefs_true, instance$noise_vars, instance$covars),
    score_model(cooperfit, target_event = 1, model =  "coxnet", instance$effects[[1]]$coefs_true, instance$noise_vars, instance$covars),
    score_model(cooperfit, target_event = 2, model =  "coxnet", instance$effects[[2]]$coefs_true, instance$noise_vars, instance$covars)
  ))
}

rfsrc_pbc_varsel <- function(data, job, instance,
                                    #mtry = 2000,
                                    nodesize = 15, splitrule = "logrank",
                                    importance = "random", cutoff_method = "vita"
) {
  # batchtools passes params as factors, need to convert
  importance <- as.character(importance)
  splitrule <- as.character(splitrule)
  cutoff_method <- as.character(cutoff_method)

  # prbly won't use permute, since random should be equivalent
  checkmate::assert_subset(importance, choices = c("random", "permute", "anti"))
  checkmate::assert_subset(cutoff_method, choices = "vita")

  rf_c1 <- randomForestSRC::rfsrc(
    Surv(time, status) ~ ., data = instance$data,
    splitrule = splitrule,
    cause = c(1, 0),
    importance = importance,
    # mtry = mtry,
    nodesize = nodesize
  )

  rf_c2 <- randomForestSRC::rfsrc(
    Surv(time, status) ~ ., data = instance$data,
    splitrule = splitrule,
    cause = c(0, 1),
    importance = importance,
    # mtry = mtry,
    nodesize = nodesize
  )

  data.table::rbindlist(list(
    score_model(rf_c1, target_event = 1, model =  "rsfrc", instance$effects[[1]]$coefs_true, instance$noise_vars, instance$covars),
    score_model(rf_c2, target_event = 2, model =  "rsfrc", instance$effects[[2]]$coefs_true, instance$noise_vars, instance$covars)
  ))
}

rsf_vita_selected <- function(rf, target_event = 1) {
  vimps <- data.table::data.table(
    variable = rownames(rf[["importance"]]),
    vi_c = rf[["importance"]][, target_event]
  )

  variable = rownames(rf[["importance"]])
  vi_c = rf[["importance"]][, target_event]

  # Apply vita/janitza shortcut
  # Take absolute value of minimal importance value as cut-off for 0-classification
  # Effectively equivalent with vita method in Degenhardt et al. (2019)
  vi_c[vi_c <= abs(min(vi_c))] <- 0

  names(vi_c) <- variable
  vi_c
}


coxboost_pbc_varsel <- function(data, job, instance,
                                cmprsk = "csh", stepno = 100, penalty = 2000
) {

  cbfit <- CoxBoost::CoxBoost(
    time = instance$data$time,
    status = instance$data$status,
    x = as.matrix(instance$data[, instance$covars, with = FALSE]),
    cmprsk = as.character(cmprsk),
    stepno = stepno,
    # Default is 9 * sum(status[subset] == 1), should be approx 9 * 250
    # since in our setting sum(instance$train$status != 0) is approx 233
    penalty = penalty
  )

  data.table::rbindlist(list(
    score_model(cbfit, target_event = 1, model =  "coxboost", instance$effects[[1]]$coefs_true, instance$noise_vars, instance$covars),
    score_model(cbfit, target_event = 2, model =  "coxboost", instance$effects[[2]]$coefs_true, instance$noise_vars, instance$covars)
  ))
}

