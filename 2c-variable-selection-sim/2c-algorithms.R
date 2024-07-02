# Variable selection ----------------------------------------------------------------------------------------------

# cooper wrapper for variable selection sim
cooper_varsel_wrapper <- function(
    data, job, instance,
    alpha = 1,
    z_method = "original",
    mt_max_iter = 2,
    t = 100,
    a = 0.5,
    thresh = 1e-7
) {

  checkmate::assert_number(alpha, lower = 0, upper = 1)
  checkmate::assert_int(mt_max_iter, lower = 1)
  checkmate::assert_number(t, lower = 0)
  checkmate::assert_number(a, lower = 0)
  checkmate::assert_number(thresh, lower = 0)

  # instance <- sim_surv_binder(n_train = 400, p = 5000)
  # instance <- sim_surv_binder(n_train = 200, p = 500, n_test = 200)
  # alpha = 1
  # z_method = "original"
  # mt_max_iter = 2
  # t = 100
  # a = 0.5
  # thresh = 1e-3

  fit <- cooper::cooper(
    instance$train,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh,
    stratify_by_status = TRUE,
    include_mt_beta_history = TRUE
  )

  p <- ncol(instance$train)
  truth <- instance$covar_true_effect
  total <- instance$covar_blocks

  # glmnet_beta1 <- fit$beta1[, 1]
  # glmnet_beta2 <- fit$beta2[, 1]
  glmnet_beta1 <- coef(fit, event = 1, use_initial_fit = TRUE)
  glmnet_beta2 <- coef(fit, event = 2, use_initial_fit = TRUE)

  # Same for cooper estimates
  cooper_beta1 <- coef(fit, event = 1)
  cooper_beta2 <- coef(fit, event = 2)

  # Get confusion matrix per block of predictors, by model, by cause.
  # I am painfully aware that this is not "nice", but time is finite and patience is, too.

  res_varsel <- data.table::rbindlist(c(
    lapply(names(total), function(x) {
      rbind(
      get_confusion(glmnet_beta1, truth, total, x,  model = "glmnet", cause = 1L),
      get_confusion(glmnet_beta2, truth, total, x,  model = "glmnet", cause = 2L),
      get_confusion(cooper_beta1, truth, total, x,  model = "cooper", cause = 1L),
      get_confusion(cooper_beta2, truth, total, x,  model = "cooper", cause = 2L)
      )
    })
  ))


 # Prediction
  scores_cmb <- data.table::rbindlist(list(
    fit_csc_coxph(instance, model = "cooper", coefs = nonzeros(cooper_beta1), cause = 1),
    fit_csc_coxph(instance, model = "cooper", coefs = nonzeros(cooper_beta2), cause = 2),
    fit_csc_coxph(instance, model = "glmnet", coefs = nonzeros(glmnet_beta1), cause = 1),
    fit_csc_coxph(instance, model = "glmnet", coefs = nonzeros(glmnet_beta2), cause = 2)
  ))

  list(
    varsel = res_varsel,
    scores = scores_cmb
  )

}


rfsrc_varselect_wrapper <- function(data, job, instance,
                                    splitrule = "logrank",
                                    importance = "random", cutoff_method = "vita"
                                    ) {
  # batchtools passes params as factors, need to convert
  importance <- as.character(importance)
  splitrule <- as.character(splitrule)
  cutoff_method <- as.character(cutoff_method)

  # prbly won't use permute, since random should be equivalent
  checkmate::assert_subset(importance, choices = c("random", "permute", "anti"))
  # not sure about other methods yet
  checkmate::assert_subset(cutoff_method, choices = "vita")

  rf_c1 <- rfsrc_tuned(
    xdat = instance$train,
    splitrule = splitrule,
    importance = importance,
    cause = 1
  )

  rf_c2 <- rfsrc_tuned(
    xdat = instance$train,
    splitrule = splitrule,
    importance = importance,
    cause = 2
  )

  rf_models = list("1" = rf_c1, "2" = rf_c2)

  vimps <- data.table::data.table(
    variable = rownames(rf_c1[["importance"]]),
    vi_c1 = rf_c1[["importance"]][, 1],
    vi_c2 = rf_c2[["importance"]][, 2]
  )

  if (cutoff_method == "vita") {
    # Apply vita/janitza shortcut
    # Take absolute value of minimal importance value as cut-off for 0-classification
    # Effectively equivalent with vita method in Degenhardt et al. (2019)
    # Selecting
    vimps[, vita_c1 := ifelse(vi_c1 <= abs(min(vi_c1)), 0, vi_c1)]
    vimps[, vita_c2 := ifelse(vi_c2 <= abs(min(vi_c2)), 0, vi_c2)]
  }

  truth <- instance$covar_true_effect
  total <- instance$covar_blocks

  reslist_varsel = lapply(names(total), function(x) {
    rbind(
      get_confusion(vimps$vita_c1, truth, total, x,  model = "rfsrc", cause = 1L),
      get_confusion(vimps$vita_c2, truth, total, x,  model = "rfsrc", cause = 2L)
    )
  })

  res_varsel <- data.table::rbindlist(reslist_varsel)

  # Prediction
  scores_cmb <- data.table::rbindlist(list(
    fit_csc_coxph(instance, model = "rfsrc", coefs = selected(rf_c1)[["1"]], cause = 1),
    fit_csc_coxph(instance, model = "rfsrc", coefs = selected(rf_c2)[["2"]], cause = 2)
  ))

  list(
    varsel = res_varsel,
    scores = scores_cmb
  )
}

coxboost_varselect_wrapper <- function(data, job, instance, cmprsk = "csh") {
  cmprsk <- checkmate::assert_subset(as.character(cmprsk), choices = c("csh", "sh", "ccsh"))
  cbfit <- coxboost_tuned(instance$train, cmprsk = cmprsk)

  # Extracts coefs at final boosting step, names list per cause (names 1, 2)
  cb_coefs <- coef(cbfit)
  truth <- instance$covar_true_effect
  total <- instance$covar_blocks

  if (cmprsk == "sh") {
    # No cause-specific coefficients with subdistribution approach, so duplicate results (kind of)
    reslist_varsel <- lapply(names(total), function(x) {
      rbind(
        get_confusion(cb_coefs, truth, total, x,  model = "coxboost", cause = 1L),
        get_confusion(cb_coefs, truth, total, x,  model = "coxboost", cause = 2L)
      )
    })

  } else {
    reslist_varsel <- lapply(names(total), function(x) {
      rbind(
        get_confusion(cb_coefs[["1"]], truth, total, x,  model = "coxboost", cause = 1L),
        get_confusion(cb_coefs[["2"]], truth, total, x,  model = "coxboost", cause = 2L)
      )
    })
  }

  res_varsel <- data.table::rbindlist(reslist_varsel)

  # Prediction
  scores_cmb <- data.table::rbindlist(list(
    fit_csc_coxph(instance, model = "coxboost", coefs = nonzeros(cb_coefs[["1"]]), cause = 1),
    fit_csc_coxph(instance, model = "coxboost", coefs = nonzeros(cb_coefs[["2"]]), cause = 2)
  ))

  list(
    varsel = res_varsel,
    scores = scores_cmb
  )
}


#' Get confusion matrix stats from vector of coefficients and true effect info
#'
#' Only to be used with `sim_surv_binder` as it returns necessary covar block structure etc.
#'
#' @param beta Vector of estimated coefficients, numeric
#' @param truth List of per-block true effects, `instance$covar_true_effect`
#' @param total List of per-block indices, `instance$covar_blocks`
#' @param block Block ID, character
#' @param model Model identifier, e.g. `"cooper"`
#' @param cause Integer cause 0 or 1
#'
#' @return A 1-row data.frame

get_confusion <- function(beta, truth, total, block = "block1", model = "glmnet", cause = 1L) {

  checkmate::assert_numeric(beta)
  checkmate::assert_list(truth, types = "integer")
  checkmate::assert_list(total, types = "integer")
  checkmate::assert_choice(block, choices = names(truth), null.ok = FALSE)
  checkmate::assert_choice(model, choices = c("glmnet", "cooper", "rfsrc", "coxboost"), null.ok = FALSE)
  checkmate::assert_choice(cause, choices = c(1L, 2L))

  # Get indices of nonzero/zero coefs, need to offset with index of block
  # e.g. block 2 starts at index 251, but which() returns indices starting at 1 -> add 251-1
  # This feels uncomfortably hacky.
  predicted <- which(beta[total[[block]]] != 0) + (total[[block]][[1]] - 1)
  predicted0 <- which(beta[total[[block]]] == 0) + (total[[block]][[1]] - 1)

  # block31 has only effect on cause 1 -> true effect index set is empty for cause 2 + vice versa
  if (cause == 1L & block == "block32") truth[["block32"]] <- integer(0)
  if (cause == 2L & block == "block31") truth[["block31"]] <- integer(0)

  # block4 has no true effects, list contains integer 0 for reasons, so this works
  total_pos <- length(truth[[block]])
  total_neg <- length(total[[block]]) - total_pos

  # intersection of predicted nonzeros & true effects gives true positives
  tp <- length(intersect(predicted, truth[[block]]))
  # predicted nonzeros, excluding true effects, leaves false positives
  fp <- length(setdiff(predicted, truth[[block]]))
  # predicted zeros, excluding true effects, leaves true negatives
  tn <- length(setdiff(predicted0, truth[[block]]))
  # intersection of predicted zeros and true effects gives false negatives
  fn <- length(intersect(predicted0, truth[[block]]))

  # Check that total matches
  checkmate::assert_true(sum(tp, fp, tn, fn) == length(total[[block]]))
  checkmate::assert_true(total_pos + total_neg == length(total[[block]]))
  checkmate::assert_true(tp + fn == total_pos)
  checkmate::assert_true(fp + tn == total_neg)
  # For cause 1 we can't have TPs / FNs in block32 and vice versa
  if ((cause == 1L & block == "block32") | (cause == 2L & block == "block31")) {
    checkmate::assert_true(tp == 0 & fn == 0)
  }
  # For block 4 and remaining noise vars there should also be no true effects
  if (block %in% c("block4", "noise")) {
    checkmate::assert_true(tp == 0 & fn == 0)
  }

  data.frame(
    model = model,
    cause = cause,
    block = block,
    tp = tp, fp = fp, tn = tn, fn = fn,
    total_pos = total_pos,
    total_neg = total_neg,
    total = length(total[[block]])
  )
}


