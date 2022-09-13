# Variable selection ----------------------------------------------------------------------------------------------

# fwelnet wrapper for variable selection sim
fwel_mt_varselect_wrapper <- function(
    data, job, instance,
    alpha = 1, z_method = "original",
    mt_max_iter = 2,
    t = 100,
    a = 0.5,
    thresh = 1e-7
) {

  checkmate::assert_numeric(alpha)
  checkmate::assert_int(mt_max_iter)
  checkmate::assert_int(t)
  checkmate::assert_numeric(a)
  checkmate::assert_numeric(thresh)

  # instance <- sim_surv_binder(n = 400, p = 5000)

  fit <- fwelnet::fwelnet_mt_cox(
    instance$data,
    mt_max_iter = mt_max_iter,
    z_method = as.character(z_method),
    alpha = alpha,
    t = t,
    a = a,
    thresh = thresh
  )

  p <- ncol(instance$data)
  truth <- instance$covar_true_effect
  total <- instance$covar_blocks

  glmnet_beta1 <- fit$beta1[, 1]
  glmnet_beta2 <- fit$beta2[, 1]

  # Same for fwelnet estimates
  fwel_beta1 <- fit$beta1[, mt_max_iter + 1]
  fwel_beta2 <- fit$beta2[, mt_max_iter + 1]

  # Get confusion matrix per block of predictors, by model, by cause.
  # I am painfully aware that this is not "nice", but time is finite and patience is, too.

  rbind(
    # glmnet, cause 1
    get_confusion(glmnet_beta1, truth, total, "block1",  model = "glmnet", cause = 1L),
    get_confusion(glmnet_beta1, truth, total, "block2",  model = "glmnet", cause = 1L),
    get_confusion(glmnet_beta1, truth, total, "block31", model = "glmnet", cause = 1L),
    get_confusion(glmnet_beta1, truth, total, "block32", model = "glmnet", cause = 1L),
    get_confusion(glmnet_beta1, truth, total, "block4",  model = "glmnet", cause = 1L),
    # Cause 2
    get_confusion(glmnet_beta2, truth, total, "block1",  model = "glmnet", cause = 2L),
    get_confusion(glmnet_beta2, truth, total, "block2",  model = "glmnet", cause = 2L),
    get_confusion(glmnet_beta2, truth, total, "block31", model = "glmnet", cause = 2L),
    get_confusion(glmnet_beta2, truth, total, "block32", model = "glmnet", cause = 2L),
    get_confusion(glmnet_beta2, truth, total, "block4",  model = "glmnet", cause = 2L),
    # fwelnet, cause 1
    get_confusion(fwel_beta1, truth, total, "block1",  model = "fwelnet", cause = 1L),
    get_confusion(fwel_beta1, truth, total, "block2",  model = "fwelnet", cause = 1L),
    get_confusion(fwel_beta1, truth, total, "block31", model = "fwelnet", cause = 1L),
    get_confusion(fwel_beta1, truth, total, "block32", model = "fwelnet", cause = 1L),
    get_confusion(fwel_beta1, truth, total, "block4",  model = "fwelnet", cause = 1L),
    # Cause 2
    get_confusion(fwel_beta2, truth, total, "block1",  model = "fwelnet", cause = 2L),
    get_confusion(fwel_beta2, truth, total, "block2",  model = "fwelnet", cause = 2L),
    get_confusion(fwel_beta2, truth, total, "block31", model = "fwelnet", cause = 2L),
    get_confusion(fwel_beta2, truth, total, "block32", model = "fwelnet", cause = 2L),
    get_confusion(fwel_beta2, truth, total, "block4",  model = "fwelnet", cause = 2L)
  )

}

get_confusion <- function(beta, truth, total, block = "block1", model = "glmnet", cause = 1L) {

  checkmate::assert_numeric(beta)
  checkmate::assert_list(truth, types = "integer")
  checkmate::assert_list(total, types = "integer")
  checkmate::assert_choice(block, choices = names(truth), null.ok = FALSE)
  checkmate::assert_choice(model, choices = c("glmnet", "fwelnet"), null.ok = FALSE)
  checkmate::assert_choice(cause, choices = c(1L, 2L))

  # Get indices of nonzero/zero coefs, need to offset with index of block
  # e.g. block 2 starts at index 251, but which() returns indices starting at 1 -> add 251-1
  # This feels uncomfortably hacky.
  predicted <- which(beta[total[[block]]] != 0) + (total[[block]][[1]] - 1)
  predicted0 <- which(beta[total[[block]]] == 0) + (total[[block]][[1]] - 1)

  # Double check that indices are set correctly. Names are x<j>, so we can check that quickly
  checkmate::assert(
    checkmate::assert_true(
      all.equal(as.integer(sub(pattern = "x", replacement = "", names(predicted))), unname(predicted))
    ),
    checkmate::assert_true(
      all.equal(as.integer(sub(pattern = "x", replacement = "", names(predicted0))), unname(predicted0))
    )
  )
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


