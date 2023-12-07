library(data.table)
library(fwelnet)

get_pbc_data <- function(job, data, na.rm = TRUE, num_noise = 0) {
  checkmate::assert_int(num_noise, lower = 0)

  ee = new.env(parent = emptyenv(), hash = FALSE)
  data(list = "Pbc3", package = "pec", envir = ee)
  xdf = ee[["Pbc3"]]

  if (na.rm) xdf <- na.omit(xdf)

  xdf$time <- xdf$days
  xdf$days <- NULL

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


# pbc_data <- get_pbc_data(na.rm = TRUE, num_noise = 1)


pbc_varsel_score <- function(job, data, instance, target_event = 1, conf_level = 0.95, mt_max_iter = 3, t = 100, thresh = 1e-7, alpha = 1) {
  pbc_data <- instance

  mod_cph <- survival::coxph(survival::Surv(time, status) ~ .,
                             data = recode_cr_status(pbc_data, target_event = target_event))

  cph_coefs <- broom::tidy(mod_cph)
  coefs_true <- cph_coefs[cph_coefs$p.value < (1 - conf_level), ][["term"]]

  cofit <- fwelnet::cooper(pbc_data, mt_max_iter = 3, stratify_by_status = TRUE, t = 100, thresh = 1e-7)

  pbc_covars = setdiff(names(pbc_data), c("time", "status"))
  noise_vars = pbc_covars[startsWith(pbc_covars, "xnoise")]

  selected_vars <- function(copper_fit, event = 1) {
    res <- coef(copper_fit, event = event)
    names(res[res != 0])
  }

  not_selected_vars <- function(copper_fit, event = 1) {
    res <- coef(copper_fit, event = event)
    names(res[res == 0])
  }

  cooper_selected_c1 <- selected_vars(cofit, event = target_event)
  cooper_not_selected_c1 <- not_selected_vars(cofit, event = target_event)

  which_positive <- pbc_covars %in% coefs_true
  which_negative <- !which_positive

  total_positive <- sum(which_positive)
  total_negative <- sum(which_negative)

  result <- data.table::data.table(
    total = length(pbc_covars),
    total_pos = sum(which_positive),
    total_neg = sum(which_negative),
    tp = sum(cooper_selected_c1 %in% pbc_covars[which_positive]),
    fp = sum(cooper_selected_c1 %in% pbc_covars[which_negative]),
    tn = sum(cooper_not_selected_c1 %in% pbc_covars[which_negative]),
    fn = sum(cooper_not_selected_c1 %in% pbc_covars[which_positive]),
    selected = list(cooper_selected_c1)
  )


  result[, tpr := tp/total_pos]
  result[, fpr := fp/total_neg]
  result[, ppv := tp/(tp + fp)]
  result[, npv := tn/(tn + fn)]
  result[, fdr := 1 - ppv]
  result[, f1 := (2 * ppv * tpr) / (ppv + tpr)]
  result[, acc := (tp+tn)/total]
  result
}
