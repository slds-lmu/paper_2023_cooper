
#' Cleanup output from riskRegression::Score to make it one long plottable data.frame
#' @param scores As returned by `riskRegression::Score()`, assumed to contain
#'   measures AUC, Brier, IBS and IPA.
cleanup_score <- function(scores) {
  result <- scores$AUC$score[scores$Brier$score, on = c("model", "times")]
  result <- data.table::melt(result , id.vars = c("model", "times"),
                             value.name = "score", variable.name = "metric",
                             meausure.vars = c("AUC", "Brier", "IBS", "IPA"))
  # Exclude some superfluous output
  result <- result[!(is.na(score) & model == "Null model"), ]
  result <- result[!(metric == "IPA" & model == "Null model"), ]
  result
}

#' Refit a GLMnet with riskRegression using the same predictors and lambda as a cooper fit.
#' @param cooperfit A model fit of class `"cooper"`
refit_glmnet <- function(cooperfit, xnew, event = 1) {
  checkmate::assert_class(cooperfit, "cooper")
  checkmate::assert_int(event, lower = 1)

  # Recode to censore other event and ensure status in 0,1
  if (event == 1) {
    xnew$status[xnew$status == 2] <- 0
  }
  if (event == 2) {
    xnew$status[xnew$status == 1] <- 0
    xnew$status[xnew$status == 2] <- 1
  }

  riskRegression::GLMnet(
    formula = stats::reformulate(cooperfit$predictors, response = "Surv(time, status)"),
    data = xnew,
    lambda = cooperfit$initial_fits[[event]]$lambda.min,
    cv = FALSE,
    standardize = TRUE, alpha = 1
  )
}
