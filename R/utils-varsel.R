#' Get selected variables for each model type
selected <- function(x, ...) {
  UseMethod("selected")
}

selected.cooper <- function(x, model = "cooper", ...) {
  checkmate::assert_list(x, any.missing = FALSE)
  checkmate::assert_subset(model, choices = c("cooper", "coxnet"))

  use_initial_fit <- model == "coxnet"

  list(
    `1` = nonzeros(coef(x, event = 1, use_initial_fit = use_initial_fit)),
    `2` = nonzeros(coef(x, event = 2, use_initial_fit = use_initial_fit))
  )
}

#' Quikcly accessing selected variables base on different model objects
#' Focusing on cause-specific models, so only one cause-specific
#' set of selected variables for each function call.
#' @param x A model object, `coxboost` `rfsrc`.
#' @return A list of cause-specific selected variables.
selected.rfsrc <- function(x, ...) {
  checkmate::assert_integerish(x$cause.wt, lower = 0, upper = 1, len = 2)
  # This is c(1, 0) for cause 1 and c(0, 1) for cause 2
  # so we get the position of the 1
  cause <- which(x$cause.wt == 1)

  vimp <- x[["importance"]][, cause]
  vimp <- vimp[vimp > abs(min(vimp))]

  res = list(vimp)
  names(res) = cause
  res
}

selected.CoxBoost <- function(x, ...) {
  list(
    `1` = nonzeros(coef(x)[[1]]),
    `2` = nonzeros(coef(x)[[2]])
  )
}


#' Utility function to access non-zero elements in a vector.
nonzeros <- function(x) {
  res <- x[x != 0]
  checkmate::assert_numeric(res, min.len = 1, finite = TRUE, any.missing = FALSE)
  res
}
