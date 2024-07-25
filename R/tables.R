#' Quick wrapper to collect variable selection scores in a table for LaTeX.
#' @param res_long A data.table with `model`, `cause`, `block`, `value` columns.
#' @param measure The measure to extract, defaults to `"PPV"`.
#' @param aggr_point,aggr_var Functions to aggregate point and variance estimates, defaults to `median` and `IQR`.
#' @param minmax Function to determine if a value is a minimum or maximum, defaults to `max`, indicating that higher is better.
#' @param kbl_format Format for `kableExtra::cell_spec()`, defaults to `"latex"`.
measure_table <- function(res_long, measure = "PPV", aggr_point = median, aggr_var = IQR, minmax = max, kbl_format = "latex") {
  res_long |>
    dplyr::filter(lambda2 == lambda1) |>
    dplyr::filter(.data$measure %in% .env$measure) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::select(model, cause, block, value) |>
    dplyr::group_by(model, cause, block) |>
    dplyr::summarize(
      aggr = 100 * aggr_point(value),
      aggrv = 100 * aggr_var(value),
      .groups = "drop"
    ) |>
    dplyr::group_by(cause, block) |>
    dplyr::mutate(
      stat_formatted = glue::glue("{round(aggr, 2)} ({round(aggrv, 2)})"),
      stat_formatted = kableExtra::cell_spec(stat_formatted, bold = aggr == minmax(aggr), format = kbl_format)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-aggr, -aggrv) |>
    tidyr::pivot_wider(names_from = model, values_from = stat_formatted) |>
    dplyr::select(Cause = cause, Block = block, CooPeR, Coxnet, CoxBoost, RSF)
}
