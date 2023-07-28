#' Simple status distribution check
#' @param xs A status variable
#' @examples
#'
#' check_status(instance$train$status)
check_status <- function(xs) {
  tab <- as.integer(table(xs))
  mtab <- rbind(events = tab, prop = 100 * round(tab/sum(tab), 2))
  colnames(mtab) <- unique(xs)
  mtab <- cbind(mtab, sum = rowSums(mtab))
  print(mtab)
}
