get_days_in_month <- function(date) {
  next_month <- seq(date, length = 2, by = "month")[2]
  as.integer(format(next_month - 1, "%d"))
}