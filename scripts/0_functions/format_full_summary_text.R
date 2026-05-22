
format_full_summary_text <- function(dates_plot, ba_f51, ba_harm, ba_s3 = NULL, stats_common = NULL) {
  
  fmt_mean <- function(x) {
    if (!any(is.finite(x))) {
      return("NA")
    }
    format_stat_value(mean(x, na.rm = TRUE), digits = 0)
  }
  
  fmt_total <- function(x) {
    if (!any(is.finite(x))) {
      return("NA")
    }
    format_stat_value(sum(x, na.rm = TRUE), digits = 0)
  }
  
  make_period_lines <- function(idx, label) {
    if (!any(idx, na.rm = TRUE)) {
      return(character(0))
    }
    
    has_s3 <- !is.null(ba_s3) && any(is.finite(ba_s3[idx]))
    
    mean_parts <- c(
      paste0("Mean F51 = ", fmt_mean(ba_f51[idx])),
      paste0("Mean Harmonised = ", fmt_mean(ba_harm[idx]))
    )
    
    total_parts <- c(
      paste0("Total F51 = ", fmt_total(ba_f51[idx])),
      paste0("Total Harmonised = ", fmt_total(ba_harm[idx]))
    )
    
    if (has_s3) {
      mean_parts <- c(paste0("Mean S3 = ", fmt_mean(ba_s3[idx])), mean_parts)
      total_parts <- c(paste0("Total S3 = ", fmt_total(ba_s3[idx])), total_parts)
    }
    
    c(
      label,
      paste(mean_parts, collapse = " | "),
      paste(total_parts, collapse = " | ")
    )
  }
  
  historical_idx <- dates_plot < min(dates_common)
  common_idx <- dates_plot >= min(dates_common)
  
  lines <- character(0)
  
  if (!is.null(stats_common) && nrow(stats_common) > 0) {
    lines <- c(lines, format_stats_text(stats_common), "")
  }
  
  lines <- c(lines, "Full-period summary:")
  
  if (any(historical_idx, na.rm = TRUE)) {
    label <- paste0(
      format(min(dates_plot[historical_idx]), "%Y"),
      "-",
      format(max(dates_plot[historical_idx]), "%Y")
    )
    lines <- c(lines, make_period_lines(historical_idx, label))
  }
  
  if (any(common_idx, na.rm = TRUE)) {
    label <- paste0(
      format(min(dates_plot[common_idx]), "%Y"),
      "-",
      format(max(dates_plot[common_idx]), "%Y")
    )
    lines <- c(lines, "", make_period_lines(common_idx, label))
  }
  
  paste(lines, collapse = "\n")
}