plot_common_timeseries <- function(
    bioma,
    safe_biome_name,
    dates_plot,
    ba_ref,
    ba_f51,
    ba_harm,
    suffix,
    title_suffix,
    stats_text = NULL
) {
  
  if (!isTRUE(SAVE_TIMESERIES_PLOTS)) {
    return(invisible(NULL))
  }
  
  max_y <- max(ba_ref, ba_f51, ba_harm, na.rm = TRUE) * 1.2
  if (!is.finite(max_y) || max_y <= 0) {
    max_y <- 1
  }
  
  jpeg_filename <- file.path(
    output_dir_plot,
    paste0(safe_biome_name, "_Time_series_", suffix, ".jpeg")
  )
  
  grDevices::jpeg(
    filename = jpeg_filename,
    width = 2000,
    height = 1200,
    res = 300
  )
  
  par(mar = c(5, 5, 6, 2))
  plot(
    dates_plot,
    ba_ref,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = "",
    ylab = expression("Burned Area (km"^2*")"),
    main = paste(bioma, "\nFireCCIS311, FireCCI51 and Harmonised", title_suffix),
    ylim = c(0, max_y),
    cex.lab = 1.2,
    cex.axis = 1.1,
    cex.main = 1.3,
    xaxt = "n"
  )
  
  lines(dates_plot, ba_f51, col = "orange", lwd = 2)
  lines(dates_plot, ba_harm, col = "darkred", lwd = 2, lty = 1)
  axis.Date(
    1,
    at = seq(min(dates_plot), max(dates_plot), by = "3 months"),
    format = "%b-%Y",
    las = 1,
    cex.axis = 0.5
  )
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  legend(
    "topright",
    legend = c("FireCCIS311", "FireCCI51", "Harmonised"),
    col = c("blue", "orange", "darkred"),
    lwd = c(2, 2, 2),
    lty = c(1, 1, 1),
    ncol = 1,
    bty = "n",
    cex = 0.5,
    xpd = TRUE,
    seg.len = 1,
    text.width = NULL,
    inset = c(0, 0)
  )
  
  if (!is.null(stats_text) && nzchar(stats_text)) {
    usr <- par("usr")
    text(
      x = usr[1],
      y = usr[4],
      labels = stats_text,
      adj = c(0, 1),
      cex = 0.42,
      col = "black",
      lineheight = 1.05
    )
  }
  
  dev.off()
  invisible(jpeg_filename)
}