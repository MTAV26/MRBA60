# ==========================================================
# 18) TREND AND SUMMARY STATISTICS FOR FIGURE TEXT
# ==========================================================
# Requires:
#   - df_monthly
#   - df_annual
#   - df_seasonal
#
# Outputs:
#   - Annual_BA_stats_MRBA60_Figure7.csv
#   - Seasonal_BA_stats_MRBA60_Figure7.csv
#   - Annual_and_Seasonal_BA_stats_MRBA60_Figure7.csv
#   - Figure7_text_numbers_summary.txt
#   - Figure7_BA_stats_MRBA60.xlsx
# ==========================================================


# ==========================================================
# 18.1) Load modified Mann-Kendall + Sen functions
# ==========================================================

mmkh_file <- "/home/miguel/Escritorio/Codigospaper/mmkh.R"
sen_file  <- "/home/miguel/Escritorio/Codigospaper/sen.R"

if (file.exists(mmkh_file)) {
  source(mmkh_file)
} else {
  stop("mmkh.R file not found: ", mmkh_file)
}

if (file.exists(sen_file)) {
  source(sen_file)
} else {
  warning("sen.R file not found: ", sen_file)
}

if (!exists("mmkh_v2")) {
  stop("Function mmkh_v2() is not available after sourcing mmkh.R.")
}

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)


# ==========================================================
# 18.2) Helper functions
# ==========================================================

pct_over_period <- function(slope, mean_y, x_min, x_max) {
  if (!is.finite(slope) || !is.finite(mean_y) || mean_y == 0) return(NA_real_)
  n_years <- x_max - x_min
  100 * slope * n_years / mean_y
}

extract_pval <- function(est) {
  
  if (is.null(est)) return(NA_real_)
  
  nm <- names(est)
  if (is.null(nm)) return(NA_real_)
  
  nm_low <- tolower(nm)
  
  cand_exact <- which(nm_low %in% c(
    "pvalue", "p-value", "p.value", "p_val", "p", "sl", "sig", "significance"
  ))
  
  if (length(cand_exact) > 0) {
    return(as.numeric(est[cand_exact[1]]))
  }
  
  cand_grep <- grep(
    "p\\s*[-._]??\\s*value|^p$|p\\.value|significance|sl",
    nm,
    ignore.case = TRUE
  )
  
  if (length(cand_grep) > 0) {
    return(as.numeric(est[cand_grep[1]]))
  }
  
  return(NA_real_)
}

extract_sen_slope <- function(est) {
  
  if (is.null(est)) return(NA_real_)
  
  nm <- names(est)
  if (is.null(nm)) return(NA_real_)
  
  if ("Sen's slope" %in% nm) {
    return(as.numeric(est["Sen's slope"]))
  }
  
  cand <- grep("sen.*slope|slope.*sen", nm, ignore.case = TRUE)
  
  if (length(cand) > 0) {
    return(as.numeric(est[cand[1]]))
  }
  
  return(NA_real_)
}

sig_star <- function(pval) {
  if (!is.finite(pval)) return("")
  if (pval < 0.01) return("***")
  if (pval < 0.05) return("**")
  if (pval < 0.10) return("*")
  return("")
}

product_period <- function(product, group_name) {
  
  # En el script actual todos los productos cubren 2003-2024.
  # Para DJF, la serie empieza en 2004 porque DJF 2003 está incompleta.
  
  start_year <- ifelse(group_name == "DJF", 2004, 2003)
  end_year   <- 2024
  
  c(start_year, end_year)
}


# ==========================================================
# 18.3) Generic statistics for one annual/seasonal series
# ==========================================================

stats_one_series <- function(df_group, group_name, product_name) {
  
  period <- product_period(product_name, group_name)
  start_year <- period[1]
  end_year   <- period[2]
  
  df_sub <- df_group %>%
    dplyr::filter(Year >= start_year, Year <= end_year) %>%
    dplyr::select(Year, dplyr::all_of(product_name)) %>%
    dplyr::rename(BA_km2 = dplyr::all_of(product_name)) %>%
    dplyr::filter(is.finite(Year), is.finite(BA_km2)) %>%
    dplyr::arrange(Year)
  
  if (nrow(df_sub) < 3) {
    return(data.frame(
      Group = group_name,
      Product = product_name,
      Period = paste0(start_year, "-", end_year),
      N_years = nrow(df_sub),
      Mean_BA_Mkm2 = NA_real_,
      Sen_slope_Mkm2_per_year = NA_real_,
      Sen_slope_pct_over_period = NA_real_,
      MK_pvalue = NA_real_,
      MK_sig_0.05 = NA,
      Year_max = NA_integer_,
      Max_BA_Mkm2 = NA_real_,
      Year_min = NA_integer_,
      Min_BA_Mkm2 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  y <- df_sub$BA_km2
  yrs <- df_sub$Year
  
  est <- mmkh_v2(y)
  
  slope_km2_per_year <- extract_sen_slope(est)
  pval <- extract_pval(est)
  
  mean_km2 <- mean(y, na.rm = TRUE)
  
  pct_period <- pct_over_period(
    slope = slope_km2_per_year,
    mean_y = mean_km2,
    x_min = min(yrs),
    x_max = max(yrs)
  )
  
  idx_max <- which.max(y)
  idx_min <- which.min(y)
  
  data.frame(
    Group = group_name,
    Product = product_name,
    Period = paste0(min(yrs), "-", max(yrs)),
    N_years = length(y),
    Mean_BA_Mkm2 = mean_km2 / 1e6,
    Sen_slope_Mkm2_per_year = slope_km2_per_year / 1e6,
    Sen_slope_pct_over_period = pct_period,
    MK_pvalue = pval,
    MK_sig_0.05 = ifelse(is.finite(pval), pval <= 0.05, NA),
    Year_max = yrs[idx_max],
    Max_BA_Mkm2 = y[idx_max] / 1e6,
    Year_min = yrs[idx_min],
    Min_BA_Mkm2 = y[idx_min] / 1e6,
    stringsAsFactors = FALSE
  )
}


# ==========================================================
# 18.4) Difference statistics: MRBA60 - FireCCI51
# ==========================================================

diff_stats_one_group <- function(df_group, group_name) {
  
  start_year <- ifelse(group_name == "DJF", 2004, 2003)
  end_year   <- 2024
  
  df_sub <- df_group %>%
    dplyr::filter(Year >= start_year, Year <= end_year) %>%
    dplyr::mutate(
      Diff_km2 = MRBA60 - FireCCI51
    ) %>%
    dplyr::filter(is.finite(Diff_km2), is.finite(FireCCI51)) %>%
    dplyr::arrange(Year)
  
  if (nrow(df_sub) == 0) {
    return(data.frame(
      Group = group_name,
      Diff_period = paste0(start_year, "-", end_year),
      Sum_Diff_MRBA60_minus_FireCCI51_Mkm2 = NA_real_,
      Mean_Diff_MRBA60_minus_FireCCI51_Mkm2 = NA_real_,
      Pct_increase_vs_FireCCI51 = NA_real_,
      Year_max_positive_diff = NA_integer_,
      Max_positive_diff_Mkm2 = NA_real_,
      Year_min_diff = NA_integer_,
      Min_diff_Mkm2 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  sum_diff_km2 <- sum(df_sub$Diff_km2, na.rm = TRUE)
  mean_diff_km2 <- mean(df_sub$Diff_km2, na.rm = TRUE)
  sum_f51_km2 <- sum(df_sub$FireCCI51, na.rm = TRUE)
  
  pct_inc <- ifelse(
    is.finite(sum_f51_km2) && sum_f51_km2 != 0,
    100 * sum_diff_km2 / sum_f51_km2,
    NA_real_
  )
  
  idx_max <- which.max(df_sub$Diff_km2)
  idx_min <- which.min(df_sub$Diff_km2)
  
  data.frame(
    Group = group_name,
    Diff_period = paste0(min(df_sub$Year), "-", max(df_sub$Year)),
    Sum_Diff_MRBA60_minus_FireCCI51_Mkm2 = sum_diff_km2 / 1e6,
    Mean_Diff_MRBA60_minus_FireCCI51_Mkm2 = mean_diff_km2 / 1e6,
    Pct_increase_vs_FireCCI51 = pct_inc,
    Year_max_positive_diff = df_sub$Year[idx_max],
    Max_positive_diff_Mkm2 = df_sub$Diff_km2[idx_max] / 1e6,
    Year_min_diff = df_sub$Year[idx_min],
    Min_diff_Mkm2 = df_sub$Diff_km2[idx_min] / 1e6,
    stringsAsFactors = FALSE
  )
}


# ==========================================================
# 18.5) Annual table
# ==========================================================

products <- c("MRBA60", "FireCCI51", "MCD64A1", "GFED5")

df_annual_stats <- do.call(
  rbind,
  lapply(products, function(pp) {
    stats_one_series(
      df_group = df_annual,
      group_name = "ANNUAL",
      product_name = pp
    )
  })
)

df_annual_diff <- diff_stats_one_group(
  df_group = df_annual,
  group_name = "ANNUAL"
)

df_annual_stats <- df_annual_stats %>%
  dplyr::left_join(df_annual_diff, by = "Group")


# ==========================================================
# 18.6) Seasonal table
# ==========================================================

season_order <- c("DJF", "MAM", "JJA", "SON")

df_seasonal_stats <- do.call(
  rbind,
  lapply(season_order, function(ss) {
    
    df_ss <- df_seasonal %>%
      dplyr::filter(Period == ss)
    
    do.call(
      rbind,
      lapply(products, function(pp) {
        stats_one_series(
          df_group = df_ss,
          group_name = ss,
          product_name = pp
        )
      })
    )
  })
)

df_seasonal_diff <- do.call(
  rbind,
  lapply(season_order, function(ss) {
    
    df_ss <- df_seasonal %>%
      dplyr::filter(Period == ss)
    
    diff_stats_one_group(
      df_group = df_ss,
      group_name = ss
    )
  })
)

df_seasonal_stats <- df_seasonal_stats %>%
  dplyr::left_join(df_seasonal_diff, by = "Group")


# ==========================================================
# 18.7) Rounding for reporting
# ==========================================================

round_stats <- function(df) {
  
  df %>%
    dplyr::mutate(
      Mean_BA_Mkm2 = round(Mean_BA_Mkm2, 3),
      Sen_slope_Mkm2_per_year = round(Sen_slope_Mkm2_per_year, 4),
      Sen_slope_pct_over_period = round(Sen_slope_pct_over_period, 2),
      MK_pvalue = signif(MK_pvalue, 3),
      Max_BA_Mkm2 = round(Max_BA_Mkm2, 3),
      Min_BA_Mkm2 = round(Min_BA_Mkm2, 3),
      Sum_Diff_MRBA60_minus_FireCCI51_Mkm2 = round(Sum_Diff_MRBA60_minus_FireCCI51_Mkm2, 3),
      Mean_Diff_MRBA60_minus_FireCCI51_Mkm2 = round(Mean_Diff_MRBA60_minus_FireCCI51_Mkm2, 3),
      Pct_increase_vs_FireCCI51 = round(Pct_increase_vs_FireCCI51, 2),
      Max_positive_diff_Mkm2 = round(Max_positive_diff_Mkm2, 3),
      Min_diff_Mkm2 = round(Min_diff_Mkm2, 3)
    )
}

df_annual_stats_out <- round_stats(df_annual_stats)
df_seasonal_stats_out <- round_stats(df_seasonal_stats)

df_all_stats_out <- rbind(
  df_annual_stats_out,
  df_seasonal_stats_out
)


# ==========================================================
# 18.8) Save CSV and XLSX
# ==========================================================

write.csv(
  df_annual_stats_out,
  file = file.path(out_dir, "Annual_BA_stats_MRBA60_Figure7.csv"),
  row.names = FALSE
)

write.csv(
  df_seasonal_stats_out,
  file = file.path(out_dir, "Seasonal_BA_stats_MRBA60_Figure7.csv"),
  row.names = FALSE
)

write.csv(
  df_all_stats_out,
  file = file.path(out_dir, "Annual_and_Seasonal_BA_stats_MRBA60_Figure7.csv"),
  row.names = FALSE
)

xlsx_file <- file.path(out_dir, "Figure7_BA_stats_MRBA60.xlsx")

openxlsx::write.xlsx(
  list(
    Annual = df_annual_stats_out,
    Seasonal = df_seasonal_stats_out,
    Annual_and_Seasonal = df_all_stats_out
  ),
  file = xlsx_file,
  overwrite = TRUE
)


# ==========================================================
# 18.9) Console summary for manuscript text
# ==========================================================

get_stat <- function(df, group, product, column) {
  df %>%
    dplyr::filter(Group == group, Product == product) %>%
    dplyr::pull(dplyr::all_of(column)) %>%
    .[1]
}

get_diff <- function(df, group, column) {
  df %>%
    dplyr::filter(Group == group) %>%
    dplyr::pull(dplyr::all_of(column)) %>%
    .[1]
}

annual_mrba_mean <- get_stat(df_annual_stats_out, "ANNUAL", "MRBA60", "Mean_BA_Mkm2")
annual_f51_mean  <- get_stat(df_annual_stats_out, "ANNUAL", "FireCCI51", "Mean_BA_Mkm2")
annual_mcd_mean  <- get_stat(df_annual_stats_out, "ANNUAL", "MCD64A1", "Mean_BA_Mkm2")
annual_gfed_mean <- get_stat(df_annual_stats_out, "ANNUAL", "GFED5", "Mean_BA_Mkm2")

annual_sum_diff  <- get_diff(df_annual_stats_out, "ANNUAL", "Sum_Diff_MRBA60_minus_FireCCI51_Mkm2")
annual_mean_diff <- get_diff(df_annual_stats_out, "ANNUAL", "Mean_Diff_MRBA60_minus_FireCCI51_Mkm2")
annual_pct_inc   <- get_diff(df_annual_stats_out, "ANNUAL", "Pct_increase_vs_FireCCI51")

annual_mrba_trend <- get_stat(df_annual_stats_out, "ANNUAL", "MRBA60", "Sen_slope_pct_over_period")
annual_f51_trend  <- get_stat(df_annual_stats_out, "ANNUAL", "FireCCI51", "Sen_slope_pct_over_period")
annual_mcd_trend  <- get_stat(df_annual_stats_out, "ANNUAL", "MCD64A1", "Sen_slope_pct_over_period")
annual_gfed_trend <- get_stat(df_annual_stats_out, "ANNUAL", "GFED5", "Sen_slope_pct_over_period")

annual_mrba_p <- get_stat(df_annual_stats_out, "ANNUAL", "MRBA60", "MK_pvalue")
annual_f51_p  <- get_stat(df_annual_stats_out, "ANNUAL", "FireCCI51", "MK_pvalue")
annual_mcd_p  <- get_stat(df_annual_stats_out, "ANNUAL", "MCD64A1", "MK_pvalue")
annual_gfed_p <- get_stat(df_annual_stats_out, "ANNUAL", "GFED5", "MK_pvalue")

summary_lines <- c(
  "==========================================================",
  "FIGURE 7 - VALUES FOR MANUSCRIPT TEXT",
  "==========================================================",
  "",
  "Annual BA:",
  paste0("MRBA60 mean annual BA: ", annual_mrba_mean, " Mkm2"),
  paste0("FireCCI51 mean annual BA: ", annual_f51_mean, " Mkm2"),
  paste0("MCD64A1 mean annual BA: ", annual_mcd_mean, " Mkm2"),
  paste0("GFED5 mean annual BA: ", annual_gfed_mean, " Mkm2"),
  "",
  "Annual MRBA60 - FireCCI51 difference:",
  paste0("Accumulated difference: ", annual_sum_diff, " Mkm2"),
  paste0("Mean annual difference: ", annual_mean_diff, " Mkm2 yr-1"),
  paste0("Increase relative to FireCCI51: ", annual_pct_inc, " %"),
  "",
  "Annual trends, accumulated relative change over period:",
  paste0("MRBA60: ", annual_mrba_trend, " %, p = ", annual_mrba_p),
  paste0("FireCCI51: ", annual_f51_trend, " %, p = ", annual_f51_p),
  paste0("MCD64A1: ", annual_mcd_trend, " %, p = ", annual_mcd_p),
  paste0("GFED5: ", annual_gfed_trend, " %, p = ", annual_gfed_p),
  "",
  "Seasonal summary:"
)

for (ss in season_order) {
  
  mrba_mean <- get_stat(df_seasonal_stats_out, ss, "MRBA60", "Mean_BA_Mkm2")
  f51_mean  <- get_stat(df_seasonal_stats_out, ss, "FireCCI51", "Mean_BA_Mkm2")
  sum_diff  <- get_diff(df_seasonal_stats_out, ss, "Sum_Diff_MRBA60_minus_FireCCI51_Mkm2")
  mean_diff <- get_diff(df_seasonal_stats_out, ss, "Mean_Diff_MRBA60_minus_FireCCI51_Mkm2")
  pct_inc   <- get_diff(df_seasonal_stats_out, ss, "Pct_increase_vs_FireCCI51")
  
  mrba_trend <- get_stat(df_seasonal_stats_out, ss, "MRBA60", "Sen_slope_pct_over_period")
  f51_trend  <- get_stat(df_seasonal_stats_out, ss, "FireCCI51", "Sen_slope_pct_over_period")
  mcd_trend  <- get_stat(df_seasonal_stats_out, ss, "MCD64A1", "Sen_slope_pct_over_period")
  gfd_trend  <- get_stat(df_seasonal_stats_out, ss, "GFED5", "Sen_slope_pct_over_period")
  
  mrba_p <- get_stat(df_seasonal_stats_out, ss, "MRBA60", "MK_pvalue")
  f51_p  <- get_stat(df_seasonal_stats_out, ss, "FireCCI51", "MK_pvalue")
  mcd_p  <- get_stat(df_seasonal_stats_out, ss, "MCD64A1", "MK_pvalue")
  gfd_p  <- get_stat(df_seasonal_stats_out, ss, "GFED5", "MK_pvalue")
  
  summary_lines <- c(
    summary_lines,
    "",
    paste0(ss, ":"),
    paste0("MRBA60 mean BA: ", mrba_mean, " Mkm2 season-1"),
    paste0("FireCCI51 mean BA: ", f51_mean, " Mkm2 season-1"),
    paste0("Accumulated MRBA60 - FireCCI51 difference: ", sum_diff, " Mkm2"),
    paste0("Mean seasonal difference: ", mean_diff, " Mkm2 season-1"),
    paste0("Increase relative to FireCCI51: ", pct_inc, " %"),
    paste0("MRBA60 trend: ", mrba_trend, " %, p = ", mrba_p),
    paste0("FireCCI51 trend: ", f51_trend, " %, p = ", f51_p),
    paste0("MCD64A1 trend: ", mcd_trend, " %, p = ", mcd_p),
    paste0("GFED5 trend: ", gfd_trend, " %, p = ", gfd_p)
  )
}

summary_txt <- file.path(out_dir, "Figure7_text_numbers_summary.txt")

writeLines(summary_lines, con = summary_txt)

cat(paste(summary_lines, collapse = "\n"))
cat("\n\nFiles saved:\n")
cat(file.path(out_dir, "Annual_BA_stats_MRBA60_Figure7.csv"), "\n")
cat(file.path(out_dir, "Seasonal_BA_stats_MRBA60_Figure7.csv"), "\n")
cat(file.path(out_dir, "Annual_and_Seasonal_BA_stats_MRBA60_Figure7.csv"), "\n")
cat(xlsx_file, "\n")
cat(summary_txt, "\n")
