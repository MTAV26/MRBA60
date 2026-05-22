# =========================
# UNCERTAINTY MONTHLY MAPS
# Proyección Robinson, sin Antártida
# PDF + JPEG
# =========================

library(ncdf4)
library(maps)
library(sf)
library(terra)

# ==========================================================
# RUTAS
# ==========================================================

# ==== rutas ====
base_dir <- "/mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/"
out_dir  <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS11_R1_unc_Robinson"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(dir.exists(base_dir))

# ==========================================================
# PROYECCIÓN Y DOMINIO
# ==========================================================

crs_ll  <- "EPSG:4326"
crs_rob <- "ESRI:54030"

# Excluir Antártida
lat_min_plot <- -60
lat_max_plot <-  90

# ==========================================================
# PALETA E INTERVALOS
# ==========================================================

cols_unc <- c(
  "#1E3A8A", # navy
  "#2F5FD9", # azul
  "#7E38E0", # morado
  "#FF5CC7", # rosa
  "#F23B43", # rojo
  "#F95B2C", # naranja-rojizo
  "#FEA319", # naranja
  "#FFD24D", # amarillo
  "#FFE680", # amarillo claro
  "#FFF5B3"  # amarillo pálido
)

labels_unc <- c(
  "< 1.5", "1.5–3.0", "3.0–4.5", "4.5–6.0",
  "6.0–7.5", "7.5–9.0", "9.0–10.5", "10.5–12.0",
  "12.0–13.5", "> 13.5"
)

breaks_unc_fixed <- c(
  0, 1.5, 3.0, 4.5, 6.0, 7.5,
  9.0, 10.5, 12.0, 13.5, Inf
)

eng_month <- c(
  "January", "February", "March", "April", "May", "June",
  "July", "August", "September", "October", "November", "December"
)

# ==========================================================
# FUNCIONES AUXILIARES
# ==========================================================

order_month_files <- function(files) {
  
  b <- basename(files)
  
  mm <- suppressWarnings(as.integer(substr(b, 5, 6)))
  
  if (all(!is.na(mm)) && all(mm >= 1 & mm <= 12)) {
    return(files[order(mm, b)])
  }
  
  m2 <- suppressWarnings(
    as.integer(regmatches(b, regexpr("(0[1-9]|1[0-2])", b)))
  )
  
  if (all(!is.na(m2))) {
    return(files[order(m2, b)])
  }
  
  files[order(b)]
}

get_nc_var <- function(nc, candidates) {
  
  for (nm in candidates) {
    if (nm %in% names(nc$var)) {
      return(ncdf4::ncvar_get(nc, nm))
    }
  }
  
  stop(
    "Ninguna de las variables existe en el NetCDF: ",
    paste(candidates, collapse = ", ")
  )
}

matrix_lonlat_to_raster <- function(mat, lon, lat, crs = "EPSG:4326") {
  
  # Asegurar lon creciente
  ord_lon <- order(lon)
  lon <- lon[ord_lon]
  mat <- mat[ord_lon, , drop = FALSE]
  
  # Asegurar lat creciente
  ord_lat <- order(lat)
  lat <- lat[ord_lat]
  mat <- mat[, ord_lat, drop = FALSE]
  
  dlon <- median(diff(lon), na.rm = TRUE)
  dlat <- median(diff(lat), na.rm = TRUE)
  
  xmin <- min(lon, na.rm = TRUE) - dlon / 2
  xmax <- max(lon, na.rm = TRUE) + dlon / 2
  ymin <- min(lat, na.rm = TRUE) - dlat / 2
  ymax <- max(lat, na.rm = TRUE) + dlat / 2
  
  # terra espera matriz [lat, lon], con primera fila al norte
  mat_lat_lon <- t(mat)
  mat_lat_lon <- mat_lat_lon[nrow(mat_lat_lon):1, , drop = FALSE]
  
  r <- terra::rast(
    mat_lat_lon,
    extent = terra::ext(xmin, xmax, ymin, ymax),
    crs = crs
  )
  
  return(r)
}

make_line_ll <- function(lon_vec, lat_vec) {
  sf::st_linestring(cbind(lon_vec, lat_vec))
}

make_parallel <- function(lat0, n = 1000) {
  lon_seq <- seq(-180, 180, length.out = n)
  sf::st_linestring(cbind(lon_seq, rep(lat0, n)))
}

# ==========================================================
# COSTAS / FRONTERAS DEL MUNDO, SIN ANTÁRTIDA
# ==========================================================

world_map <- maps::map(
  "world",
  plot = FALSE,
  fill = TRUE
)

world_sf <- sf::st_as_sf(world_map)
sf::st_crs(world_sf) <- crs_ll

crop_poly_ll <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = -180,
      xmax =  180,
      ymin = lat_min_plot,
      ymax = lat_max_plot
    ),
    crs = sf::st_crs(crs_ll)
  )
)

world_sf_crop <- suppressWarnings(
  sf::st_intersection(
    sf::st_make_valid(world_sf),
    crop_poly_ll
  )
)

world_rob <- sf::st_transform(world_sf_crop, crs_rob)

# ==========================================================
# BORDE EXTERIOR DEL DOMINIO
# Meridianos -180 / 180 y paralelos -60 / 90
# Densificados antes de proyectar para evitar diagonales
# ==========================================================

n_edge <- 1000

edge_bottom <- make_line_ll(
  lon_vec = seq(-180, 180, length.out = n_edge),
  lat_vec = rep(lat_min_plot, n_edge)
)

edge_top <- make_line_ll(
  lon_vec = seq(-180, 180, length.out = n_edge),
  lat_vec = rep(lat_max_plot, n_edge)
)

edge_left <- make_line_ll(
  lon_vec = rep(-180, n_edge),
  lat_vec = seq(lat_min_plot, lat_max_plot, length.out = n_edge)
)

edge_right <- make_line_ll(
  lon_vec = rep(180, n_edge),
  lat_vec = seq(lat_min_plot, lat_max_plot, length.out = n_edge)
)

domain_lines_ll <- sf::st_sfc(
  edge_bottom,
  edge_top,
  edge_left,
  edge_right,
  crs = crs_ll
)

domain_lines_rob <- sf::st_transform(domain_lines_ll, crs_rob)

domain_bbox <- sf::st_bbox(domain_lines_rob)

xlim_rob <- c(domain_bbox["xmin"], domain_bbox["xmax"])
ylim_rob <- c(domain_bbox["ymin"], domain_bbox["ymax"])

# ==========================================================
# PARALELOS: 40S, 0, 40N
# ==========================================================

parallels_ll <- sf::st_sfc(
  make_parallel(-40),
  make_parallel(0),
  make_parallel(40),
  crs = crs_ll
)

parallels_rob <- sf::st_transform(parallels_ll, crs_rob)

# ==========================================================
# LOCALIZAR SUBCARPETAS DE AÑOS YYYY
# ==========================================================

cand1 <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
cand1 <- cand1[file.info(cand1)$isdir]

year_dirs <- cand1[grepl("[/\\][0-9]{4}$", cand1)]

if (length(year_dirs) == 0) {
  cand2 <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
  year_dirs <- cand2[grepl("[/\\][0-9]{4}$", cand2)]
}

if (length(year_dirs) == 0) {
  stop("No encontré carpetas de años dentro de: ", base_dir)
}

years <- as.integer(sub(".*[/\\\\]([0-9]{4})$", "\\1", year_dirs))

ord_y <- order(years)
year_dirs <- year_dirs[ord_y]
years     <- years[ord_y]

message("Años detectados: ", paste(years, collapse = ", "))

# ==========================================================
# BUCLE PRINCIPAL POR AÑO Y MES
# ==========================================================

for (i in seq_along(years)) {
  
  y <- years[i]
  year_path <- year_dirs[i]
  
  files_y <- list.files(
    year_path,
    pattern = "\\.nc$",
    full.names = TRUE
  )
  
  if (length(files_y) == 0) {
    message("Sin .nc para el año ", y, " en ", year_path)
    next
  }
  
  files_y <- order_month_files(files_y)
  
  for (f in files_y) {
    
    bname <- basename(f)
    
    # ------------------------------------------------------
    # Obtener mes desde el nombre del archivo
    # ------------------------------------------------------
    
    m_chr <- substr(bname, 5, 6)
    m_int <- suppressWarnings(as.integer(m_chr))
    
    if (is.na(m_int) || m_int < 1 || m_int > 12) {
      m_alt <- suppressWarnings(
        as.integer(regmatches(bname, regexpr("(0[1-9]|1[0-2])", bname)))
      )
      if (!is.na(m_alt)) m_int <- m_alt
    }
    
    if (is.na(m_int)) {
      m_int <- ((which(files_y == f) - 1) %% 12) + 1
    }
    
    m <- sprintf("%02d", m_int)
    month_lab <- paste(eng_month[m_int], y)
    
    message("Procesando: ", bname, " — ", month_lab)
    
    # ------------------------------------------------------
    # Leer NetCDF y variable de incertidumbre
    # ------------------------------------------------------
    
    nc <- ncdf4::nc_open(f)
    
    lon <- ncdf4::ncvar_get(nc, "lon")
    lat <- ncdf4::ncvar_get(nc, "lat")
    
    unc <- get_nc_var(nc, c("uncertainty"))
    
    ncdf4::nc_close(nc)
    
    # ------------------------------------------------------
    # Orientar a [lon, lat] si viene [lat, lon]
    # ------------------------------------------------------
    
    if (nrow(unc) == length(lat) && ncol(unc) == length(lon)) {
      unc <- t(unc)
    }
    
    if (!(nrow(unc) == length(lon) && ncol(unc) == length(lat))) {
      stop(
        "Dimensiones no compatibles en archivo: ", f,
        "\nDim unc: ", paste(dim(unc), collapse = " x "),
        "\nlon: ", length(lon),
        "\nlat: ", length(lat)
      )
    }
    
    # ------------------------------------------------------
    # Convertir a km² y ocultar ceros
    # ------------------------------------------------------
    
    unc_km2 <- unc / 1e6
    unc_km2[unc_km2 == 0] <- NA
    
    if (!any(is.finite(unc_km2))) {
      message("Todo NA en ", y, "-", m, ". Omito el mapa.")
      next
    }
    
    # ------------------------------------------------------
    # Crear raster lon/lat
    # ------------------------------------------------------
    
    r_ll <- matrix_lonlat_to_raster(
      mat = unc_km2,
      lon = lon,
      lat = lat,
      crs = crs_ll
    )
    
    # ------------------------------------------------------
    # Recortar Antártida antes de proyectar
    # ------------------------------------------------------
    
    r_ll_crop <- terra::crop(
      r_ll,
      terra::ext(-180, 180, lat_min_plot, lat_max_plot)
    )
    
    # ------------------------------------------------------
    # Reproyectar a Robinson
    # method = near para mantener clases/valores sin suavizado
    # ------------------------------------------------------
    
    r_rob <- terra::project(
      r_ll_crop,
      crs_rob,
      method = "near"
    )
    
    # ======================================================
    # Función interna para dibujar el mapa en Robinson
    # ======================================================
    
    plot_monthly_unc_robinson <- function() {
      
      par(
        mar = c(0.8, 0.8, 1.2, 0.8),
        xaxs = "i",
        yaxs = "i"
      )
      
      plot(
        r_rob,
        col = cols_unc,
        breaks = breaks_unc_fixed,
        legend = FALSE,
        axes = FALSE,
        box = FALSE,
        main = "",
        asp = NA,
        maxcell = Inf,
        xlim = xlim_rob,
        ylim = ylim_rob
      )
      
      # Paralelos: 40S, 0, 40N
      plot(
        sf::st_geometry(parallels_rob),
        add = TRUE,
        col = "grey75",
        lwd = 0.5,
        lty = 1
      )
      
      # Costas / fronteras del mundo
      plot(
        sf::st_geometry(world_rob),
        add = TRUE,
        border = "grey35",
        lwd = 0.35,
        col = NA
      )
      
      # Borde exterior del dominio Robinson
      plot(
        sf::st_geometry(domain_lines_rob),
        add = TRUE,
        col = "grey45",
        lwd = 0.55
      )
      
      # No usar box(), para evitar recuadro negro externo
      
      # ----------------------------------------------------
      # Leyenda interna, 2 columnas de 5, bajos arriba
      # ----------------------------------------------------
      
      usr <- par("usr")
      xr  <- diff(usr[1:2])
      yr  <- diff(usr[3:4])
      
      x0 <- usr[1] + 0.055 * xr
      y0 <- usr[3] + 0.115 * yr
      
      dx  <- 0.030 * xr
      dy  <- 0.030 * yr
      pad <- dy * 0.45
      
      gap_cols <- 0.025 * xr
      txt_pad  <- 0.006 * xr
      
      idx1 <- 1:5
      idx2 <- 6:10
      
      tw   <- max(strwidth(labels_unc, cex = 0.90))
      w1   <- dx + txt_pad + tw
      w2   <- dx + txt_pad + tw
      wTot <- w1 + gap_cols + w2
      
      xC <- x0 + wTot / 2
      
      yTop <- y0 +
        max(length(idx1), length(idx2)) * (dy + pad) +
        pad * 1.2
      
      title_gap    <- dy * 1.0
      subtitle_gap <- dy * 0.6
      
      rect(
        x0 - 0.008 * xr,
        y0 - pad * 0.7,
        x0 + wTot + 0.008 * xr,
        yTop + pad * 0.7 + title_gap + subtitle_gap + dy,
        col = adjustcolor("white", 0.88),
        border = NA
      )
      
      sub_txt <- if (y <= 2018) {
        expression("Root mean square error (km"^2*")")
      } else {
        expression("Standard error (km"^2*")")
      }
      
      text(
        xC,
        yTop + subtitle_gap + title_gap,
        month_lab,
        cex = 1.05,
        font = 2
      )
      
      text(
        xC,
        yTop + subtitle_gap,
        sub_txt,
        cex = 0.95
      )
      
      # Columna 1
      n1 <- length(idx1)
      top1 <- y0 + n1 * (dy + pad) - dy
      
      for (k in seq_along(idx1)) {
        
        ii <- idx1[k]
        yy <- top1 - (k - 1) * (dy + pad)
        
        rect(
          x0, yy,
          x0 + dx, yy + dy,
          col = cols_unc[ii],
          border = "black",
          lwd = 0.5
        )
        
        text(
          x0 + dx + txt_pad,
          yy + dy / 2,
          labels_unc[ii],
          adj = c(0, 0.5),
          cex = 0.90
        )
      }
      
      # Columna 2
      n2 <- length(idx2)
      top2 <- y0 + n2 * (dy + pad) - dy
      
      x1 <- x0 + w1 + gap_cols
      
      for (k in seq_along(idx2)) {
        
        ii <- idx2[k]
        yy <- top2 - (k - 1) * (dy + pad)
        
        rect(
          x1, yy,
          x1 + dx, yy + dy,
          col = cols_unc[ii],
          border = "black",
          lwd = 0.5
        )
        
        text(
          x1 + dx + txt_pad,
          yy + dy / 2,
          labels_unc[ii],
          adj = c(0, 0.5),
          cex = 0.90
        )
      }
    }
    
    # ======================================================
    # Guardar PDF
    # ======================================================
    
    # pdf_file <- file.path(out_dir, sprintf("UNC_%d-%s_Robinson.pdf", y, m))
    # 
    # pdf(
    #   file = pdf_file,
    #   width = 11,
    #   height = 6.5,
    #   useDingbats = FALSE
    # )
    # 
    # plot_monthly_unc_robinson()
    # dev.off()
    # 
    # ======================================================
    # Guardar JPEG
    # ======================================================
    
    jpeg_file <- file.path(out_dir, sprintf("UNC_%d-%s_Robinson.jpeg", y, m))
    
    jpeg(
      filename = jpeg_file,
      width = 13,
      height = 6.5,
      units = "in",
      res = 300,
      quality = 95
    )
    
    plot_monthly_unc_robinson()
    dev.off()
    
    message("PDF guardado:  ", pdf_file)
    message("JPEG guardado: ", jpeg_file)
  }
}

message("Listo. PDFs y JPEGs en proyección Robinson guardados en: ", out_dir)