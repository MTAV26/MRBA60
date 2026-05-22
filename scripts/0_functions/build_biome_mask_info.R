
build_biome_mask_info <- function(bioma,
                                  biomas_shp,
                                  lon_range,
                                  lat_range,
                                  lon_mat,
                                  lat_mat,
                                  status_matrix2_tot) {
  bioma_sel <- biomas_shp %>% filter(cont_bm == bioma)
  
  if (nrow(bioma_sel) == 0) return(NULL)
  
  bbox <- st_bbox(bioma_sel)
  
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  
  if (length(lon_idx) == 0 || length(lat_idx) == 0) return(NULL)
  
  lon_mat_crop <- lon_mat[lon_idx, lat_idx, drop = FALSE]
  lat_mat_crop <- lat_mat[lon_idx, lat_idx, drop = FALSE]
  status_crop  <- status_matrix2_tot[lon_idx, lat_idx, , drop = FALSE]
  
  grid_points_crop <- st_as_sf(
    data.frame(
      lon = as.vector(lon_mat_crop),
      lat = as.vector(lat_mat_crop)
    ),
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  inter <- st_intersects(grid_points_crop, biomas_shp)
  
  bioma_asignado <- sapply(inter, function(i) {
    if (length(i) == 0) return("Ninguno")
    biomas_shp$cont_bm[i[1]]
  })
  
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  presencia_en_serie <- apply(status_crop, c(1, 2), function(x) any(x == 1, na.rm = TRUE))
  grid_points_biomas$presencia_serie <- as.vector(presencia_en_serie)
  
  puntos_sin_bioma_con_presencia <- grid_points_biomas %>%
    filter(bioma_final == "Ninguno" & presencia_serie)
  
  puntos_con_bioma <- grid_points_biomas %>%
    filter(bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    nearest_idx <- st_nearest_feature(puntos_sin_bioma_con_presencia, puntos_con_bioma)
    
    idx_replace <- which(grid_points_biomas$bioma_final == "Ninguno" &
                           grid_points_biomas$presencia_serie)
    
    grid_points_biomas$bioma_final[idx_replace] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  mask <- grid_points_biomas$bioma_final == bioma
  mask_matrix <- matrix(
    mask,
    nrow = length(lon_idx),
    ncol = length(lat_idx),
    byrow = FALSE
  )
  
  list(
    bioma = bioma,
    safe_biome_name = sanitize_biome_name(bioma),
    bioma_sel = bioma_sel,
    lon_idx = lon_idx,
    lat_idx = lat_idx,
    lon_mat_crop = lon_mat_crop,
    lat_mat_crop = lat_mat_crop,
    mask_matrix = mask_matrix
  )
}
