Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es


README - 4-model
================

This folder contains the modelling scripts used to build the MRBA60 harmonised burned area product. The modelling workflow is divided into two main stages:

  4.1 Predictor pre-selection by biome and month.
  4.2 Final harmonisation model using the selected predictors.

The target product is a monthly global burned area dataset on a 0.25-degree grid covering 2003-2024. The model uses FireCCI51 as the long-term input product and FireCCIS311 as the reference product during the common period 2019-2024.


Scripts included
================

  4.1-pre-model-auto-rfe.R

  4.2-Harmonised_MRBA60_second_part_parallel3_noLOYO_stats_S3FULL_FINAL.R


General workflow
================

The workflow is designed as a two-step process.

First, 4.1-pre-model-auto-rfe.R selects the best predictors for each biome and month. It applies an autocorrelation filter followed by Recursive Feature Elimination (RFE). The output is a set of CSV files named:

  SelectedPredictors_<biome>_COMMON.csv

Second, 4.2-Harmonised_MRBA60_second_part_parallel3_noLOYO_stats_S3FULL_FINAL.R reads those CSV files and fits the final monthly random forest models by biome. The script then applies the fitted models retrospectively to the full 2003-2024 period and generates harmonised burned area arrays, diagnostics, SHAP outputs and plots.

The LOYO validation code is preserved in the second script but is commented out. The final operational run uses the common model trained with all available common-period data.


Model name and main directories
===============================

Both scripts use:

  Modelo = "B1-MRBA60-2003-2024"

Main input directory:

  /mnt/disco6tb/MRBA60/data/A3_ADJ

Main output directory:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024

Main output subdirectories:

  csv
  plot
  plot_rle
  plot_scatter
  RData
  logs_selected_predictors
  logs_harmonisation


Input datasets
==============

The model uses the following main input datasets:

Burned area and fire-detection variables:

  FireCCI51_2003_2024_0.25degree.RData
  FireCCIS311_2019_2024_0.25degree.RData
  MODIS-AFcount_conf30_angle30-200301-202412-025.RData
  MODIS-FRPsum_conf30_angle30-200301-202412-025.RData
  MODIS-FRPmedian_conf30_angle30-200301-202412-025.RData

Auxiliary predictors:

  ERA5-TEMP-MEAN-2003-2024-MONTLY-025_ADJ.nc
  ERA5-TOT-PREC-2003-2024-MONTLY-025_ADJ.nc
  ERA5-WIND-SPEED-2003-2024-MONTLY-025_ADJ.nc
  ERA5-TOT-CLOUD-2003-2024-MONTLY-025_ADJ.nc
  NDVI-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  VPD-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  SMs-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN.RData

Spatial inputs:

  longitude.RData
  latitude.RData
  continental-biomes_dinerstein_V10.shp
  FireMask_AF3030F.RData

The scripts assume that all predictor products have already been remapped to the common 0.25-degree grid and that missing values over terrestrial fire pixels have already been handled in the previous processing block.


Temporal coverage
=================

The full period is:

  2003-01 to 2024-12

The common calibration period is:

  2019-01 to 2024-12

FireCCI51 covers the full period 2003-2024.

FireCCIS311 covers the common period 2019-2024 and is expanded internally to a full 2003-2024 array, with NA values outside the common period.


Predictor variables
===================

The candidate predictors used in the model are:

  f5                 FireCCI51 burned area
  count_ActiveFire   MODIS active fire count
  prec               ERA5 precipitation
  temp               ERA5 temperature
  FRPsum             MODIS FRP sum
  FRPmedian          MODIS FRP median
  NDVI               MODIS NDVI
  FWI                Number of days with FWI above monthly P95
  wind               ERA5 wind speed
  lat                Latitude
  lon                Longitude
  cloud              ERA5 total cloud cover
  vpd                Vapor pressure deficit
  soil               Surface soil moisture

The response variable is:

  f3                 FireCCIS311 burned area

GPP is explicitly excluded and is not loaded or used in the final model.


4.1 - Predictor selection
=========================

Script:

  4.1-pre-model-auto-rfe.R

Purpose:

  Select the predictor subset to be used later by the final harmonisation model.

Main operations:

  1. Load all global input arrays.
  2. Load the biome shapefile.
  3. Load FireMask_AF3030F and restrict fire-related predictors to the valid fire domain.
  4. Process each biome independently.
  5. For each biome and each month, build a training dataset using the common period.
  6. Use a three-month moving window around each target month.
  7. Remove predictors with zero variance.
  8. Apply a Spearman-correlation clustering filter to reduce collinearity.
  9. Apply Recursive Feature Elimination using random forest.
  10. Save the selected predictors by biome and month.

The three-month window means that, for a target month m, the training data include the previous month, the target month and the following month, constrained to the available common-period dates.

The script runs biomes in parallel using:

  N_CORES <- min(3, max(1, parallel::detectCores() - 1))

The number of processed biomes is controlled by:

  idx_biomas_pendientes

Set this object to NULL to process all biomes, or to a vector of biome indices to process only selected biomes.

Main outputs:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/csv/SelectedPredictors_<biome>_COMMON.csv

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/logs_selected_predictors/LOG_SelectedPredictors_<biome>.txt

The CSV files contain, by biome and month:

  initial predictors
  predictors retained after the autocorrelation filter
  formula after autocorrelation filtering
  predictors retained after RFE
  formula after RFE
  number of training rows
  number of valid years with at least 30 samples


4.2 - Final harmonisation model
===============================

Script:

  4.2-Harmonised_MRBA60_second_part_parallel3_noLOYO_stats_S3FULL_FINAL.R

Purpose:

  Fit the final harmonisation models using the predictor subsets generated in step 4.1 and apply them to the full 2003-2024 period.

Main operations:

  1. Load the selected-predictor CSV files.
  2. Load all global burned area, active fire and auxiliary predictor arrays.
  3. Load the biome shapefile and FireMask_AF3030F.
  4. Process each biome independently in parallel.
  5. For each biome and month, read the selected predictors from the corresponding CSV file.
  6. Fit a random forest model using all common-period training data.
  7. Estimate SHAP values by biome and month.
  8. Apply the model to the common period and to the full 2003-2024 period.
  9. Apply physical constraints to avoid negative values and values above the grid-cell area.
  10. Save local biome-level outputs.
  11. Reconstruct global mosaics sequentially after parallel processing.
  12. Export RData arrays, CSV diagnostics and plots.

The script uses:

  N_CORES <- min(4, max(1, parallel::detectCores() - 1))
  SEED <- 123
  N_SHAP <- 100
  SAVE_TIMESERIES_PLOTS <- TRUE
  SAVE_MODELS <- FALSE

The number of processed biomes is controlled by:

  idx_biomas_pendientes

Set this object to NULL to process all biomes, or to a vector of biome indices to process selected biomes.


Selected-predictor CSV location
===============================

The final model reads selected-predictor CSV files from:

  selected_predictors_dir2 = "/mnt/disco6tb/MRBA60/results/B1-MRBA60-2003-2024/csv-1/"

This path should be checked before running the script. If the CSV files generated by 4.1 are stored in a different folder, this variable must be updated accordingly.

Expected CSV naming pattern:

  SelectedPredictors_<biome>_COMMON.csv


Random forest configuration
===========================

The final model uses random forest regression with:

  ntree = 250
  mtry = min(number_of_predictors, max(1, round(number_of_predictors / 3)))

A separate model is fitted for each biome and month.

The response is FireCCIS311 burned area during the common period.

The predictors are those selected by the previous RFE step for the corresponding biome and month.


Physical constraints
====================

The prediction step applies two physical constraints:

  1. Predictions lower than zero are set to zero.
  2. Predictions larger than the estimated area of the 0.25-degree grid cell are capped to the cell area.

The grid-cell area is estimated using latitude-dependent cell area:

  area = base_cell_area * cos(latitude)

This prevents physically impossible burned area values at the pixel-month scale.


SHAP outputs
============

The final model estimates SHAP values for interpretability using fastshap.

SHAP values are computed by biome and month using:

  N_SHAP = 100

The script stores SHAP-related outputs and uses predictor labels for clearer interpretation in plots and summaries.


Main outputs from 4.2
=====================

The exact set of output files depends on the selected options and the number of processed biomes. The main output categories are:

RData outputs:

  Global harmonised burned area arrays.
  Biome-level local arrays.
  Model outputs and diagnostic objects, depending on SAVE_MODELS.
  SHAP outputs.
  Final merged mosaics.

CSV outputs:

  Predictor tables actually used by biome and month.
  Monthly model evaluation logs.
  Maximum training values by biome and month.
  Time-series and performance statistics.

Plot outputs:

  Common-period time series comparing FireCCIS311, FireCCI51 and harmonised burned area.
  Full-period time series.
  Scatter and diagnostic plots.
  SHAP and predictor-importance visualisations.

Log outputs:

  logs_harmonisation/LOG_Harmonisation_<biome>.txt


Execution order
===============

The recommended execution order is:

  1. Run 4.1-pre-model-auto-rfe.R
  2. Check the generated SelectedPredictors_<biome>_COMMON.csv files.
  3. Confirm that selected_predictors_dir2 in 4.2 points to the correct CSV folder.
  4. Run 4.2-Harmonised_MRBA60_second_part_parallel3_noLOYO_stats_S3FULL_FINAL.R
  5. Check the logs_harmonisation folder.
  6. Check the global RData outputs and diagnostic plots.

The final harmonisation script should not be run before the selected-predictor CSV files exist.


Required previous processing blocks
===================================

Before running this modelling block, the following previous steps must have been completed:

  1. Burned area RData generation for FireCCI51 and FireCCIS311.
  2. MODIS active fire and FRP preprocessing.
  3. ERA5, NDVI, VPD, SMs and FWI preprocessing to the 0.25-degree grid.
  4. FireMask_AF3030F generation.
  5. Fire-adjusted land/sea mask generation.
  6. KNN filling of auxiliary predictors over valid terrestrial fire pixels.


Software requirements
=====================

The scripts require R and the following main packages:

  dplyr
  lubridate
  sf
  ncdf4
  randomForest
  caret
  parallel
  fastshap
  tidyr
  tibble
  ggplot2
  RColorBrewer
  scales
  grid

The 4.1 script requires caret for Recursive Feature Elimination.

The 4.2 script requires fastshap for SHAP-based model interpretation.


Parallel execution notes
========================

The scripts explicitly limit internal BLAS/OpenMP threading to avoid over-parallelisation:

  OMP_NUM_THREADS = 1
  OPENBLAS_NUM_THREADS = 1
  MKL_NUM_THREADS = 1
  VECLIB_MAXIMUM_THREADS = 1
  NUMEXPR_NUM_THREADS = 1

Parallelisation is applied across biomes, not within individual model fits.

For HPC execution, assign enough memory because each worker loads or receives large 3D arrays. If memory becomes limiting, reduce N_CORES or process fewer biomes at once using idx_biomas_pendientes.


Important checks before running
===============================

Before running 4.1:

  Confirm that all input predictor arrays exist in /mnt/disco6tb/MRBA60/data/A3_ADJ.
  Confirm that FireMask_AF3030F.RData exists.
  Confirm that the biome shapefile exists and contains the field cont_bm.
  Set idx_biomas_pendientes as needed.

Before running 4.2:

  Confirm that the selected-predictor CSV files exist.
  Confirm that selected_predictors_dir2 points to the correct folder.
  Confirm that the selected-predictor CSV files contain the Month column.
  Confirm that predictor names in the CSV files match the internal predictor names.
  Set idx_biomas_pendientes as needed.
  Review SAVE_TIMESERIES_PLOTS and SAVE_MODELS.


Notes
=====

The modelling framework is biome-specific and month-specific. This design allows the relationship between FireCCI51, FireCCIS311 and the auxiliary predictors to vary across ecological regions and seasons.

The calibration target is FireCCIS311 during 2019-2024. The trained relationship is then applied retrospectively to FireCCI51 and auxiliary predictors for the full period 2003-2024.

The final model does not recalculate the RFE step. It depends on the CSV files produced by 4.1.

LOYO validation code remains in the script for traceability but is commented out in this final operational version.

End of README.
