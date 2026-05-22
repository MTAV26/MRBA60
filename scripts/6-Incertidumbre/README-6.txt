Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es


README - 6-Incertidumbre
========================

This folder contains the scripts used to estimate and assemble uncertainty layers for the MRBA60 burned area product.

The uncertainty workflow combines several components:

  1. Uncertainty associated with the final harmonised MRBA60 product.
  2. Uncertainty associated with FireCCI51.
  3. Masks identifying biome-month combinations or extreme-event cells where the final product should inherit FireCCI51 uncertainty.
  4. Standard error layers from FireCCIS311/Sentinel-3 SYN burned area files.

The target period is 2003-2024, with monthly data on the global 0.25-degree grid.


Scripts included
================

  9.0-Incertidumbre.R

  9.1-Incertidumbre-FireCCI51.R

  9.2_maskcaraBiomasNoHarmonised-Incertidumbre.R

  10-estandar_error_S3.R


General workflow
================

The recommended execution order is:

  1. 10-estandar_error_S3.R

     Reads monthly FireCCIS311 NetCDF files and extracts the standard error variable for 2019-2024.

  2. 9.0-Incertidumbre.R

     Estimates uncertainty for the MRBA60 harmonised product using the common-period relationship between the harmonised product and FireCCIS311.

  3. 9.1-Incertidumbre-FireCCI51.R

     Estimates uncertainty for FireCCI51 using the common-period relationship between FireCCI51 and FireCCIS311.

  4. 9.2_maskcaraBiomasNoHarmonised-Incertidumbre.R

     Builds masks for areas or biome-month combinations where uncertainty from FireCCI51 should be used or combined, and generates the final uncertainty layer.

The final objective is to produce an uncertainty field for MRBA60 that reflects both the harmonised model uncertainty and the additional uncertainty associated with cells or periods where the final product is effectively inherited from FireCCI51 or corrected using FireCCI51 thresholds.


Model name and main directories
===============================

The scripts use:

  Modelo = "B1-MRBA60-2003-2024"

Main input directory:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/

Main output directory:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024

Main output subdirectories:

  csv
  plot
  RData
  plot_QM001M


Common input files
==================

The scripts use the following main inputs:

  longitude.RData
  latitude.RData
  FireCCI51_2003_2024_0.25degree.RData
  FireCCIS311_2019_2024_0.25degree.RData
  BA_MRBA60.RData
  BA_harmonised_correctedByF51Tope_B1-MRBA60-2003-2024.RData
  continental-biomes_dinerstein_V10.shp
  F51_maskAboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData

The burned area inputs are converted to km2 using division by 1e6 where required.


Temporal coverage
=================

Full MRBA60 period:

  2003-01 to 2024-12

Common FireCCIS311 calibration/evaluation period:

  2019-01 to 2024-12

FireCCIS311 standard error period read by the script:

  2019-01 to 2024-12

Historical period without FireCCIS311:

  2003-01 to 2018-12


10-estandar_error_S3.R
======================

Purpose
-------

This script reads the monthly FireCCIS311/Sentinel-3 SYN burned area NetCDF files and extracts the standard error variable into a single monthly array.

Main operations
---------------

  1. Define the annual directories for FireCCIS311 PSD_Grid files from 2019 to 2024.

  2. Build the expected monthly file names from January 2019 to December 2024.

  3. Use the special internal version name for July and August 2022:

       fv1.2internal

     Other months use:

       fv1.1

  4. Check missing NetCDF files.

  5. Open the first available NetCDF file to detect the standard-error variable automatically.

  6. Read longitude and latitude.

  7. Read the standard error layer month by month.

  8. Convert fill values and missing values to NA.

  9. Reorient the latitude dimension if needed so that latitude is stored from south to north.

  10. Save the resulting standard error array.

Main output
-----------

The script builds:

  SE_S3

with dimensions:

  [lon, lat, time]

covering 72 monthly layers from 2019-01 to 2024-12.

The output is saved under:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/

Notes
-----

The script detects the standard error variable using variable names matching patterns such as standard error, std err, se or rmse. If multiple candidates exist, the first detected candidate is used. This should be checked in the console output.


9.0-Incertidumbre.R
===================

Purpose
-------

This script estimates uncertainty for the harmonised MRBA60 product.

Main operations
---------------

  1. Load FireCCIS311 for 2019-2024.

  2. Load FireCCI51 for 2003-2024.

  3. Load the final MRBA60 product from:

       BA_MRBA60.RData

  4. Build a full-period FireCCIS311 array by inserting 2019-2024 into the 2003-2024 temporal grid.

  5. Build a common information mask where any of the products has burned area information.

  6. Use the final MRBA60 burned area product as the predictor for uncertainty estimation.

  7. Loop over biomes.

  8. For each biome, compare MRBA60 and FireCCIS311 during the common period.

  9. Bin predicted burned area values into quantile bins.

  10. Estimate absolute RMSE and relative RMSE as a function of predicted burned area.

  11. Fit a polynomial relationship between predicted burned area and RMSE.

  12. Apply the fitted uncertainty relationship to the full 2003-2024 period.

  13. Merge biome-level uncertainty estimates into global arrays.

Main outputs
------------

The script generates global uncertainty arrays, including:

  BA_Incertidumbre_HARMONISED_abs.RData

and associated relative uncertainty outputs, depending on the save commands active in the script.

The core output object is generally:

  global_rmse_abs_full

representing absolute uncertainty for the harmonised product.


9.1-Incertidumbre-FireCCI51.R
=============================

Purpose
-------

This script estimates uncertainty for FireCCI51, using FireCCIS311 as the common-period reference.

Main operations
---------------

  1. Load FireCCIS311 for 2019-2024.

  2. Load FireCCI51 for 2003-2024.

  3. Load MRBA60 for consistency and to build the common information domain.

  4. Build a full-period FireCCIS311 array.

  5. Construct a status mask where FireCCI51, FireCCIS311 or MRBA60 has burned area information.

  6. Loop over biomes.

  7. For each biome, compare FireCCI51 and FireCCIS311 during the common period.

  8. Bin FireCCI51 burned area values into quantile bins.

  9. Estimate absolute and relative RMSE relationships.

  10. Apply the uncertainty model to the full 2003-2024 FireCCI51 period.

  11. Merge biome-level uncertainty estimates into global arrays.

Main outputs
------------

The script generates FireCCI51 uncertainty layers, including:

  BA_Incertidumbre_FireCCI51.RData

The core output object is generally:

  global_rmse_abs_full

representing absolute uncertainty for FireCCI51.


9.2_maskcaraBiomasNoHarmonised-Incertidumbre.R
==============================================

Purpose
-------

This script builds masks for cells or biome-month combinations where the final MRBA60 uncertainty should be adjusted using FireCCI51 uncertainty. It then assembles the final MRBA60 uncertainty field.

Main operations
---------------

  1. Load FireCCIS311, FireCCI51 and the harmonised product corrected by the FireCCI51 threshold step.

  2. Build the common information mask.

  3. Define selected biome-month combinations for which the harmonised product should be treated as not fully harmonised or should inherit FireCCI51 uncertainty.

  4. Build a global TRUE/FALSE mask for those selected biome-month combinations.

  5. Save this mask as:

       mask_truefalse_biomas_meses_B1-MRBA60-2003-2024.RData

  6. Load the FireCCI51 uncertainty layer:

       BA_Incertidumbre_FireCCI51.RData

  7. Load the harmonised-product uncertainty layer:

       BA_Incertidumbre_HARMONISED_abs.RData

  8. Load the FireCCI51 threshold-exceedance mask:

       F51_maskAboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData

  9. Combine the threshold-exceedance mask with the biome-month mask.

  10. Create a union mask:

       Mask_union_global

  11. Save the union mask as:

       MASK_NOHARMONISED.RData

  12. Replace or augment harmonised uncertainty with FireCCI51 uncertainty where the union mask is TRUE.

  13. Save the final MRBA60 uncertainty layer.

Selected biome-month mask
-------------------------

The script explicitly defines biome-month combinations where uncertainty treatment differs. The selected entries include:

  Australia-Tropical Broadleaf Forests:
    February, March, April

  Eurasia-Tundra:
    January, February, December

  Europe-Boreal Forests/Taiga:
    January, February, December

  North America-Boreal Forests/Taiga:
    February

  North America-Mediterranean Forests, Woodlands & Scrub:
    January

  North America-Tundra:
    January, February, March, April, October, November, December

  South America-Temperate Broadleaf & Mixed Forests:
    July, August, September

Main outputs
------------

  mask_truefalse_biomas_meses_B1-MRBA60-2003-2024.RData

  MASK_NOHARMONISED.RData

  Final uncertainty array for MRBA60, commonly saved as:

    BA_Incertidumbre_FireCCI60.RData

or an equivalent final uncertainty RData file depending on the active save command.


Uncertainty method
==================

The uncertainty scripts estimate uncertainty empirically from product differences during the FireCCIS311 common period.

For each biome:

  1. The product to be evaluated is compared with FireCCIS311.

  2. Product burned area values are grouped into quantile bins.

  3. For each bin, the RMSE between the product and FireCCIS311 is computed.

  4. A polynomial relationship is fitted between predicted/product burned area and RMSE.

  5. The fitted relationship is applied to estimate uncertainty over the full 2003-2024 period.

The method produces both absolute and relative uncertainty arrays, although the final MRBA60 workflow mainly uses the absolute uncertainty layer.


Physical and spatial assumptions
================================

The scripts assume:

  The spatial grid is global 0.25 degrees.

  Arrays are stored as [lon, lat, time].

  Burned area values are in km2 after conversion.

  The biome shapefile uses the field cont_bm.

  FireCCIS311 is the reference product for the common period.

  The full monthly sequence contains 264 layers from 2003-01 to 2024-12.

The scripts also compute a latitude-dependent grid-cell area:

  area = (110.57 * 0.25) * (111.32 * 0.25) * cos(latitude)

This area matrix is used for diagnostics and for consistency with previous MRBA60 processing blocks.


Software requirements
=====================

The scripts require R and the following packages:

  sf
  dplyr
  lubridate
  terra
  raster
  ggplot2
  ggtext
  ncdf4
  sp
  fields
  maps
  RColorBrewer
  rworldmap
  graticule
  rnaturalearth
  rnaturalearthdata
  viridis
  caret
  gplots
  dendextend
  corrplot
  randomForest
  tidyr
  tibble
  cowplot
  fastshap
  qmap
  scales

The standard-error script mainly requires:

  ncdf4
  lubridate


Recommended checks after execution
==================================

After running 10-estandar_error_S3.R:

  Check that SE_S3 has dimensions [lon, lat, 72].
  Check that latitude has been reoriented correctly if needed.
  Plot several monthly standard-error layers.

After running 9.0-Incertidumbre.R:

  Check that BA_Incertidumbre_HARMONISED_abs.RData exists.
  Confirm that global_rmse_abs_full has dimensions [1440, 720, 264].
  Check that uncertainty values are non-negative.

After running 9.1-Incertidumbre-FireCCI51.R:

  Check that BA_Incertidumbre_FireCCI51.RData exists.
  Compare FireCCI51 uncertainty patterns against the harmonised uncertainty.

After running 9.2_maskcaraBiomasNoHarmonised-Incertidumbre.R:

  Check the proportion of TRUE values in mask_truefalse.
  Check the union mask Mask_union_global.
  Confirm that FireCCI51 uncertainty is used in the intended cells/months.
  Plot the final uncertainty layer for several months.


Important notes
===============

The uncertainty layer is not a formal pixel-level probabilistic confidence interval. It is an empirical uncertainty estimate derived from product disagreement during the common FireCCIS311 period and extrapolated to the full MRBA60 period.

Cells and months where MRBA60 is restored from or strongly influenced by FireCCI51 should not use only the harmonised-model uncertainty. The mask logic in 9.2 explicitly addresses those cases.

The standard-error layer from FireCCIS311 provides an additional uncertainty input for the Sentinel-3 period, but the empirical RMSE-based uncertainty remains the main basis for the full 2003-2024 MRBA60 uncertainty field.

End of README.
