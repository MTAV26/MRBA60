README.txt
Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es

Last updated: 2026-05-22

============================================================
1. PURPOSE
============================================================

This folder contains the scripts used to build and post-process the MRBA60 fire-detection mask and probability-based active-fire filtering workflow.

The workflow combines the following sources on a global 0.25-degree grid:

  - FireCCI51 burned area, 2003-2024
  - FireCCIS311 burned area, 2019-2024
  - MODIS MCD14ML active-fire count filtered at confidence >= 30 and scan angle <= 30 degrees
  - Continental biome polygons
  - ATLAS/reference land-sea mask

The main objective is to derive a consistent fire mask for MRBA60 by combining burned-area products and active-fire information, estimating monthly biome-specific probabilities, filtering likely false active-fire detections using ROC/Youden thresholds, and adjusting the land-sea mask where valid fire detections occur in cells originally classified as ocean.

The scripts are designed for the MRBA60-2003-2024-V1 processing chain.

============================================================
2. SCRIPT OVERVIEW
============================================================

A1-MHBA60-PROB-2003-2024-CL30-R30.R
  Builds probability layers for the harmonised burned-area/fire-detection framework.

  Main role:

    - Loads FireCCI51, FireCCIS311 and MODIS active-fire count.
    - Uses FireCCI51 and active-fire count as predictors.
    - Uses FireCCIS311 as the response during the common period 2019-2024.
    - Processes data by biome and month.
    - Produces probability arrays for the common period and for the full 2003-2024 period.

  Main input directory:

    /mnt/disco6tb/MRBA60/data/A3_ADJ

  Required inputs:

    longitude.RData
    latitude.RData
    FireCCI51_2003_2024_0.25degree.RData
    FireCCIS311_2019_2024_0.25degree.RData
    MODIS-AFcount_conf30_angle30-200301-202412-025.RData
    /mnt/disco6tb/MRBA60/data/A1_RAW/MBC/continental-biomes_dinerstein_V10.shp

  Main output directory:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/RData

  Main outputs:

    BA_PROB_MRBA60-2003-2024-V1_<BIOME>_FireHarmonized_Common.RData
    BA_PROB_MRBA60-2003-2024-V1_<BIOME>_FireHarmonized_Full.RData
    BA_MRBA60-2003-2024-V1global_BA_PROB_FireHarmonized_Common.RData
    BA_MRBA60-2003-2024-V1global_BA_PROB_FireHarmonized_Full.RData

  Notes:

    - FireCCIS311 is expanded internally to the full 2003-2024 timeline, with valid values only for 2019-2024.
    - The common calibration period is 2019-2024.
    - Values equal to zero are initially converted to NA and then restored to zero inside valid fire-status cells.
    - The active-fire metric currently used is AFcount.

------------------------------------------------------------

A2-MHBA60-POSTPROCES-ROC-cell.R
  Applies the probability-based postprocessing step using ROC-derived thresholds.

  Main role:

    - Loads FireCCI51, FireCCIS311, active-fire count and the probability output from A1.
    - Builds a combined valid-status mask from active fires, FireCCI51 and FireCCIS311.
    - Computes ROC statistics by biome and month.
    - Uses the Youden index to select probability thresholds.
    - Filters likely false active-fire-only detections.
    - Saves the filtered active-fire mask and ROC threshold tables.

  Main input directory:

    /mnt/disco6tb/MRBA60/data/A3_ADJ

  Additional required input from A1:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/RData/BA_MRBA60-2003-2024-V1global_BA_PROB_FireHarmonized_Full.RData

  Main output directory:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC

  Main outputs:

    RData/MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData
    csv/ROC_Thresholds_ByBioma_Month_MASK.csv
    RData/ROC_Thresholds_ByBioma_Month_MASK.RData

  Notes:

    - The final filtered object is saved as count_ActiveFire_tot_filtered.
    - The filtering rule sets active-fire-only detections to zero when their estimated probability is below the biome-month threshold and there is no independent burned-area detection.
    - For 2003-2018, independent detection is assessed using FireCCI51.
    - For 2019-2024, independent detection is assessed using FireCCI51 or FireCCIS311.

------------------------------------------------------------

A3-MHBA60-RF-PROB-MASK.R
  Builds the final monthly FireMask_AF3030F after ROC filtering and produces annual status maps.

  Main role:

    - Loads FireCCI51 and FireCCIS311 burned area.
    - Loads the filtered active-fire mask produced by A2.
    - Expands FireCCIS311 to the full 2003-2024 timeline.
    - Builds FireMask_AF3030F as a monthly binary mask.
    - Produces annual maps showing agreement and disagreement among Active Fire, FireCCI51 and FireCCIS311.

  Required input from A2:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData

  Main output directory:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC

  Main outputs:

    RData/FireMask_AF3030F.RData
    MAPS_ANNUAL_STATUS/mapa_all_<YEAR>_WGS84_025_FILTER.jpg

  FireMask_AF3030F rule:

    FireMask_AF3030F = 1 if at least one of the following is true for a given cell and month:

      - filtered Active Fire > 0
      - FireCCI51 > 0
      - FireCCIS311 > 0

    FireMask_AF3030F = 0 otherwise.

  Annual status-map classes:

    0 = only Active Fire
    1 = only FireCCIS311
    2 = only FireCCI51
    3 = Active Fire + FireCCIS311
    4 = FireCCI51 + FireCCIS311
    5 = Active Fire + FireCCI51
    6 = all three datasets

------------------------------------------------------------

A4-MHBA60-RF-PROB-plot-NOfilter.R
  Produces annual diagnostic status maps before applying the ROC filtering.

  Main role:

    - Loads raw FireCCI51, FireCCIS311 and active-fire count.
    - Expands FireCCIS311 to 2003-2024.
    - Builds annual presence/absence combinations without using the filtered active-fire mask.
    - Saves annual diagnostic maps for comparison with the filtered outputs from A3.

  Main output directory:

    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/plot_mask_NoFilter

  Main outputs:

    mapa_all_<YEAR>_WGS84_025_NOFILTER.jpg

  Notes:

    - This script is diagnostic only.
    - It should be used to compare the spatial distribution of detections before and after the ROC-based filtering.

------------------------------------------------------------

A5_mask_land_sea_ATLAS.R
  Adjusts the binary land-sea mask using the final fire mask.

  Main role:

    - Loads the ATLAS/reference land-sea mask at 0.25 degrees.
    - Loads FireMask_AF3030F from A3.
    - Identifies cells classified as ocean that have fire at least once during 2003-2024.
    - Reclassifies those ocean cells as land.
    - Saves the adjusted land-sea mask as both RData and NetCDF.

  Required inputs:

    /mnt/disco6tb/MRBA60/data/A2_TEMP/land_sea_mask_025degree_binary_1440x720.nc
    /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/FireMask_AF3030F.RData

  Main outputs:

    /mnt/disco6tb/MRBA60/data/A3_ADJ/land_sea_mask_025degree_binary_1440x720_fire_adjusted.RData
    /mnt/disco6tb/MRBA60/data/A3_ADJ/land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc

  Output convention:

    sftlf = 0 means ocean
    sftlf = 1 means land

  Notes:

    - The adjustment prevents valid fire cells near coastlines or small islands from being excluded later because of the original land-sea classification.
    - The output NetCDF includes summary attributes documenting how many cells were modified.

============================================================
3. RECOMMENDED EXECUTION ORDER
============================================================

The scripts should be run in the following order:

  1. A1-MHBA60-PROB-2003-2024-CL30-R30.R
     Build probability layers using FireCCI51, FireCCIS311 and active-fire count.

  2. A2-MHBA60-POSTPROCES-ROC-cell.R
     Estimate ROC/Youden thresholds and filter low-probability active-fire-only detections.

  3. A3-MHBA60-RF-PROB-MASK.R
     Build the final FireMask_AF3030F and generate annual status maps after filtering.

  4. A4-MHBA60-RF-PROB-plot-NOfilter.R
     Optional diagnostic step. Generate annual status maps before filtering.

  5. A5_mask_land_sea_ATLAS.R
     Adjust the land-sea mask using the final FireMask_AF3030F.

A4 can be run before or after A3, because it does not depend on the filtered mask. It is included as a diagnostic comparison against the filtered maps.

============================================================
4. MAIN INPUTS REQUIRED BEFORE RUNNING THIS WORKFLOW
============================================================

The following files should already exist before running A1-A5:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/FireCCI51_2003_2024_0.25degree.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/FireCCIS311_2019_2024_0.25degree.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/MODIS-AFcount_conf30_angle30-200301-202412-025.RData
  /mnt/disco6tb/MRBA60/data/A1_RAW/MBC/continental-biomes_dinerstein_V10.shp
  /mnt/disco6tb/MRBA60/data/A2_TEMP/land_sea_mask_025degree_binary_1440x720.nc

The active-fire count file should be generated beforehand from the MODIS MCD14ML preprocessing workflow.

============================================================
5. MAIN OUTPUTS OF THE WORKFLOW
============================================================

Probability outputs from A1:

  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/RData/BA_MRBA60-2003-2024-V1global_BA_PROB_FireHarmonized_Common.RData
  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/RData/BA_MRBA60-2003-2024-V1global_BA_PROB_FireHarmonized_Full.RData

ROC and filtered active-fire outputs from A2:

  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/csv/ROC_Thresholds_ByBioma_Month_MASK.csv
  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/ROC_Thresholds_ByBioma_Month_MASK.RData
  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData

Final fire mask from A3:

  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/FireMask_AF3030F.RData

Annual diagnostic maps:

  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/MAPS_ANNUAL_STATUS/mapa_all_<YEAR>_WGS84_025_FILTER.jpg
  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/plot_mask_NoFilter/mapa_all_<YEAR>_WGS84_025_NOFILTER.jpg

Adjusted land-sea mask from A5:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/land_sea_mask_025degree_binary_1440x720_fire_adjusted.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc

============================================================
6. SOFTWARE DEPENDENCIES
============================================================

R packages used across the scripts include:

  caret
  corrplot
  cowplot
  data.table
  dendextend
  dplyr
  fastshap
  fields
  ggplot2
  ggpmisc
  ggpointdensity
  ggtext
  gplots
  grid
  lubridate
  MASS
  maps
  ncdf4
  pheatmap
  pROC
  randomForest
  raster
  RColorBrewer
  rnaturalearth
  rnaturalearthdata
  rworldmap
  scales
  sf
  sp
  stringr
  terra
  tibble
  tidyr
  viridis

Not every script uses every package. The most important packages for this workflow are sf, terra/raster, ncdf4, pROC, dplyr, fields and maps.

============================================================
7. KEY TEMPORAL AND SPATIAL SETTINGS
============================================================

Spatial grid:

  Global regular lon-lat grid at 0.25 degrees.

Expected dimensions:

  longitude: 1440 cells
  latitude: 720 cells
  time: 264 monthly layers for 2003-01 to 2024-12

Main temporal periods:

  Full period: 2003-01 to 2024-12
  Common FireCCI51-FireCCIS311 period: 2019-01 to 2024-12

Active-fire filter represented in filenames:

  conf30_angle30

This indicates confidence >= 30 and scan angle <= 30 degrees.

============================================================
8. QUALITY-CONTROL CHECKS INCLUDED IN THE SCRIPTS
============================================================

The scripts include multiple internal checks:

  - Verification that all required input files exist.
  - Verification that burned-area, active-fire and probability arrays have three dimensions.
  - Verification that the temporal dimension matches the expected number of months.
  - Automatic detection or standardisation of spatial order to [lon, lat, time].
  - Expansion of FireCCIS311 from 2019-2024 to the full 2003-2024 timeline.
  - Verification that FireMask_AF3030F matches the presence of at least one detection source.
  - Verification that, after the land-sea adjustment, no cell with fire occurrence remains classified as ocean.

============================================================
9. IMPORTANT NOTES
============================================================

1. The current model name is:

     MRBA60-2003-2024-V1

2. The current active-fire metric is:

     AFcount

3. The final mask name is:

     FireMask_AF3030F

4. The label AF3030F refers to active-fire data filtered using:

     confidence >= 30
     scan angle <= 30 degrees
     final postprocessing applied

5. The A2 output name still contains the string ONLYFRP:

     MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData

   However, in the current configuration the active-fire metric is AFcount. This name is retained by the script and should be interpreted as the filtered active-fire mask used in the current workflow.

6. The scripts assume that the required FireCCI51, FireCCIS311, active-fire and coordinate files were already created by the previous preprocessing workflow.

7. The A5 script should be run only after FireMask_AF3030F has been generated by A3.

============================================================
10. MINIMAL COMMAND EXAMPLES
============================================================

From an R-enabled environment:

  Rscript A1-MHBA60-PROB-2003-2024-CL30-R30.R
  Rscript A2-MHBA60-POSTPROCES-ROC-cell.R
  Rscript A3-MHBA60-RF-PROB-MASK.R
  Rscript A4-MHBA60-RF-PROB-plot-NOfilter.R
  Rscript A5_mask_land_sea_ATLAS.R

If running on an HPC system, these commands can be included in a Slurm job script, making sure that the required R modules and geospatial libraries are loaded beforehand.

============================================================
11. FINAL PRODUCTS GENERATED BY THIS BLOCK
============================================================

The two most important products generated by this workflow are:

  1. FireMask_AF3030F.RData

     Monthly binary fire mask for 2003-2024, where each cell-month is set to 1 if there is evidence of fire in the filtered active-fire product, FireCCI51 or FireCCIS311.

  2. land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc

     Fire-adjusted binary land-sea mask in which ocean cells with fire occurrence during 2003-2024 are reclassified as land.

These outputs are intended to support the subsequent MRBA60 harmonisation and masking steps.
