Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es


README - 3-apply-mask
=====================

This folder contains the scripts used to apply the final fire/land mask to the auxiliary monthly predictor datasets used in the MRBA60 workflow. The scripts fill missing predictor values only in valid terrestrial fire pixels, using a K-nearest-neighbour (KNN) approach constrained by the adjusted land/sea mask and the monthly FireMask_AF3030F mask.

The target period is 2003-01 to 2024-12, with all products expected on a global 0.25-degree grid.

The scripts included in this step are:

  D2-NDVI-MASK_025.R
  E2-VPD-MASK_025.R
  F2-SMs-MASK_025.R
  G2-FWI-MASK_025.R


Purpose of this step
====================

The aim of this processing block is to ensure that the auxiliary predictors used later in the MRBA60 modelling framework are spatially and temporally consistent with the fire mask.

For each monthly predictor cube, the scripts identify pixels that meet all of the following conditions:

  1. FireMask_AF3030F == 1 for that pixel and month.
  2. The predictor value is missing.
  3. The pixel is classified as land in the fire-adjusted land/sea mask.

Those missing values are then filled using nearby valid terrestrial pixels from the same monthly layer. Ocean pixels are not filled.


Common input files
==================

All scripts use the following common inputs:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData

  /mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/FireMask_AF3030F.RData

  /mnt/disco6tb/MRBA60/data/A3_ADJ/land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc

The fire mask must have already been generated in the previous processing block. The adjusted land/sea mask must also already exist. This adjusted mask reclassifies as land any ocean cell that had fire detections at least once during 2003-2024.


KNN configuration
=================

The KNN configuration is consistent across the scripts:

  k_nn = 8
  max_dist_cells = 8
  grid_res_deg = 0.25, where explicitly defined

The maximum distance corresponds approximately to 2 degrees on the 0.25-degree grid.

For NDVI, the script fills missing values only when valid neighbours are found within the maximum distance. For VPD, SMs and FWI95d, the scripts also include a fallback strategy: if no valid terrestrial neighbour is available within the maximum distance, the nearest available valid terrestrial neighbours are used and the distance is recorded in the log.

For FWI95d, the filled values are additionally constrained to the physically meaningful monthly range from 0 to the number of days in the month. The script is configured to round the filled values to integer day counts.


Script descriptions
===================

D2-NDVI-MASK_025.R
------------------

This script processes monthly MOD13C2 NDVI for 2003-2024.

Main input:

  /mnt/disco6tb/MRBA60/data/A2_TEMP/MOD13C2_NDVI_2003_2024_025.nc

Main output:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/NDVI-2003_2024-MONTHLY-025-mask-landsea-KNN.RData

Processing summary:

  The script reads the NDVI NetCDF file, detects the NDVI variable automatically, checks that the spatial and temporal dimensions match FireMask_AF3030F, and fills missing NDVI values only where fire was detected and the adjusted land/sea mask indicates land.

  Missing values are filled using the mean of up to eight valid neighbouring land pixels within the maximum distance of eight grid cells. The script also creates a monthly fill log in memory and reports the number of pixels filled, not filled, and remaining after the KNN procedure.


E2-VPD-MASK_025.R
-----------------

This script processes monthly TerraClimate VPD for 2003-2024.

Main input:

  /mnt/disco6tb/MRBA60/data/A2_TEMP/vpd_2003-2024_0.25deg_bil.nc

Main outputs:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/VPD-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-with-original.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-fill_log.csv
  /mnt/disco6tb/MRBA60/data/A3_ADJ/VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-summary.RData

Processing summary:

  The script reads the VPD NetCDF file, identifies the VPD variable automatically, stores a copy of the original cube, and applies KNN filling only to missing terrestrial fire pixels.

  It records whether each filled value was obtained within the maximum search distance or using the fallback nearest-neighbour strategy beyond the initial radius. The script also writes a CSV fill log and an RData summary.


F2-SMs-MASK_025.R
-----------------

This script processes monthly GLEAM surface soil moisture (SMs) for 2003-2024.

Main input:

  /mnt/disco6tb/MRBA60/data/A2_TEMP/SMs_2003_2024_GLEAM_v4.2b_MO_025deg_bil.nc

Main outputs:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/SMs-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/SMs-2003_2024-MONTHLY-025-mask-landsea-KNN-with-original.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/SMs-2003_2024-MONTHLY-025-mask-landsea-KNN-fill_log.csv
  /mnt/disco6tb/MRBA60/data/A3_ADJ/SMs-2003_2024-MONTHLY-025-mask-landsea-KNN-summary.RData

Processing summary:

  The script reads the SMs NetCDF file, detects the soil moisture variable automatically, keeps the original cube, and fills missing terrestrial fire pixels using KNN. As in the VPD script, cases filled beyond the initial maximum search radius are explicitly logged.


G2-FWI-MASK_025.R
-----------------

This script processes monthly FWI95d for 2003-2024. FWI95d represents the number of days per month with FWI above the monthly 95th percentile.

Main input:

  /mnt/disco6tb/MRBA60/data/A2_TEMP/FWI/fwi-era5_count_exceed95_200301-202412_0.25.nc

Main outputs:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN-with-original.RData
  /mnt/disco6tb/MRBA60/data/A3_ADJ/FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN-fill_log.csv
  /mnt/disco6tb/MRBA60/data/A3_ADJ/FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN-summary.RData

Processing summary:

  The script reads the monthly FWI95d NetCDF file, detects the FWI/count variable automatically, stores the original cube, and fills missing terrestrial fire pixels using KNN.

  Because FWI95d is a monthly count of days, filled values are constrained between 0 and the number of days in the corresponding month. The default configuration rounds the filled values to integer counts.


Expected order of execution
===========================

The recommended execution order is:

  1. D2-NDVI-MASK_025.R
  2. E2-VPD-MASK_025.R
  3. F2-SMs-MASK_025.R
  4. G2-FWI-MASK_025.R

The scripts are independent once the common masks and input NetCDF files exist. They can therefore be run separately or submitted as independent jobs, provided that enough memory is available.


Required previous steps
=======================

Before running these scripts, the following files must already exist:

  1. FireMask_AF3030F.RData

     Generated during the mask-building and post-processing block. It is the monthly 2003-2024 mask indicating whether a pixel/month is considered part of the valid fire domain.

  2. land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc

     Generated by the land/sea adjustment script. It is used to prevent KNN filling over ocean pixels while preserving near-coastal fire pixels that were reclassified as land.

  3. Predictor NetCDF files already harmonised to 0.25 degrees:

     MOD13C2_NDVI_2003_2024_025.nc
     vpd_2003-2024_0.25deg_bil.nc
     SMs_2003_2024_GLEAM_v4.2b_MO_025deg_bil.nc
     fwi-era5_count_exceed95_200301-202412_0.25.nc


Software requirements
=====================

The scripts require R and the following packages:

  ncdf4
  RANN
  ggplot2
  viridis
  patchwork
  rnaturalearth
  sf
  grid
  ggtext

The NDVI script only requires ncdf4 and RANN for the core processing, although plotting and diagnostics may use base R graphics.


Main assumptions
================

The scripts assume that:

  The spatial grid is global 0.25 degrees.

  The temporal period is monthly from January 2003 to December 2024.

  Arrays are stored as [lon, lat, time] or are already compatible with the dimensions of FireMask_AF3030F and the adjusted land/sea mask.

  FireMask_AF3030F and the auxiliary predictors have the same spatial and temporal dimensions.

  The adjusted land/sea mask uses:
    0 = ocean
    1 = land

  Missing values are represented as NA after reading the NetCDF variables.


Quality-control outputs
=======================

The VPD, SMs and FWI95d scripts generate CSV and RData logs summarising the KNN filling process. These logs include, depending on the script:

  date
  number of missing terrestrial fire pixels
  number of valid terrestrial donor pixels
  number of filled pixels
  number of pixels filled within the maximum distance
  number of pixels filled beyond the maximum distance
  minimum, mean and maximum distances for fallback cases
  number of remaining missing pixels

These files should be checked after execution to verify whether any predictor still contains missing values over terrestrial fire pixels.


Notes
=====

The KNN filling is performed independently for each month. Donor pixels are taken only from the same monthly layer, so there is no temporal interpolation.

The scripts do not fill ocean pixels. Pixels with fire detections that were originally classified as ocean should already have been handled by the fire-adjusted land/sea mask.

For VPD, SMs and FWI95d, values filled beyond the initial search radius are not silently accepted. They are explicitly recorded in the output logs, which allows later quality control.

For FWI95d, the filled values are treated as monthly counts rather than continuous meteorological values. This is why the script constrains values to the valid range of each month and rounds them by default.

End of README.
