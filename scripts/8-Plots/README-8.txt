Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es


README - 8-Plots paper
======================

This folder contains the R scripts used to generate the main and supplementary figures, tables and diagnostic visualisations for the MRBA60 paper.

The scripts use the final MRBA60 product, external burned area reference products and validation datasets to produce publication-ready plots. Most outputs are written under:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/

The scripts assume that all previous processing blocks have already been completed, especially:

  1. Final MRBA60 burned area RData product.
  2. Final MRBA60 uncertainty RData product.
  3. Final monthly NetCDF files.
  4. FireCCI51 and FireCCIS311 input products.
  5. External validation products such as ONFIRE, MapBiomas, Landsat/S2 validation layers, MCD64A1 and GFED5.


Scripts included
================

Main figures:

  Figure2_R1.R

  Figure4_0_ONFIRE.R
  Figure4_1BRASIL.R
  Figure4_2_FL10.R

  Figure5_Africa_R1.R
  Figure5_Brasil_R1.R

  Figure6_R1.R

  Figure7_R1(1).R
  Figure7.1_R1.R

  Figure8_R1(1).R

Supplementary figures and tables:

  Figure_S2_R1.R
  Figure_S3_R1.R
  Figure_S4_R1.R
  Figure_S8_R1.R
  Figure_S11a_R1.R
  Figure_S11b_R1.R
  Figure_S14_R1.R
  Figure_S15_R1.R
  Table2_and_Fig_S7_R1.R


General input data
==================

Most scripts use one or more of the following inputs:

MRBA60 final product:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_BA_m2_monthly_2003_2024.RData

Main MRBA60 model folder:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/

Input adjusted data folder:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/

Common core inputs:

  longitude.RData
  latitude.RData
  FireCCI51_2003_2024_0.25degree.RData
  FireCCIS311_2019_2024_0.25degree.RData

Final NetCDF product folder:

  /mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/

External comparison products:

  MCD64A1 monthly burned area
  GFED5 burned area
  ONFIRE regional validation data
  MapBiomas Brazil burned area
  Sentinel-2/Landsat validation layers
  GFED/Giglio region masks


Output root
===========

The paper figures are mainly written to:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/

Several scripts create figure-specific subfolders, for example:

  Figure4_R1
  Figure5_R1
  Figure6_R1
  Figure7_R1
  Figure8_R1
  FigureS11_R1_Robinson
  FigureS11_R1_unc_Robinson
  FigureS14_R1
  FigureS15_R1

Some older or supplementary scripts may also write to:

  /mnt/disco6tb/MRBA60/results/D1-Plots/

Check each script before execution if a consistent output root is required.


Figure2_R1.R
============

Purpose
-------

Generates a global comparison figure based on MRBA60, FireCCI51 and FireCCIS311 time series and/or spatial summaries.

Main inputs
-----------

  BA_MRBA60.RData
  FireCCIS311_2019_2024_0.25degree.RData
  FireCCI51_2003_2024_0.25degree.RData
  longitude.RData
  latitude.RData

Main output folder
------------------

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/

Notes
-----

The script reads MRBA60, FireCCIS311 and FireCCI51, converts burned area to km2 where needed, expands FireCCIS311 into the full 2003-2024 temporal grid, and prepares the data used for Figure 2.


Figure4 scripts
===============

The Figure 4 scripts compare MRBA60 against regional or independent validation datasets.

Figure4_0_ONFIRE.R
------------------

Purpose:

  Compares MRBA60 with ONFIRE regional burned area datasets.

Main inputs:

  /mnt/disco6tb/ONFIRE/*.nc
  /mnt/disco6tb/ONFIRE/domain/domain.shp
  MRBA60 final burned area product

Temporal handling:

  The script reads ONFIRE data from 2003 to 2021, with region-specific trimming:

    Europe: 2003-2015
    Canada-NBAC: 2003-2020
    Other regions: 2003-2021

Processing summary:

  ONFIRE burned area is converted from m2 to km2. Regional time series and comparison statistics are generated for the validation figure.

Figure4_1BRASIL.R
-----------------

Purpose:

  Compares MRBA60, FireCCI51, MCD64A1 and GFED5 against MapBiomas Brazil.

Main inputs:

  MRBA60_BA_m2_monthly_2003_2024.RData
  FireCCI51_2003_2024_0.25degree.RData
  MapBiomas Brazil annual burned area layers
  MCD64A1
  GFED5

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure4_R1/

Figure4_2_FL10.R
----------------

Purpose:

  Produces another Figure 4 validation/comparison panel, using the same product group:

    MRBA60
    FireCCI51
    MCD64A1
    GFED5

Processing summary:

  Monthly products are aggregated to annual burned area, then compared over the corresponding validation domain.

Notes
-----

The Figure 4 scripts use a common colour map:

  MRBA60: blue
  FireCCI51: orange
  MCD64A1: brown
  GFED5: grey/black


Figure5 scripts
===============

Figure5_Africa_R1.R
-------------------

Purpose:

  Generates a regional validation figure for Africa using independent reference data, MRBA60 and FireCCI51.

Main processing:

  The script loads validation burned area data, reads MRBA60 and FireCCI51, converts units to km2, aligns products to the validation grid and produces regional comparison panels.

Figure5_Brasil_R1.R
-------------------

Purpose:

  Generates a Brazil-focused validation figure using MapBiomas.

Main inputs:

  MapBiomas annual burned area layers
  MRBA60_BA_m2_monthly_2003_2024.RData
  FireCCI51_2003_2024_0.25degree.RData

Key settings:

  year_target = 2005
  buffer_dist_m = 25000
  CRS for Brazil metric operations = EPSG:5880

Processing summary:

  The script compares annual burned area in Brazil while avoiding boundary effects using a buffer-based approach.


Figure6_R1.R
============

Purpose
-------

Generates global maps of annual mean burned area and product differences.

Products included
-----------------

  MRBA60
  FireCCI51
  MCD64A1
  GFED5

Main outputs
------------

  Publication maps of:

    Mean annual burned area.
    Difference maps between MRBA60 and comparison products.

Important plotting behaviour
----------------------------

  Zero and NA values are left without colour.

  Mean burned area maps use discrete burned-area classes.

  Difference maps use discrete difference classes.

  The script is configured to write to:

    /mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure6_R1/


Figure7 scripts
===============

Figure7_R1(1).R
---------------

Purpose:

  Generates the main global burned area time-series figure.

Products:

  MRBA60
  FireCCI51
  MCD64A1
  GFED5

Output:

  One multi-panel figure with:

    a) Annual burned area
    b) DJF
    c) MAM
    d) JJA
    e) SON

Design:

  The figure uses one column and five rows. Lines show burned area by product. Grey bars show the difference MRBA60 - FireCCI51. The legend is placed below the seasonal panels.

Main input:

  MRBA60_BA_m2_monthly_2003_2024.RData

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure7_R1/

Figure7.1_R1.R
---------------

Purpose:

  Generates trend and summary statistics for the Figure 7 text.

Main outputs:

  Annual_BA_stats_MRBA60_Figure7.csv
  Seasonal_BA_stats_MRBA60_Figure7.csv
  Annual_and_Seasonal_BA_stats_MRBA60_Figure7.csv
  Figure7_text_numbers_summary.txt
  Figure7_BA_stats_MRBA60.xlsx

Methods:

  The script uses modified Mann-Kendall and Sen slope functions loaded from external files:

    mmkh.R
    sen.R

The generated tables provide annual and seasonal trend statistics used in the manuscript text.


Figure8_R1(1).R
===============

Purpose
-------

Generates a regional analysis using GFED/Giglio regions.

Products
--------

  MRBA60
  GFED5
  MCD64A1

Common period
-------------

  2003-2024

Main analyses
-------------

  Regional monthly climatology.

  Spatial Spearman correlation by region and month:

    MRBA60 vs GFED5
    MRBA60 vs MCD64A1

  Zeros are included in the correlation calculations.

Main outputs
------------

  PDF figure.

  CSV, RDS and Excel tables.

Main output folders:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure8_R1
  /mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure8_R1/tablas

Input region file:

  GFED5_Beta_monthly_2002.nc


Supplementary figures
=====================

Figure_S2_R1.R
--------------

Purpose:

  Summarises the selected predictors from the SelectedPredictors_*.csv files.

Main input:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/csv/SelectedPredictors_*.csv

Processing summary:

  Reads the predictor-selection CSV files, cleans biome names and summarises predictor occurrence after the autocorrelation filtering and after RFE.

Figure_S3_R1.R
--------------

Purpose:

  Generates a predictor frequency plot based on the selected-predictor CSV files.

Main input:

  SelectedPredictors_*.csv

Key setting:

  N_total = 579

This is used as the denominator for the total number of biome-month combinations considered.

Figure_S4_R1.R
--------------

Purpose:

  Generates a biome-month predictor selection summary from the RFE selected predictors.

Main input:

  SelectedPredictors_*.csv

Processing summary:

  Reads RFE predictors from each biome CSV, renames predictors using publication-friendly labels and prepares a plot of predictor selection patterns.

Figure_S8_R1.R
--------------

Purpose:

  Generates a Brazil regional zoom figure.

Products:

  FireCCI51
  MRBA60
  MapBiomas
  Difference maps

Design:

  Individual panels a-d plus a 1 x 4 combined layout.

  Includes an inset overview of Brazil with a zoom rectangle.

  Uses a fixed zoom window approximately between 15 degrees S and 5 degrees S.

  Uses a continuous inferno palette.

Main output folder:

  /mnt/disco6tb/MRBA60/results/D1-Plots/Figure_S8_R1

Figure_S11a_R1.R
----------------

Purpose:

  Generates monthly burned area maps from the final MRBA60 NetCDF product.

Input:

  /mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/

Output:

  Monthly PDF and JPEG maps.

Projection:

  Robinson projection.

Domain:

  Antarctica excluded using lat > -60.

Variable:

  burned_area

Units in plotting:

  km2

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS11_R1_Robinson

Figure_S11b_R1.R
----------------

Purpose:

  Generates monthly uncertainty maps from the final MRBA60 NetCDF product.

Input:

  /mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/

Output:

  Monthly PDF and JPEG uncertainty maps.

Projection:

  Robinson projection.

Domain:

  Antarctica excluded using lat > -60.

Variable:

  uncertainty

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS11_R1_unc_Robinson

Figure_S14_R1.R
---------------

Purpose:

  Generates a map of the mean annual difference between FireCCIS311 and FireCCI51 for 2019-2024.

Difference:

  FireCCIS311 - FireCCI51

Units:

  km2 yr-1 mean over 2019-2024

Projection:

  Robinson projection.

Domain:

  Antarctica excluded using lat > -60.

Important plotting behaviour:

  Differences exactly equal to zero are left without colour.

  The central class from -5 to 5 is shown in light yellow.

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS14_R1/

Figure_S15_R1.R
---------------

Purpose:

  Performs a sensitivity analysis comparing:

    MRBA60H, trained with six common years and covering 2003-2024.
    FireCCI60, trained with four common years and covering 2003-2022.
    FireCCI51, original ESA product.

Evaluation period:

  2003-2018

Main output folder:

  /mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS15_R1/

The script writes outputs under subfolders:

  csv
  plot
  RData


Table2_and_Fig_S7_R1.R
======================

Purpose
-------

Generates Table 2 and Supplementary Figure S7.

Main input
----------

  MRBA60_BA_m2_monthly_2003_2024.RData

Processing summary
------------------

The script loads MRBA60, splits the historical harmonised period and the operational/Sentinel-3 period, computes product summaries and prepares the table and figure outputs.

Temporal structure:

  2003-2018:
    Historical MRBA60 harmonised period.

  2019-2024:
    MRBA60 operational/Sentinel-3 period.

The script is configured so that both historical and operational phases are labelled consistently as MRBA60.


Software requirements
=====================

Across the full plotting folder, the scripts use the following R packages:

  ncdf4
  lubridate
  fields
  terra
  raster
  maps
  ggplot2
  gridExtra
  Metrics
  dplyr
  tidyr
  purrr
  tibble
  sf
  rnaturalearth
  rnaturalearthdata
  scales
  patchwork
  viridis
  viridisLite
  RColorBrewer
  cowplot
  ggtext
  readr
  openxlsx
  forcats
  stringr
  abind

Not every package is used by every script. Several scripts load broad package sets to preserve compatibility with earlier plotting workflows.


Recommended execution order
===========================

A practical execution order is:

  1. Run Figure2_R1.R after the final MRBA60 product is available.

  2. Run validation figures:

       Figure4_0_ONFIRE.R
       Figure4_1BRASIL.R
       Figure4_2_FL10.R
       Figure5_Africa_R1.R
       Figure5_Brasil_R1.R

  3. Run global product comparison maps:

       Figure6_R1.R

  4. Run global and seasonal time-series figures and statistics:

       Figure7_R1(1).R
       Figure7.1_R1.R

  5. Run regional/GFED analysis:

       Figure8_R1(1).R

  6. Run supplementary predictor-selection figures:

       Figure_S2_R1.R
       Figure_S3_R1.R
       Figure_S4_R1.R

  7. Run supplementary maps and sensitivity tests:

       Figure_S8_R1.R
       Figure_S11a_R1.R
       Figure_S11b_R1.R
       Figure_S14_R1.R
       Figure_S15_R1.R
       Table2_and_Fig_S7_R1.R

The scripts are mostly independent once all input products exist. Some scripts require objects generated or loaded in previous plotting scripts, especially Figure_S8_R1.R, which explicitly expects several objects already in memory. If running scripts independently, check the required-object section at the top of each script.


Important checks before running
===============================

Before running the plotting scripts, check:

  The final MRBA60 RData product exists.

  The final NetCDF folder exists if running monthly map scripts.

  External validation datasets are available at the paths hard-coded in the scripts.

  Output directories point to the intended paper-figure folder.

  Unit conversions are correct:
    MRBA60 NetCDF/RData products are often loaded in m2 and converted to km2 for plotting.
    FireCCI51 and FireCCIS311 raw RData products are divided by 1e6 where needed.

  File paths using older folders such as /mnt/disco6tb/MRBA60/results should be reviewed if all outputs should go to /mnt/disco6tb/MRBA60-2/results.


Notes
=====

These scripts are figure-production scripts, not core processing scripts. They are designed to reproduce the paper figures and associated tables using the final MRBA60 outputs.

Several scripts include hard-coded paths to external datasets and validation products. Those paths must exist on the execution machine.

Several scripts use publication-specific styling, colour palettes, panel labels, axis limits and file naming. Avoid changing them unless the manuscript figure layout is being revised.

Some scripts still refer to historical product names such as FireCCI60 or MRBA60H in comments or variable names. In the final paper-facing labels, the scripts generally aim to use MRBA60 consistently.

End of README.
