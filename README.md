# Harmonised Medium Resolution Burned Area product version 6.0.0 (MRBA60)

**Miguel Ángel Torres-Vázquez**  
miguela.torres@uah.es  

This repository contains the scripts required to reproduce the analyses and results presented in the study **“MRBA60: Long-term consistent global burned area dataset from Sentinel-3 and MODIS products”**.

MRBA60 is a long-term monthly global burned area product at 0.25-degree spatial resolution. It combines the temporal coverage of MODIS-based FireCCI51 with the Sentinel-3-based FireCCIS311/MRBA60 reference period through a biome- and month-specific harmonisation framework.

---

## 1. Repository purpose

This repository documents the complete processing chain used to generate, evaluate and export the MRBA60 product version 6.0.0.

The workflow covers:

1. Preprocessing of burned area, active fire and environmental predictor datasets.
2. Construction of the fire mask and fire-adjusted land/sea mask.
3. Application of the fire/land mask to auxiliary predictors.
4. Biome- and month-specific model training.
5. Post-processing of harmonised burned area.
6. Estimation and assembly of uncertainty layers.
7. Export of final monthly NetCDF files.
8. Generation of paper figures, supplementary figures and tables.

The scripts are designed around a global monthly grid of:

```text
longitude x latitude x time = 1440 x 720 x 264
```

covering:

```text
2003-01 to 2024-12
```

The common FireCCI51–FireCCIS311 calibration period is:

```text
2019-01 to 2024-12
```

---

## 2. Product overview

The final product is exported as monthly NetCDF files with:

```text
burned_area  units: m2
uncertainty  units: m2
```

Final NetCDF file naming convention:

```text
YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc
```

Example:

```text
20190801-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc
```

The final temporal composition is:

```text
2003-2018  MRBA60 harmonised historical period
2019-2024  FireCCIS311 / Sentinel-3 reference period
```

---

## 3. Main workflow structure

The repository is organised into eight major processing blocks.

```text
1-preprocess
2-build-mask
3-apply-mask
4-model
5-post
6-Incertidumbre
7-NetCDF
8-Plots paper
```

Each block has its own detailed README file. This main README provides the global reproducibility scheme and the expected order of execution.

---

## 4. Reproducibility scheme

The complete workflow should be executed in the following order.

### Step 1. Preprocess input datasets

Folder:

```text
1-preprocess
```

Purpose:

Prepare all burned area, active fire and auxiliary predictor datasets on the common 0.25-degree grid.

Main operations:

- Build FireCCI51 and FireCCIS311 burned area RData arrays.
- Rasterise MODIS MCD14ML active-fire detections and FRP metrics.
- Join monthly active-fire and FRP GeoTIFFs into 3D R arrays.
- Regrid ERA5 predictors.
- Convert and aggregate MOD13C2 NDVI.
- Regrid TerraClimate VPD.
- Regrid GLEAM surface soil moisture.
- Build monthly FWI95d as the number of days with FWI above the monthly P95.

Main outputs:

```text
longitude.RData
latitude.RData
FireCCI51_2003_2024_0.25degree.RData
FireCCIS311_2019_2024_0.25degree.RData
MODIS-AFcount_conf30_angle30-200301-202412-025.RData
MODIS-FRPsum_conf30_angle30-200301-202412-025.RData
MODIS-FRPmedian_conf30_angle30-200301-202412-025.RData
ERA5-*-2003-2024-MONTLY-025_ADJ.nc
MOD13C2_NDVI_2003_2024_025.nc
vpd_2003-2024_0.25deg_bil.nc
SMs_2003_2024_GLEAM_v4.2b_MO_025deg_bil.nc
fwi-era5_count_exceed95_200301-202412_0.25.nc
```

Detailed instructions:

```text
README-1.txt
```

---

### Step 2. Build fire mask and adjusted land/sea mask

Folder:

```text
2-build-mask
```

Purpose:

Build the monthly MRBA60 fire mask and the fire-adjusted land/sea mask.

Main operations:

- Build probability layers from FireCCI51 and active-fire count using FireCCIS311 as the common-period target.
- Apply ROC/Youden thresholding to filter low-probability active-fire-only detections.
- Generate the final monthly `FireMask_AF3030F`.
- Produce diagnostic annual status maps before and after filtering.
- Adjust the ATLAS/reference land/sea mask by reclassifying ocean cells with observed fire as land.

Main outputs:

```text
FireMask_AF3030F.RData
MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData
ROC_Thresholds_ByBioma_Month_MASK.csv
land_sea_mask_025degree_binary_1440x720_fire_adjusted.RData
land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc
```

Detailed instructions:

```text
README-2.txt
```

---

### Step 3. Apply mask to auxiliary predictors

Folder:

```text
3-apply-mask
```

Purpose:

Ensure that environmental predictors are complete over the valid terrestrial fire domain.

Main operations:

- Apply `FireMask_AF3030F`.
- Apply the fire-adjusted land/sea mask.
- Fill missing predictor values over terrestrial fire pixels using K-nearest neighbours.
- Process NDVI, VPD, SMs and FWI95d.

KNN configuration:

```text
k_nn = 8
max_dist_cells = 8
grid_res_deg = 0.25
```

Main outputs:

```text
NDVI-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
VPD-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
SMs-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN.RData
```

Detailed instructions:

```text
README-3.txt
```

---

### Step 4. Model MRBA60

Folder:

```text
4-model
```

Purpose:

Train and apply the biome- and month-specific harmonisation models.

The modelling block has two stages.

#### 4.1 Predictor selection

Script:

```text
4.1-pre-model-auto-rfe.R
```

Operations:

- Build training datasets by biome and month.
- Use a three-month moving window around each target month.
- Remove collinear predictors using Spearman-correlation clustering.
- Apply Recursive Feature Elimination using random forest.
- Save selected predictors by biome and month.

Main output:

```text
SelectedPredictors_<biome>_COMMON.csv
```

#### 4.2 Final harmonisation

Script:

```text
4.2-Harmonised_MRBA60_second_part_parallel3_LOYO_stats_S3FULL_FINAL.R
```

Operations:

- Read selected predictors from CSV.
- Train final random forest models by biome and month.
- Use FireCCIS311 as the response during 2019-2024.
- Apply models retrospectively to 2003-2024.
- Apply physical constraints to avoid negative values and values above grid-cell area.
- Estimate SHAP values.
- Save biome-level harmonised outputs.

Main outputs:

```text
BA_<Modelo>_<biome>_FireHarmonized_Loyo.RData
BA_<Modelo>_<biome>_FireHarmonized_Common.RData
BA_<Modelo>_<biome>_FireHarmonized_Full.RData
Max_Tope_COMMON_ByMonth_<biome>.csv
```

Detailed instructions:

```text
README-4.txt
```

---

### Step 5. Post-process harmonised burned area

Folder:

```text
5-post
```

Purpose:

Merge biome-level model outputs, restore historical extreme events and apply final empirical quantile mapping.

Main operations:

- Merge biome-level harmonised outputs into global arrays.
- Restore FireCCI51 historical extreme events above biome-month thresholds.
- Apply empirical quantile mapping using positive burned area values only.
- Apply minimum and maximum physical constraints.

Main outputs:

```text
BA_B1-MRBA60-2003-2024global_BA_FireHarmonized_Common.RData
BA_B1-MRBA60-2003-2024global_BA_FireHarmonized_Full.RData
BA_harmonised_correctedByF51Tope_B1-MRBA60-2003-2024.RData
BA_MRBA60.RData
F51_maskAboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData
```

Detailed instructions:

```text
README_5_POST.txt
```

---

### Step 6. Estimate uncertainty

Folder:

```text
6-Incertidumbre
```

Purpose:

Build the final MRBA60 uncertainty layer.

Main operations:

- Extract FireCCIS311 standard error from monthly NetCDF files.
- Estimate uncertainty for the harmonised product from common-period disagreement with FireCCIS311.
- Estimate FireCCI51 uncertainty from common-period disagreement with FireCCIS311.
- Build masks for cells and biome-month combinations where FireCCI51 uncertainty should be used.
- Assemble final MRBA60 uncertainty.

Main outputs:

```text
FireCCIS311_S3_SE_monthly_2019_2024.RData
BA_Incertidumbre_HARMONISED_abs.RData
BA_Incertidumbre_FireCCI51.RData
mask_truefalse_biomas_meses_B1-MRBA60-2003-2024.RData
MASK_NOHARMONISED.RData
BA_Incertidumbre_MRBA60.RData
```

Detailed instructions:

```text
README-6.txt
```

---

### Step 7. Export final NetCDF files

Folder:

```text
7-NetCDF
```

Purpose:

Prepare final burned area and uncertainty arrays in square metres and write the monthly NetCDF files.

Main operations:

- Build final `BA_MRBA60` in m2.
- Build final `Unc_MRBA60` in m2.
- Replace historical values with FireCCI51 where required by the final mask.
- Use FireCCIS311 burned area and standard error for 2019-2024.
- Apply physical minimum and maximum constraints.
- Write one NetCDF file per month with CF/CCI-style metadata.

Main RData outputs:

```text
MRBA60_BA_m2_monthly_2003_2024.RData
MRBA60_Unc_m2_monthly_2003_2024.RData
```

Final NetCDF output directory:

```text
/mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/
```

Detailed instructions:

```text
README-7.txt
```

---

### Step 8. Generate figures and tables for the paper

Folder:

```text
8-Plots paper
```

Purpose:

Generate the main figures, supplementary figures and tables used in the MRBA60 paper.

Main figure groups:

- Global product comparison.
- Validation against ONFIRE, MapBiomas, S2/Landsat-type validation datasets.
- Comparison with FireCCI51, MCD64A1 and GFED5.
- Global maps of mean annual burned area and product differences.
- Annual and seasonal time series.
- Regional GFED/Giglio analyses.
- Predictor-selection summaries.
- Monthly maps of burned area and uncertainty.
- Sensitivity analyses.

Main output root:

```text
/mnt/disco6tb/MRBA60-2/results/D1-Plots/
```

Detailed instructions:

```text
README-8.txt
```

---

## 5. Main directory structure

The scripts use the following project structure.

```text
/mnt/disco6tb/MRBA60/
├── data/
│   ├── A1_RAW/
│   ├── A2_TEMP/
│   └── A3_ADJ/
└── results/
    └── A1-Built-mask/

/mnt/disco6tb/MRBA60-2/
└── results/
    ├── B1-MRBA60-2003-2024/
    │   ├── csv/
    │   ├── plot/
    │   ├── plot_rle/
    │   ├── plot_scatter/
    │   ├── RData/
    │   ├── logs_selected_predictors/
    │   └── logs_harmonisation/
    ├── C1-MRBA60-2003-2024-NetCDF/
    └── D1-Plots/
```

Some scripts still contain older or alternative path roots, especially:

```text
/mnt/disco6tb/MHBA60/
/mnt/disco6tb/MRBA60/results/
/mnt/disco6tb/Dropbox/
```

Before running the full workflow, check whether these paths should be replaced by the current MRBA60 project root.

---

## 6. Required external datasets

The workflow depends on several public and project-specific datasets.

Core burned area datasets:

- FireCCI51, monthly, 2003-2024.
- FireCCIS311, monthly, 2019-2024.

Active fire data:

- MODIS MCD14ML active-fire detections, 2003-2024.

Environmental predictors:

- ERA5 temperature, precipitation, wind speed and cloud cover.
- MOD13C2 NDVI.
- TerraClimate VPD.
- GLEAM surface soil moisture.
- ERA5-derived FWI.

Spatial layers:

- Continental biome polygons.
- ATLAS/reference land/sea mask.

Validation and comparison datasets:

- ONFIRE regional datasets.
- MapBiomas Brazil.
- MCD64A1.
- GFED5.
- S2/Landsat-type validation layers.
- GFED/Giglio region masks.

These datasets are not all included in the repository and must be downloaded or made available locally before execution.

---

## 7. Software requirements

The workflow uses R, shell scripts and CDO/GDAL utilities.

Core command-line tools:

```text
R
CDO
GDAL
bash
```

Main R packages used across the workflow:

```text
ncdf4
terra
raster
sf
sp
dplyr
tidyr
tibble
lubridate
data.table
ggplot2
ggtext
fields
maps
rnaturalearth
rnaturalearthdata
RColorBrewer
viridis
viridisLite
caret
randomForest
fastshap
qmap
RANN
pROC
readr
readxl
openxlsx
patchwork
cowplot
scales
gridExtra
Metrics
forcats
stringr
abind
```

Some scripts use broader package lists than strictly necessary. Install missing packages before running each block.

---

## 8. Reproducibility checklist

Before running the workflow:

- Confirm that all input datasets are available locally.
- Confirm that all directory paths match your local or HPC environment.
- Confirm that `longitude.RData` and `latitude.RData` define a 1440 x 720 global 0.25-degree grid.
- Confirm that all arrays follow the expected order:

```text
[lon, lat, time]
```

- Confirm that the monthly date sequence has 264 layers:

```text
2003-01 to 2024-12
```

- Confirm that FireCCIS311 has 72 layers:

```text
2019-01 to 2024-12
```

- Run each block in the order described above.
- Check logs after each block.
- Check dimensions after every major RData output.
- Plot at least one layer per major output to confirm orientation.
- Check that burned area and uncertainty values do not exceed the 0.25-degree grid-cell area.
- Check that final NetCDF files contain 12 files per year and all expected metadata.

---

## 9. Expected final outputs

The key reproducible outputs are:

```text
MRBA60_BA_m2_monthly_2003_2024.RData
MRBA60_Unc_m2_monthly_2003_2024.RData
YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc
```

Each monthly NetCDF contains:

```text
burned_area
uncertainty
lon
lat
time
lon_bounds
lat_bounds
time_bounds
crs
```

The final monthly NetCDF dataset is stored by year:

```text
C1-MRBA60-2003-2024-NetCDF/
├── 2003/
├── 2004/
├── ...
└── 2024/
```

---

## 10. Validation and paper outputs

The `8-Plots paper` block reproduces the outputs used for the manuscript, including:

- Main figures.
- Supplementary figures.
- Table 2.
- Predictor-selection figures.
- Validation plots.
- Monthly burned area maps.
- Monthly uncertainty maps.
- Sensitivity analyses.

The main figure output root is:

```text
/mnt/disco6tb/MRBA60-2/results/D1-Plots/
```

---

## 11. Recommended citation

Please cite the associated paper when using the MRBA60 scripts or product:

```text
Torres-Vázquez M.A. (2026). MRBA60: Long-term consistent global burned area dataset from Sentinel-3 and MODIS products. MTAV26 GitHub[code]. Date. 
```

Add the final journal citation and DOI here once available.

---

## 12. Contact

For questions about the scripts or processing chain:

```text
Miguel Ángel Torres-Vázquez
miguela.torres@uah.es
```

---

## 13. Notes on reproducibility

The workflow was developed for a Linux/HPC environment with large local storage. Many scripts use absolute paths. For use on a different system, update paths before execution.

The repository is intended to support transparent reproduction of the MRBA60 product and associated analyses. It is not a minimal software package. The scripts are organised to preserve the full research workflow, including diagnostic outputs, intermediate products and manuscript figures.
