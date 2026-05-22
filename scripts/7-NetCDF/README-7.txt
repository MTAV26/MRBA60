Harmonised Medium Resolution Burned Area Grid product version 6.0.0 (MRBA60)

Miguel Ángel Torres-Vázquez
miguela.torres@uah.es


README - 7-NetCDF
=================

This folder contains the final export scripts used to prepare the MRBA60 burned area and uncertainty arrays and write the monthly NetCDF product files for 2003-2024.

The final NetCDF files follow the ESA CCI-style monthly naming convention and include burned area, uncertainty, coordinate bounds, time bounds, CRS information and product-level metadata.


Scripts included
================

  11_BA_2003_2024.R

  11_Incertidumbre_2003_2024.R

  12-Save-NetCDF-Metadata-2003-2024.R


General workflow
================

The workflow is composed of three sequential steps:

  1. Build the final burned area array in square metres.

  2. Build the final uncertainty array in square metres.

  3. Write one NetCDF file per month from January 2003 to December 2024.

The recommended execution order is:

  1. 11_BA_2003_2024.R
  2. 11_Incertidumbre_2003_2024.R
  3. 12-Save-NetCDF-Metadata-2003-2024.R


Model name and main directories
===============================

The scripts use:

  Modelo = "B1-MRBA60-2003-2024"

Main input directory:

  /mnt/disco6tb/MRBA60/data/A3_ADJ/

Main RData output directory:

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/

Final NetCDF output directory:

  /mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/


11_BA_2003_2024.R
=================

Purpose
-------

This script builds the final monthly burned area array for MRBA60 in square metres.

Main operations
---------------

  1. Load the final mask identifying cells/months that should use FireCCI51 values:

       MASK_NOHARMONISED.RData

  2. Load the MRBA60 burned area product:

       BA_MRBA60.RData

  3. Load FireCCI51 for 2003-2024:

       FireCCI51_2003_2024_0.25degree.RData

  4. Replace MRBA60 values with FireCCI51 values where Mask_union_global is TRUE.

  5. Load the FireCCI51 threshold-exceedance mask and threshold-exceedance values:

       F51_maskAboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData
       F51_aboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData

  6. Reapply FireCCI51 threshold-exceedance values for the historical period 2003-2018.

  7. Check whether burned area exceeds the latitude-dependent 0.25-degree grid-cell area.

  8. Cap burned area values larger than the physical maximum cell area.

  9. Set burned area values below the latitude-dependent minimum threshold to zero.

  10. Convert the historical MRBA60 period from km2 to m2.

  11. Append the original FireCCIS311 2019-2024 burned area array.

  12. Save the final burned area array.

Main output
-----------

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_BA_m2_monthly_2003_2024.RData

The saved object is:

  BA_MRBA60

with dimensions:

  [lon, lat, time] = [1440, 720, 264]

Units:

  m2

Temporal structure:

  2003-01 to 2018-12:
    MRBA60 historical harmonised burned area after mask replacement, threshold restoration and physical filtering.

  2019-01 to 2024-12:
    FireCCIS311 burned area.


Physical constraints
--------------------

The script applies two latitude-dependent physical constraints:

  Maximum burned area:
    burned area cannot exceed the area of the 0.25-degree grid cell.

  Minimum burned area:
    very small positive values below the minimum threshold are set to zero.

The maximum cell area is computed as:

  (110.57 * 0.25) * (111.32 * 0.25) * cos(latitude)

The minimum threshold is based on:

  amin_px_eq = 0.09 km2

scaled by cos(latitude).


11_Incertidumbre_2003_2024.R
============================

Purpose
-------

This script builds the final monthly uncertainty array for MRBA60 in square metres.

Main operations
---------------

  1. Load the MRBA60 uncertainty layer:

       BA_Incertidumbre_MRBA60.RData

  2. Load the FireCCIS311 standard error array for 2019-2024:

       FireCCIS311_S3_SE_monthly_2019_2024.RData

  3. Extract the historical uncertainty period 2003-2018 from the MRBA60 uncertainty array.

  4. Convert the historical uncertainty from km2 to m2.

  5. Append the FireCCIS311 standard error layers for 2019-2024.

  6. Convert the combined array back to km2 temporarily for physical filtering.

  7. Cap uncertainty values larger than the latitude-dependent grid-cell area.

  8. Set uncertainty values below the latitude-dependent minimum threshold to zero.

  9. Convert the final uncertainty array to m2.

  10. Save the final uncertainty array.

Main output
-----------

  /mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_Unc_m2_monthly_2003_2024.RData

The saved object is:

  Unc_MRBA60

with dimensions:

  [lon, lat, time] = [1440, 720, 264]

Units:

  m2

Temporal structure:

  2003-01 to 2018-12:
    MRBA60 empirical uncertainty.

  2019-01 to 2024-12:
    FireCCIS311 standard error.


12-Save-NetCDF-Metadata-2003-2024.R
===================================

Purpose
-------

This script writes the final MRBA60 monthly NetCDF files for the complete 2003-2024 period.

Main operations
---------------

  1. Load:

       MRBA60_BA_m2_monthly_2003_2024.RData
       MRBA60_Unc_m2_monthly_2003_2024.RData

  2. Assign:

       burned_area = BA_MRBA60
       uncertainty = Unc_MRBA60

  3. Replace NA values with zero before writing to NetCDF.

  4. Define the 0.25-degree global coordinate vectors:

       lon = -179.875 to 179.875
       lat = -89.875 to 89.875

  5. Define monthly dates from 2003-01-01 to 2024-12-01.

  6. Write one NetCDF file per month.

  7. Include coordinate bounds, time bounds and CRS metadata.

  8. Include global product metadata following CF/CCI conventions.

  9. Perform a quick read/plot check for one output NetCDF file.

Output directory
----------------

  /mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/

The script creates one subdirectory per year:

  2003/
  2004/
  ...
  2024/

NetCDF naming convention
------------------------

Each monthly file is named:

  YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc

Example:

  20190801-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc


NetCDF variables
================

The monthly NetCDF files contain the following variables:

  lon

    Longitude coordinate.

  lat

    Latitude coordinate.

  time

    Monthly time coordinate.

  lon_bounds

    Longitude cell bounds.

  lat_bounds

    Latitude cell bounds.

  time_bounds

    Monthly start and end date bounds.

  crs

    Coordinate reference system variable.

  burned_area

    Monthly burned area.

    Units:
      m2

    Cell method:
      time: sum

    Valid range:
      0 to 769288944

  uncertainty

    Uncertainty of burned area, calculated as Root Mean Square Error.

    Units:
      m2

    Valid range:
      0 to 769288944


Coordinate system
=================

The output grid is global WGS84 latitude-longitude at 0.25-degree resolution.

Longitude:

  -179.875 to 179.875

Latitude:

  -89.875 to 89.875

Grid size:

  1440 x 720

The CRS variable includes a WGS84 WKT definition and the affine transform attribute:

  0.25,0.0,0.0,-0.25,-180.0,90.0


Time encoding
=============

The time coordinate uses:

  days since 1970-01-01 00:00:00

Calendar:

  standard

Each monthly file contains one time step.

The time bounds represent the first and last day of the corresponding month.

Global metadata
===============

The NetCDF files include global metadata such as:

  title
  institution
  source
  history
  references
  tracking_id
  Conventions
  product_version
  format_version
  summary
  keywords
  id
  naming_authority
  doi
  cdm_data_type
  date_created
  creator_name
  creator_url
  creator_email
  contact
  developer_email
  project
  geospatial bounds
  time coverage
  license
  platform
  sensor
  spatial_resolution
  key_variables

The product version written in the NetCDF files is:

  v6.0.0

The CF convention attribute is:

  CF-1.7

The CCI format version attribute is:

  CCI Data Standards v2.3


Important notes
===============

The 11_BA_2003_2024.R script produces the final burned area array in m2. This array is required by the NetCDF writer.

The 11_Incertidumbre_2003_2024.R script produces the final uncertainty array in m2. This array is also required by the NetCDF writer.

The 12-Save-NetCDF-Metadata-2003-2024.R script expects both RData files to exist before execution.

The output NetCDF files replace NA values with zero for both burned_area and uncertainty before writing.

The final product combines:

  2003-2018:
    MRBA60 harmonised historical product.

  2019-2024:
    FireCCIS311 burned area and FireCCIS311 standard error.

Recommended checks
==================

After running 11_BA_2003_2024.R:

  Check that MRBA60_BA_m2_monthly_2003_2024.RData exists.
  Confirm that BA_MRBA60 has dimensions 1440 x 720 x 264.
  Confirm that units are m2.
  Plot several months from both the historical and Sentinel-3 periods.
  Confirm that no burned area value exceeds the grid-cell area.

After running 11_Incertidumbre_2003_2024.R:

  Check that MRBA60_Unc_m2_monthly_2003_2024.RData exists.
  Confirm that Unc_MRBA60 has dimensions 1440 x 720 x 264.
  Confirm that uncertainty values are non-negative.
  Confirm that no uncertainty value exceeds the grid-cell area.

After running 12-Save-NetCDF-Metadata-2003-2024.R:

  Check that yearly folders from 2003 to 2024 were created.
  Confirm that each year contains 12 monthly NetCDF files.
  Open a sample file with nc_open or cdo sinfo.
  Check variables burned_area and uncertainty.
  Check lon_bounds, lat_bounds and time_bounds.
  Plot a sample burned_area layer to confirm orientation.

Example validation commands
===========================

Using CDO:

  cdo sinfo YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc

  cdo showname YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc

  cdo showunit YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc

Using R:

  nc <- ncdf4::nc_open("YYYYMM01-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc")
  names(nc$var)
  ncdf4::nc_close(nc)


End of README.
