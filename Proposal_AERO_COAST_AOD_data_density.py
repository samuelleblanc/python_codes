#!/usr/bin/env python
# coding: utf-8

# # Coastal California AOD Data Density Analysis
# 
# **Name:** coastal_CA_AOD_data_density  
# 
# **Purpose:**  
# Characterise the spatial data density and aerosol optical depth (AOD) sampling coverage
# along the California coast for satellite and airborne platforms used in proposal figures.
# Computes the signed distance of each valid retrieval to the nearest coastline
# (positive = over water, negative = over land) and generates:
# - Maps of satellite pixel locations and aircraft tracks (California overview + Monterey Bay zoom)
# - A histogram of valid-AOD observations as a function of distance to coast
# - Binned aerosol property profiles (AOD, Angstrom exponent) vs. distance to coast
# 
# **Input:**  
# - `/data2/Proposals/Atmosphere_2025/data/PACE/PACE_OCI.*.L2.OC_AOP.V3_2.nc` — PACE OCI L2 ocean color/AOP granules  
# - `/data2/Proposals/Atmosphere_2025/data/VIIRS/AERDB_L2_VIIRS_NOAA20.*.nc` — VIIRS NOAA-20 Deep Blue L2 aerosol  
# - `/data2/Proposals/Atmosphere_2025/data/SpexOne/PACE_SPEXONE.*.L2.RTAP_*.nc` — PACE SPEXone RemoTAP L2 aerosol (land + ocean)  
# - `/data2/Proposals/Atmosphere_2025/data/PACE_PAX/*.ict` — PACE-PAX ER2, CIRPAS Twin Otter, AirSHARP Twin Otter MetNav ICARTT files  
# 
# **Output:** Maps and histograms displayed inline and saved as PNG files  
# 
# **Dependencies:** `numpy`, `scipy`, `matplotlib`, `cartopy`, `shapely`, `geopandas`, `netCDF4`, `pandas`  
# 
# **History:**  
# Written: Samuel LeBlanc, Santa Cruz, CA, 2026-09-01  
# Updated: Added PACE SPEXone RemoTAP data (land + ocean), 2026-09-02

# # Prepare the Python Environment

# In[1]:


import os
import glob
import warnings
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
warnings.filterwarnings('ignore')


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# # Paths and Configuration

# In[3]:


BASE_DIR    = '/data2/Proposals/Atmosphere_2025/data'
PACE_DIR    = os.path.join(BASE_DIR, 'PACE')
VIIRS_DIR   = os.path.join(BASE_DIR, 'VIIRS')
NAV_DIR     = os.path.join(BASE_DIR, 'PACE_PAX')
SPEX_DIR    = os.path.join(BASE_DIR, 'SpexOne')

# California coast bounding box for spatial filtering
CA_LON_MIN, CA_LON_MAX = -130.0, -113.0
CA_LAT_MIN, CA_LAT_MAX =   30.0,   44.0

# Monterey Bay zoom box
MB_LON_MIN, MB_LON_MAX = -123.5, -120.5
MB_LAT_MIN, MB_LAT_MAX =   35.5,   38.0

# PACE OCI L2 bitmask — flag bits to EXCLUDE from ocean-AOD retrieval
# flag_meanings order: ATMFAIL LAND PRODWARN HIGLINT HILT HISATZEN COASTZ SPARE
#   STRAYLIGHT CLDICE COCCOLITH TURBIDW HISOLZEN ... NAVFAIL FILTER ... PRODFAIL
PACE_BAD_FLAGS = (
    0x1        |  # ATMFAIL  (bit 0)
    0x2        |  # LAND     (bit 1)
    0x8        |  # HIGLINT  (bit 3)
    0x10       |  # HILT     (bit 4)
    0x20       |  # HISATZEN (bit 5)
    0x100      |  # STRAYLIGHT (bit 8)
    0x200      |  # CLDICE   (bit 9)
    0x1000     |  # HISOLZEN (bit 12)
    0x2000000  |  # NAVFAIL  (bit 25)
    0x40000000    # PRODFAIL (bit 30)
)
PACE_LAND_BIT = 0x2  # LAND flag used for sign determination

# Aircraft missing/fill value in ICT files
ICT_FILL = -9999.0

# SPEXone RemoTAP fill value
SPEX_FILL = -32767.0

# Colours and styles per dataset
COLORS = {
    'PACE':       '#0077BB',
    'VIIRS':      '#EE7733',
    'SPEXone':    '#228833',
    'ER2':        '#AA3377',
    'CIRPAS_TO':  '#009988',
    'AirSHARP_TO':'#CC3311',
}
LABELS = {
    'PACE':       'PACE OCI AOT$_{865}$',
    'VIIRS':      'VIIRS DB AOT$_{550}$',
    'SPEXone':    'PACE SPEXone RemoTAP AOT$_{550}$',
    'ER2':        'ER-2 HSRL AOD$_{532}$',
    'CIRPAS_TO':  'CIRPAS Twin Otter aerosol ext$_{532}$',
    'AirSHARP_TO':'AirSHARP Twin Otter 4STAR AOD',
}


# # Load files

# ## ICARTT (NASA Ames 1001) File Reader
# 
# All aircraft navigation files use the NASA Ames ICARTT format (FFI 1001).  
# The first line gives `NLHEAD, FFI`; data begins on line `NLHEAD + 1`.
# The last header line contains comma-separated column names.

# In[4]:


def read_ict(filepath, fill_val=ICT_FILL):
    """Read a NASA Ames ICARTT 1001 (.ict) file into a pandas DataFrame.
    
    Header lines (NLHEAD) are skipped.  The second-to-last skipped line is
    the column-name row.  Fill values are replaced with NaN.
    """
    with open(filepath, 'r') as fh:
        first = fh.readline()
    nlhead = int(first.split(',')[0].strip())
    df = pd.read_csv(filepath, skiprows=nlhead - 1, sep=',',
                     skipinitialspace=True, low_memory=False)
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    # Replace fill values with NaN (allow for both integer and float fill)
    df.replace(fill_val, np.nan, inplace=True)
    df.replace(int(fill_val), np.nan, inplace=True)
    return df


# ## Load Aircraft Navigation Data
# 
# Three platforms:
# - **ER-2** (NASA high-altitude): `PACEPAX-MetNav_ER2_*.ict`  
# - **CIRPAS Twin Otter** (low-altitude, PACE-PAX): `PACEPAX-MetNav_CIRPAS-TO_*.ict`  
# - **AirSHARP Twin Otter** (low-altitude, AirSHARP): `AirSHARP-MetNav-1Hz_AirSHARP-TO_*.ict`  
# 
# Files with `_L1` / `_L2` suffixes are separate flight legs for the same day — all are loaded.
# Common columns extracted: `Latitude`, `Longitude`, `Time_Start`.

# In[5]:


def load_aircraft_group(pattern, label):
    """Load all ICT files matching glob pattern; return concatenated DataFrame."""
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'  {label}: no files found for {pattern}')
        return pd.DataFrame(columns=['Latitude', 'Longitude', 'Time_Start', 'filename'])
    dfs = []
    for fp in files:
        df = read_ict(fp)
        df['filename'] = os.path.basename(fp)
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)
    # Filter to valid lat/lon in the California domain
    valid = (
        combined['Latitude'].between(CA_LAT_MIN, CA_LAT_MAX) &
        combined['Longitude'].between(CA_LON_MIN, CA_LON_MAX)
    )
    combined = combined[valid].copy().reset_index(drop=True)
    print(f'  {label}: {len(files)} files → {len(combined):,} valid records')
    return combined

print('Loading aircraft navigation data...')
df_er2     = load_aircraft_group(os.path.join(NAV_DIR, 'PACEPAX-MetNav_ER2_*.ict'), 'ER-2')
df_cirpas  = load_aircraft_group(os.path.join(NAV_DIR, 'PACEPAX-MetNav_CIRPAS-TO_*.ict'), 'CIRPAS Twin Otter')
df_airsharp= load_aircraft_group(os.path.join(NAV_DIR, 'AirSHARP-MetNav-1Hz_AirSHARP-TO_2024*.ict'), 'AirSHARP Twin Otter')


# ## Load the CIRPAS optical

# In[6]:


from scipy import interpolate


# In[7]:


def load_aircraft_group_NOLAT(pattern, label):
    """Load all ICT files matching glob pattern; return concatenated DataFrame."""
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'  {label}: no files found for {pattern}')
        return pd.DataFrame(columns=['Latitude', 'Longitude', 'Time_Start', 'filename'])
    dfs = []
    for fp in files:
        df = read_ict(fp)
        df['filename'] = os.path.basename(fp)
        # check the navmet to interpolate
        try:
            dn = read_ict(fp.replace('LARGE-OPTICAL','MetNav'))
            ffxlat = interpolate.interp1d(dn['Time_Start'], dn['Latitude'], bounds_error=False, fill_value=np.nan)
            ffxlon = interpolate.interp1d(dn['Time_Start'], dn['Longitude'], bounds_error=False, fill_value=np.nan)
            df['Latitude'] = ffxlat(df['Time_Mid'])
            df['Longitude'] = ffxlon(df['Time_Mid'])
        except:
            print(f'Error interpolating for file: {fp}')
        dfs.append(df)
    combined = pd.concat(dfs, ignore_index=True)
    # Filter to valid lat/lon in the California domain
    #valid = (
    #    combined['Latitude'].between(CA_LAT_MIN, CA_LAT_MAX) &
    #    combined['Longitude'].between(CA_LON_MIN, CA_LON_MAX)
    #)
    #combined = combined[valid].copy().reset_index(drop=True)
    print(f'  {label}: {len(files)} files → {len(combined):,} valid records')
    return combined


# In[8]:


#PACEPAX-LARGE-OPTICAL_CIRPAS-TO_20240927_R0.ict
df_optical  = load_aircraft_group_NOLAT(os.path.join(NAV_DIR, 'PACEPAX-LARGE-OPTICAL_CIRPAS-TO_*.ict'), 'CIRPAS Twin Otter')


# In[9]:


df_optical['Latitude'][~np.isfinite(df_optical['Ext532_total_amb'])]=np.nan
df_optical['Longitude'][~np.isfinite(df_optical['Ext532_total_amb'])]=np.nan


# In[10]:


df_optical.keys()


# ## Load the 4STAR netcdf aircraft data

# In[11]:


STAR_DIR   = '/data/sunsat/AirSHARP_2024/data_archival/AOD_R0_nc'

star_files = sorted(glob.glob(os.path.join(STAR_DIR, '*.nc')))
print(f'Found {len(star_files)} 4STAR files')

star_lat_all  = []
star_lon_all  = []
star_aot_all  = []
star_ang_all  = []   # Angstrom exponent 443–865 nm
star_lat_withbad  = []
star_lon_withbad  = []

for fp in star_files:
    with nc.Dataset(fp) as ds:
        lat  = ds['Latitude'][:]
        lon  = ds['Longitude'][:]
        aot  = ds['AOD'][:]   # masked array
        flags= np.array(ds['qual_flag'][:], dtype=np.int8)
        ang  = ds['AOD_angstrom_470_865'][:]  # masked array

    # Spatial filter to California region
    in_box = (
        (lat >= CA_LAT_MIN) & (lat <= CA_LAT_MAX) &
        (lon >= CA_LON_MIN) & (lon <= CA_LON_MAX)
    )
    if not in_box.any():
        continue

    # Quality filter: exclude pixels with any bad flag set
    bad = (flags).astype(bool)
    # Also exclude masked (fill) values in aot
    aot_fill = np.ma.getmaskarray(aot).any(axis=0)
    good = in_box & ~bad & ~aot_fill

    if good.any():
        star_lat_all.append(np.array(lat[good]))
        star_lon_all.append(np.array(lon[good]))
        star_aot_all.append(np.ma.filled(aot[:,good], np.nan))
        star_ang_all.append(np.ma.filled(ang[good], np.nan))
    if in_box.any():
        star_lat_withbad.append(np.array(lat[in_box]))
        star_lon_withbad.append(np.array(lon[in_box]))

star_lat  = np.concatenate(star_lat_all)  if star_lat_all  else np.array([])
star_lon  = np.concatenate(star_lon_all)  if star_lon_all  else np.array([])
star_latwb  = np.concatenate(star_lat_withbad)  if star_lat_all  else np.array([])
star_lonwb  = np.concatenate(star_lon_withbad)  if star_lon_all  else np.array([])
star_aot  = np.concatenate(star_aot_all,axis=1)  if star_aot_all  else np.array([])
star_ang  = np.concatenate(star_ang_all)  if star_ang_all  else np.array([])

print(f'4STAR valid pixels in CA domain: {len(star_lat):,} with the bad points: {len(star_latwb):,}')
print(f'  AOT_865 range: {np.nanmin(star_aot):.3f} – {np.nanmax(star_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(star_ang):.3f} – {np.nanmax(star_ang):.3f}')


# In[12]:


ds = nc.Dataset(fp)


# In[13]:


ds['wavelength'][:][192]


# ## Load the ER2 data (for aod flags)

# In[14]:


HSRL_DIR = '/data2/PACEPAX/ER2'
hsrl_files = sorted(glob.glob(os.path.join(HSRL_DIR, '*.h5')))
print(f'Found {len(hsrl_files)} HSRL files')

hsrl_lat_all  = []
hsrl_lon_all  = []
hsrl_aot_all  = []
hsrl_ang_all  = []   # Angstrom exponent 443–865 nm
hsrl_lat_withbad  = []
hsrl_lon_withbad  = []

for fp in hsrl_files:
    with nc.Dataset(fp) as ds:
        lat  = ds['Nav_Data/gps_lat'][:]
        lon  = ds['Nav_Data/gps_lon'][:]
        aot  = ds['DataProducts/532_AOT_from_bsc'][:]   # masked array
        flags= ds['DataUncertainty/532_AOT_hi_unc'][:]
        ang  = ds['DataProducts/Angstrom_532_355'][:].mean(axis=1)  # masked array

    # Spatial filter to California region
    in_box = (
        (lat >= CA_LAT_MIN) & (lat <= CA_LAT_MAX) &
        (lon >= CA_LON_MIN) & (lon <= CA_LON_MAX)
    )
    if not in_box.any():
        continue

    # Quality filter: exclude pixels with any bad flag set
    bad = flags.mask
    # Also exclude masked (fill) values in aot
    aot_fill = np.ma.getmaskarray(aot)
    good = in_box & ~bad & ~aot_fill

    if good.any():
        hsrl_lat_all.append(np.array(lat[good]))
        hsrl_lon_all.append(np.array(lon[good]))
        hsrl_aot_all.append(np.ma.filled(aot[good], np.nan))
        hsrl_ang_all.append(np.ma.filled(ang[good[:,0]], np.nan))
    if in_box.any():
        hsrl_lat_withbad.append(np.array(lat[in_box]))
        hsrl_lon_withbad.append(np.array(lon[in_box]))

hsrl_lat  = np.concatenate(hsrl_lat_all)  if hsrl_lat_all  else np.array([])
hsrl_lon  = np.concatenate(hsrl_lon_all)  if hsrl_lon_all  else np.array([])
hsrl_latwb  = np.concatenate(hsrl_lat_withbad)  if hsrl_lat_all  else np.array([])
hsrl_lonwb  = np.concatenate(hsrl_lon_withbad)  if hsrl_lon_all  else np.array([])
hsrl_aot  = np.concatenate(hsrl_aot_all)  if hsrl_aot_all  else np.array([])
hsrl_ang  = np.concatenate(hsrl_ang_all)  if hsrl_ang_all  else np.array([])

print(f'HSRL valid pixels in CA domain: {len(hsrl_lat):,} with the bad points: {len(hsrl_latwb):,}')
print(f'  AOT_532 range: {np.nanmin(hsrl_aot):.3f} – {np.nanmax(hsrl_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(hsrl_ang):.3f} – {np.nanmax(hsrl_ang):.3f}')


# ## Load PACE OCI L2 Satellite Data
# 
# Variables used:
# - `navigation_data/latitude`, `navigation_data/longitude`
# - `geophysical_data/aot_865` — AOT at 865 nm (masked array, scale already applied by netCDF4)
# - `geophysical_data/l2_flags` — bitmask; apply `PACE_BAD_FLAGS` to exclude poor-quality pixels
# 
# Only ocean pixels outside bad-flag conditions are retained.  
# Spatial subset: California bounding box.

# In[15]:


pace_files = sorted(glob.glob(os.path.join(PACE_DIR, 'PACE_OCI.*.L2.OC_AOP.V3_2.nc')))
print(f'Found {len(pace_files)} PACE OCI granules')

pace_lat_all  = []
pace_lon_all  = []
pace_aot_all  = []
pace_land_all = []   # True if pixel is over land (for sign of distance)
pace_ang_all  = []   # Angstrom exponent 443–865 nm
pace_lat_withbad  = []
pace_lon_withbad  = []
pace_land_withbad = [] 

for fp in pace_files:
    with nc.Dataset(fp) as ds:
        lat  = ds['navigation_data/latitude'][:]
        lon  = ds['navigation_data/longitude'][:]
        aot  = ds['geophysical_data/aot_865'][:]   # masked array
        flags= np.array(ds['geophysical_data/l2_flags'][:], dtype=np.int32)
        ang  = ds['geophysical_data/angstrom'][:]  # masked array

    # Spatial filter to California region
    in_box = (
        (lat >= CA_LAT_MIN) & (lat <= CA_LAT_MAX) &
        (lon >= CA_LON_MIN) & (lon <= CA_LON_MAX)
    )
    if not in_box.any():
        continue

    # Quality filter: exclude pixels with any bad flag set
    bad = (flags & PACE_BAD_FLAGS).astype(bool)
    # Also exclude masked (fill) values in aot
    aot_fill = np.ma.getmaskarray(aot)
    good = in_box & ~bad & ~aot_fill

    if good.any():
        pace_lat_all.append(np.array(lat[good]))
        pace_lon_all.append(np.array(lon[good]))
        pace_aot_all.append(np.ma.filled(aot[good], np.nan))
        pace_ang_all.append(np.ma.filled(ang[good], np.nan))
        is_land = (flags[good] & PACE_LAND_BIT).astype(bool)
        pace_land_all.append(is_land)
    if in_box.any():
        pace_lat_withbad.append(np.array(lat[in_box]))
        pace_lon_withbad.append(np.array(lon[in_box]))
        is_landwb = (flags[in_box] & PACE_LAND_BIT).astype(bool)
        pace_land_withbad.append(is_landwb)

pace_lat  = np.concatenate(pace_lat_all)  if pace_lat_all  else np.array([])
pace_lon  = np.concatenate(pace_lon_all)  if pace_lon_all  else np.array([])
pace_latwb  = np.concatenate(pace_lat_withbad)  if pace_lat_all  else np.array([])
pace_lonwb  = np.concatenate(pace_lon_withbad)  if pace_lon_all  else np.array([])
pace_aot  = np.concatenate(pace_aot_all)  if pace_aot_all  else np.array([])
pace_ang  = np.concatenate(pace_ang_all)  if pace_ang_all  else np.array([])
pace_land = np.concatenate(pace_land_all) if pace_land_all else np.array([], dtype=bool)
pace_landwb = np.concatenate(pace_land_withbad) if pace_land_withbad else np.array([], dtype=bool)

print(f'PACE valid pixels in CA domain: {len(pace_lat):,} with the bad points: {len(pace_latwb):,}')
print(f'  AOT_865 range: {np.nanmin(pace_aot):.3f} – {np.nanmax(pace_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(pace_ang):.3f} – {np.nanmax(pace_ang):.3f}')


# ## Load VIIRS NOAA-20 Deep Blue L2 Satellite Data
# 
# Variables used:
# - `Latitude`, `Longitude`
# - `Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate` — already QA-filtered (QA ≥ 2 land, QA = 3 ocean)
# - `Angstrom_Exponent_Land_Ocean_Best_Estimate`
# - `Aerosol_Optical_Thickness_550_Ocean` — used to determine sign (valid → over water)
# 
# Fill value: −999.  Spatial subset: California bounding box.

# In[16]:


viirs_files = sorted(glob.glob(os.path.join(VIIRS_DIR, 'AERDB_L2_VIIRS_NOAA20.*.nc')))
print(f'Found {len(viirs_files)} VIIRS granules')

viirs_lat_all   = []
viirs_lon_all   = []
viirs_lat_withbad   = []
viirs_lon_withbad   = []
viirs_aot_all   = []
viirs_ang_all   = []
viirs_ocean_all = []  # True if pixel is over ocean (sign determination)
viirs_ocean_withbad = []

VIIRS_FILL = -999.0

for fp in viirs_files:
    with nc.Dataset(fp) as ds:
        lat  = np.array(ds['Latitude'][:])
        lon  = np.array(ds['Longitude'][:])
        aot  = np.array(ds['Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'][:])
        ang  = np.array(ds['Angstrom_Exponent_Land_Ocean_Best_Estimate'][:])
        aot_ocn = np.array(ds['Aerosol_Optical_Thickness_550_Ocean'][:])
        qa_land = np.array(ds['Aerosol_Optical_Thickness_QA_Flag_Land'][:])
        qa_ocn  = np.array(ds['Aerosol_Optical_Thickness_QA_Flag_Ocean'][:])

    # Replace fills with NaN
    aot    = np.where(aot    <= VIIRS_FILL, np.nan, aot)
    ang    = np.where(ang    <= VIIRS_FILL, np.nan, ang)
    aot_ocn= np.where(aot_ocn<= VIIRS_FILL, np.nan, aot_ocn)

    # Spatial + valid filter
    in_box = (
        (lat >= CA_LAT_MIN) & (lat <= CA_LAT_MAX) &
        (lon >= CA_LON_MIN) & (lon <= CA_LON_MAX)
    )
    good = in_box & np.isfinite(aot)

    if good.any():
        viirs_lat_all.append(lat[good])
        viirs_lon_all.append(lon[good])
        viirs_aot_all.append(aot[good])
        viirs_ang_all.append(ang[good])
        # Ocean pixel if ocean QA is valid (1 = poor but retrieved, 3 = good)
        is_ocean = np.isfinite(aot_ocn[good])
        viirs_ocean_all.append(is_ocean)
    if in_box.any():
        viirs_lat_withbad.append(lat[in_box])
        viirs_lon_withbad.append(lon[in_box])
        is_oceanwb = np.isfinite(aot_ocn[in_box])
        viirs_ocean_withbad.append(is_oceanwb)

viirs_lat   = np.concatenate(viirs_lat_all)   if viirs_lat_all   else np.array([])
viirs_lon   = np.concatenate(viirs_lon_all)   if viirs_lon_all   else np.array([])
viirs_latwb = np.concatenate(viirs_lat_withbad)   if viirs_lat_withbad   else np.array([])
viirs_lonwb = np.concatenate(viirs_lon_withbad)   if viirs_lon_withbad   else np.array([])
viirs_aot   = np.concatenate(viirs_aot_all)   if viirs_aot_all   else np.array([])
viirs_ang   = np.concatenate(viirs_ang_all)   if viirs_ang_all   else np.array([])
viirs_ocean = np.concatenate(viirs_ocean_all) if viirs_ocean_all else np.array([], dtype=bool)
viirs_oceanwb = np.concatenate(viirs_ocean_withbad) if viirs_ocean_withbad else np.array([], dtype=bool)


print(f'VIIRS valid pixels in CA domain: {len(viirs_lat):,}, with bad: {len(viirs_latwb):,}')
print(f'  AOT_550 range: {np.nanmin(viirs_aot):.3f} – {np.nanmax(viirs_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(viirs_ang):.3f} – {np.nanmax(viirs_ang):.3f}')


# ## Load PACE SPEXone RemoTAP L2 Data
# 
# PACE SPEXone multi-angle polarimeter retrievals from the RemoTAP algorithm.
# Separate netCDF files are provided for land (`RTAP_LAND`) and ocean (`RTAP_OCEAN`)
# retrievals — the surface type is therefore known directly from the filename.
# 
# Variables used:
# - `geolocation_data/latitude`, `geolocation_data/longitude`
# - `geophysical_data/aot550` — AOT at 550 nm; fill = −32767
# - `geophysical_data/angstrom_440_870` — Ångström exponent 440–870 nm
# 
# No explicit QA flag variable; invalid pixels have fill value.
# Sign: RTAP_LAND → negative (over land), RTAP_OCEAN → positive (over ocean).

# In[17]:


spex_lat_all   = []
spex_lon_all   = []
spex_lat_withbad   = []
spex_lon_withbad   = []
spex_aot_all   = []
spex_ang_all   = []
spex_land_all  = []  # True = over land (sign = negative)
spex_land_withbad  = []

for surf, is_land_flag in [('OCEAN', False), ('LAND', True)]:
    files = sorted(glob.glob(os.path.join(SPEX_DIR, f'PACE_SPEXONE.*.RTAP_{surf}.*.nc')))
    print(f'SPEXone {surf}: {len(files)} files')
    for fp in files:
        with nc.Dataset(fp) as ds:
            lat = np.ma.filled(ds['geolocation_data/latitude'][:],  np.nan)
            lon = np.ma.filled(ds['geolocation_data/longitude'][:], np.nan)
            aot = np.where(np.array(ds['geophysical_data/aot550'][:])          <= SPEX_FILL,
                           np.nan, np.array(ds['geophysical_data/aot550'][:]))
            ang = np.where(np.array(ds['geophysical_data/angstrom_440_870'][:]) <= SPEX_FILL,
                           np.nan, np.array(ds['geophysical_data/angstrom_440_870'][:]))

        in_box = (
            (lat >= CA_LAT_MIN) & (lat <= CA_LAT_MAX) &
            (lon >= CA_LON_MIN) & (lon <= CA_LON_MAX)
        )
        good = in_box & np.isfinite(aot)
        if good.any():
            spex_lat_all.append(lat[good])
            spex_lon_all.append(lon[good])
            spex_aot_all.append(aot[good])
            spex_ang_all.append(ang[good])
            spex_land_all.append(np.full(good.sum(), is_land_flag, dtype=bool))
        if in_box.any():
            spex_lat_withbad.append(lat[in_box])
            spex_lon_withbad.append(lon[in_box])
            spex_land_withbad.append(np.full(in_box.sum(), is_land_flag, dtype=bool))

spex_lat  = np.concatenate(spex_lat_all)  if spex_lat_all  else np.array([])
spex_lon  = np.concatenate(spex_lon_all)  if spex_lon_all  else np.array([])
spex_latwb  = np.concatenate(spex_lat_withbad)  if spex_lat_withbad  else np.array([])
spex_lonwb  = np.concatenate(spex_lon_withbad)  if spex_lon_withbad  else np.array([])
spex_aot  = np.concatenate(spex_aot_all)  if spex_aot_all  else np.array([])
spex_ang  = np.concatenate(spex_ang_all)  if spex_ang_all  else np.array([])
spex_land = np.concatenate(spex_land_all) if spex_land_all else np.array([], dtype=bool)
spex_landwb = np.concatenate(spex_land_withbad) if spex_land_withbad else np.array([], dtype=bool)

print(f'\nSPEXone total valid pixels in CA domain: {len(spex_lat):,}, with bad: {len(spex_latwb):,}')
print(f'  over ocean: {(~spex_land).sum():,}  |  over land: {spex_land.sum():,}')
if len(spex_aot) > 0:
    print(f'  AOT_550 range: {np.nanmin(spex_aot):.3f} – {np.nanmax(spex_aot):.3f}')
    print(f'  Angstrom (440-870) range: {np.nanmin(spex_ang):.3f} – {np.nanmax(spex_ang):.3f}')


# # Build Distance-to-Coast Tools
# 
# **Magnitude** — KD-tree built from Natural Earth 10 m coastline, densified to ≈0.005°
# spacing (≈500 m), converted to km using a flat-Earth approximation centred on the
# California domain.  
# 
# **Sign** — Land mask rasterised at 0.02° (≈2 km) onto the California grid using
# `matplotlib.path.Path.contains_points`.  Sign is −1 if the point is over land,
# +1 if over water.  The sign from the satellite flag bits is used directly for PACE
# and VIIRS pixels (faster and consistent with each product's own surface classification).

# In[18]:


def build_distance_tools(bounds=(CA_LON_MIN, CA_LON_MAX, CA_LAT_MIN, CA_LAT_MAX),
                         resolution='10m'):
    """Build KD-tree and land-mask raster for signed distance-to-coast computation.
    
    Returns
    -------
    kd        : scipy cKDTree, queryable in km-projected coordinates
    scale     : array([km_per_deg_lon, km_per_deg_lat]) for coordinate conversion
    is_land_grid : bool 2-D array (lat × lon) at 0.02° resolution
    lon_r, lat_r : 1-D coordinate arrays for the land-mask grid
    """
    lon_min, lon_max, lat_min, lat_max = bounds
    lat0  = (lat_min + lat_max) / 2.0
    scale = np.array([111.32 * np.cos(np.radians(lat0)), 111.32])  # km / degree
    box_geom = sgeom.box(lon_min - 1, lat_min - 1, lon_max + 1, lat_max + 1)

    # --- Coastline → densified points → KD-tree ---
    print('  Building coastline KD-tree...')
    coast_path = shpreader.natural_earth(resolution, 'physical', 'coastline')
    raw_segs = []
    for geom in shpreader.Reader(coast_path).geometries():
        if not geom.intersects(box_geom):
            continue
        segs = [geom] if geom.geom_type == 'LineString' else list(geom.geoms)
        for seg in segs:
            coords = np.array(seg.coords)
            # Densify each segment to ~0.005° spacing
            for k in range(len(coords) - 1):
                d_deg = np.hypot(*(coords[k + 1] - coords[k]))
                n = max(2, int(d_deg / 0.005) + 1)
                raw_segs.append(np.linspace(coords[k], coords[k + 1], n))
    if not raw_segs:
        raise RuntimeError('No coastline found in bounds')
    coast_pts = np.vstack(raw_segs)  # shape (N, 2) lon/lat
    kd = cKDTree(coast_pts * scale)  # in km
    print(f'    {len(coast_pts):,} coastline points in KD-tree')

    # --- Land polygons → raster mask ---
    print('  Building land-mask raster...')
    lon_r = np.arange(lon_min, lon_max + 0.02, 0.02)
    lat_r = np.arange(lat_min, lat_max + 0.02, 0.02)
    LO, LA = np.meshgrid(lon_r, lat_r)
    grid_pts = np.column_stack([LO.ravel(), LA.ravel()])
    is_land_flat = np.zeros(len(grid_pts), dtype=bool)

    land_path = shpreader.natural_earth(resolution, 'physical', 'land')
    for geom in shpreader.Reader(land_path).geometries():
        if not geom.intersects(box_geom):
            continue
        polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
        for poly in polys:
            ext = np.array(poly.exterior.coords)
            p = MplPath(ext)
            is_land_flat |= p.contains_points(grid_pts)

    is_land_grid = is_land_flat.reshape(LA.shape)
    print(f'    Land mask: {LA.shape} grid ({is_land_flat.sum():,} / {len(grid_pts):,} land cells)')
    return kd, scale, is_land_grid, lon_r, lat_r


print('Building distance tools (one-time setup)...')
kd_coast, scale_km, is_land_grid, lon_r, lat_r = build_distance_tools()
print('Done.')


# In[19]:


def dist_to_coast_km(lats, lons, sign_from_flag=None):
    """Compute signed distance to nearest coastline in km.
    
    Parameters
    ----------
    lats, lons       : array-like of geographic coordinates
    sign_from_flag   : bool array, True = over land (negative distance).
                       If None, determined from the land-mask raster.
    
    Returns
    -------
    dist : float array, km.  Positive = over water, negative = over land.
    """
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)
    valid = np.isfinite(lats) & np.isfinite(lons)

    dist = np.full(len(lats), np.nan)
    if valid.any():
        pts_km = np.column_stack([lons[valid], lats[valid]]) * scale_km
        dist[valid], _ = kd_coast.query(pts_km)

    # Determine sign
    if sign_from_flag is not None:
        is_land = np.asarray(sign_from_flag, dtype=bool)
    else:
        # Lookup in raster via nearest grid index
        lon_idx = np.searchsorted(lon_r, lons, side='right') - 1
        lat_idx = np.searchsorted(lat_r, lats, side='right') - 1
        in_grid = (
            valid &
            (lon_idx >= 0) & (lon_idx < len(lon_r)) &
            (lat_idx >= 0) & (lat_idx < len(lat_r))
        )
        is_land = np.zeros(len(lats), dtype=bool)
        is_land[in_grid] = is_land_grid[lat_idx[in_grid], lon_idx[in_grid]]

    dist[is_land] *= -1
    return dist


# # Compute Distance to Coast for All Datasets

# In[20]:


print('Computing distance to coast...')

# PACE: use LAND flag bit for sign (consistent with own retrieval surface class)
pace_dist = dist_to_coast_km(pace_lat, pace_lon, sign_from_flag=pace_land)
pace_distwb = dist_to_coast_km(pace_latwb, pace_lonwb, sign_from_flag=pace_landwb)
print(f'  PACE:     {np.sum(np.isfinite(pace_dist)):,} valid, '
      f'range {np.nanmin(pace_dist):.1f} – {np.nanmax(pace_dist):.1f} km')

# VIIRS: ocean pixels → positive, land pixels → negative
viirs_sign_land = ~viirs_ocean  # land if ocean retrieval is not valid
viirs_dist = dist_to_coast_km(viirs_lat, viirs_lon, sign_from_flag=viirs_sign_land)
viirs_distwb = dist_to_coast_km(viirs_latwb, viirs_lonwb, sign_from_flag=~viirs_oceanwb)
print(f'  VIIRS:    {np.sum(np.isfinite(viirs_dist)):,} valid, '
      f'range {np.nanmin(viirs_dist):.1f} – {np.nanmax(viirs_dist):.1f} km')

# SPEXone: sign comes directly from RTAP_LAND / RTAP_OCEAN filename
spex_dist = dist_to_coast_km(spex_lat, spex_lon, sign_from_flag=spex_land)
spex_distwb = dist_to_coast_km(spex_latwb, spex_lonwb, sign_from_flag=spex_landwb)
if len(spex_dist) > 0:
    print(f'  SPEXone:  {np.sum(np.isfinite(spex_dist)):,} valid, '
          f'range {np.nanmin(spex_dist):.1f} – {np.nanmax(spex_dist):.1f} km')

# Aircraft: use land-mask raster (no instrument flag available for surface type)
#er2_dist     = dist_to_coast_km(df_er2['Latitude'],      df_er2['Longitude'])

cirpas_dist  = dist_to_coast_km(df_cirpas['Latitude'],   df_cirpas['Longitude'])
#airsharp_dist= dist_to_coast_km(df_airsharp['Latitude'], df_airsharp['Longitude'])
#airsharp_distwb= dist_to_coast_km(star_latwb, star_lonwb)
print(f'  ER-2:        {np.sum(np.isfinite(er2_dist)):,} valid')
print(f'  CIRPAS TO:   {np.sum(np.isfinite(cirpas_dist)):,} valid')
print(f'  AirSHARP TO: {np.sum(np.isfinite(airsharp_dist)):,} valid')

print('Done.')


# In[21]:


er2_dist     = dist_to_coast_km(hsrl_lat,      hsrl_lon)
er2_distwb     = dist_to_coast_km(hsrl_latwb,      hsrl_lonwb)
airsharp_dist= dist_to_coast_km(star_lat, star_lon)
airsharp_distwb= dist_to_coast_km(star_latwb, star_lonwb)
print(f'  CIRPAS TO:   {np.sum(np.isfinite(cirpas_dist)):,} valid')


# In[ ]:


cirpas_dist  = dist_to_coast_km(df_optical['Latitude'],   df_optical['Longitude'])
cirpas_distwb  = dist_to_coast_km(df_cirpas['Latitude'],   df_cirpas['Longitude'])


# # Make the figures

# ## Figure 1 – California Coast Overview Map
# 
# Satellite valid-AOD pixel centroids (subsampled for clarity) and all aircraft tracks
# overlaid on a map of the California coast.

# In[ ]:


matplotlib.rcParams['text.usetex'] = False 


# In[ ]:


# Subsample satellite pixels for plotting (every N-th point)
PACE_STRIDE  = max(1, len(pace_lat)  // 50000)
VIIRS_STRIDE = max(1, len(viirs_lat) // 50000)
SPEX_STRIDE  = max(1, len(spex_lat)  // 10000)

fig = plt.figure(figsize=(11, 10))
ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([CA_LON_MIN, CA_LON_MAX, CA_LAT_MIN, CA_LAT_MAX], crs=ccrs.PlateCarree())

# Background
ax.add_feature(cfeature.OCEAN.with_scale('10m'), color='#d0e8f5', zorder=0)
ax.add_feature(cfeature.LAND.with_scale('10m'),  color='#e8e0d5', zorder=1)
ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='gray', linewidth=0.6, zorder=2)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', linewidth=0.8, zorder=3)

# PACE pixels
if len(pace_lat) > 0:
    ax.scatter(pace_lon[::PACE_STRIDE], pace_lat[::PACE_STRIDE],
               c=COLORS['PACE'], s=1, alpha=0.25, transform=ccrs.PlateCarree(),
               zorder=4, rasterized=True)

# VIIRS pixels
if len(viirs_lat) > 0:
    ax.scatter(viirs_lon[::VIIRS_STRIDE], viirs_lat[::VIIRS_STRIDE],
               c=COLORS['VIIRS'], s=3, alpha=0.25, marker='s',
               transform=ccrs.PlateCarree(), zorder=4, rasterized=True)

# SPEXone pixels
if len(spex_lat) > 0:
    ax.scatter(spex_lon[::SPEX_STRIDE], spex_lat[::SPEX_STRIDE],
               c=COLORS['SPEXone'], s=5, alpha=0.5, marker='^',
               transform=ccrs.PlateCarree(), zorder=4, rasterized=True)

# Aircraft tracks
if len(df_er2) > 0:
    ax.plot(df_er2['Longitude'].values, df_er2['Latitude'].values,
            color=COLORS['ER2'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), zorder=5)
if len(df_cirpas) > 0:
    ax.plot(df_cirpas['Longitude'].values, df_cirpas['Latitude'].values,
            color=COLORS['CIRPAS_TO'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), zorder=5)
if len(df_airsharp) > 0:
    ax.plot(df_airsharp['Longitude'].values, df_airsharp['Latitude'].values,
            color=COLORS['AirSHARP_TO'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), zorder=5)

# Grid lines
gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
gl.top_labels = False
gl.right_labels = False

# Legend
handles = []
if len(pace_lat)  > 0: handles.append(Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['PACE'],     markersize=6, label=LABELS['PACE']))
if len(viirs_lat) > 0: handles.append(Line2D([0],[0], marker='s', color='w', markerfacecolor=COLORS['VIIRS'],    markersize=6, label=LABELS['VIIRS']))
if len(spex_lat)  > 0: handles.append(Line2D([0],[0], marker='^', color='w', markerfacecolor=COLORS['SPEXone'],  markersize=6, label=LABELS['SPEXone']))
if len(df_er2)    > 0: handles.append(Line2D([0],[0], color=COLORS['ER2'],        lw=1.5, label=LABELS['ER2']))
if len(df_cirpas) > 0: handles.append(Line2D([0],[0], color=COLORS['CIRPAS_TO'],  lw=1.5, label=LABELS['CIRPAS_TO']))
if len(df_airsharp)>0: handles.append(Line2D([0],[0], color=COLORS['AirSHARP_TO'],lw=1.5, label=LABELS['AirSHARP_TO']))
ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.85)

ax.set_title('Satellite AOD Pixel Coverage and Aircraft Tracks\nCalifornia Coast', fontsize=13)
plt.tight_layout()
plt.savefig('coastal_CA_map_overview.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_map_overview.png')


# In[ ]:


print(matplotlib.__file__)


# ## Figure 2 – Monterey Bay Zoom Map

# In[ ]:


fig = plt.figure(figsize=(9, 8))
ax  = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([MB_LON_MIN, MB_LON_MAX, MB_LAT_MIN, MB_LAT_MAX], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.OCEAN.with_scale('10m'), color='#d0e8f5', zorder=0)
ax.add_feature(cfeature.LAND.with_scale('10m'),  color='#e8e0d5', zorder=1)
ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='gray', linewidth=0.8, zorder=2)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', linewidth=1.0, zorder=3)

# Spatial filter to Monterey Bay
def in_mb(lat_arr, lon_arr):
    return (
        (lat_arr >= MB_LAT_MIN) & (lat_arr <= MB_LAT_MAX) &
        (lon_arr >= MB_LON_MIN) & (lon_arr <= MB_LON_MAX)
    )

# PACE pixels
if len(pace_lat) > 0:
    m = in_mb(pace_lat, pace_lon)
    if m.any():
        ax.scatter(pace_lon[m], pace_lat[m],
                   c=COLORS['PACE'], s=4, alpha=0.4,
                   transform=ccrs.PlateCarree(), zorder=4, rasterized=True)

# VIIRS pixels
if len(viirs_lat) > 0:
    m = in_mb(viirs_lat, viirs_lon)
    if m.any():
        ax.scatter(viirs_lon[m], viirs_lat[m],
                   c=COLORS['VIIRS'], s=8, alpha=0.4, marker='s',
                   transform=ccrs.PlateCarree(), zorder=4, rasterized=True)

# SPEXone pixels
if len(spex_lat) > 0:
    m = in_mb(spex_lat, spex_lon)
    if m.any():
        ax.scatter(spex_lon[m], spex_lat[m],
                   c=COLORS['SPEXone'], s=12, alpha=0.7, marker='^',
                   transform=ccrs.PlateCarree(), zorder=4, rasterized=True)

# Aircraft tracks
for df, key in [(df_er2, 'ER2'), (df_cirpas, 'CIRPAS_TO'), (df_airsharp, 'AirSHARP_TO')]:
    if len(df) == 0:
        continue
    m = in_mb(df['Latitude'].values, df['Longitude'].values)
    if m.any():
        ax.plot(df['Longitude'].values[m], df['Latitude'].values[m],
                color=COLORS[key], lw=1.0, alpha=0.85,
                transform=ccrs.PlateCarree(), zorder=5)

gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
gl.top_labels = False
gl.right_labels = False

handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['PACE'],     markersize=7, label=LABELS['PACE']),
    Line2D([0],[0], marker='s', color='w', markerfacecolor=COLORS['VIIRS'],    markersize=7, label=LABELS['VIIRS']),
    Line2D([0],[0], marker='^', color='w', markerfacecolor=COLORS['SPEXone'],  markersize=7, label=LABELS['SPEXone']),
    Line2D([0],[0], color=COLORS['ER2'],         lw=2, label=LABELS['ER2']),
    Line2D([0],[0], color=COLORS['CIRPAS_TO'],   lw=2, label=LABELS['CIRPAS_TO']),
    Line2D([0],[0], color=COLORS['AirSHARP_TO'], lw=2, label=LABELS['AirSHARP_TO']),
]
ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.85)
ax.set_title('Satellite AOD Pixel Coverage and Aircraft Tracks\nMonterey Bay', fontsize=13)

plt.tight_layout()
plt.savefig('coastal_CA_map_monterey_bay.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_map_monterey_bay.png')


# ## Figure 3 – Sample Satellite Pixel Grid
# 
# One representative PACE granule and one VIIRS granule showing the pixel centroid
# grid pattern over the Monterey Bay area, colour-coded by AOT.

# In[ ]:


# Pick a single granule per instrument with the most pixels in the Monterey Bay box
def best_granule_for_mb(files, reader_fn):
    """Return (filepath, lat, lon, aot, ang) for the granule with most MB pixels."""
    best = (0, None, None, None, None, None)
    for fp in files:
        lat, lon, aot, ang = reader_fn(fp)
        m = in_mb(lat, lon) & np.isfinite(aot)
        if m.sum() > best[0]:
            best = (m.sum(), fp, lat, lon, aot, ang)
    return best[1:]

def read_pace_granule(fp):
    with nc.Dataset(fp) as ds:
        lat   = np.array(ds['navigation_data/latitude'][:])
        lon   = np.array(ds['navigation_data/longitude'][:])
        aot   = np.ma.filled(ds['geophysical_data/aot_865'][:], np.nan)
        ang   = np.ma.filled(ds['geophysical_data/angstrom'][:], np.nan)
        flags = np.array(ds['geophysical_data/l2_flags'][:], dtype=np.int32)
    bad = (flags & PACE_BAD_FLAGS).astype(bool)
    aot[bad] = np.nan
    return lat.ravel(), lon.ravel(), aot.ravel(), ang.ravel()

def read_viirs_granule(fp):
    with nc.Dataset(fp) as ds:
        lat = np.array(ds['Latitude'][:])
        lon = np.array(ds['Longitude'][:])
        aot = np.where(np.array(ds['Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'][:]) <= VIIRS_FILL,
                       np.nan, np.array(ds['Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'][:]))
        ang = np.where(np.array(ds['Angstrom_Exponent_Land_Ocean_Best_Estimate'][:]) <= VIIRS_FILL,
                       np.nan, np.array(ds['Angstrom_Exponent_Land_Ocean_Best_Estimate'][:]))
    return lat.ravel(), lon.ravel(), aot.ravel(), ang.ravel()

def read_spexone_granule(fp):
    with nc.Dataset(fp) as ds:
        lat = np.ma.filled(ds['geolocation_data/latitude'][:],  np.nan)
        lon = np.ma.filled(ds['geolocation_data/longitude'][:], np.nan)
        aot = np.where(np.array(ds['geophysical_data/aot550'][:]) <= SPEX_FILL,
                       np.nan, np.array(ds['geophysical_data/aot550'][:]))
        ang = np.where(np.array(ds['geophysical_data/angstrom_440_870'][:]) <= SPEX_FILL,
                       np.nan, np.array(ds['geophysical_data/angstrom_440_870'][:]))
    return lat.ravel(), lon.ravel(), aot.ravel(), ang.ravel()

# All SPEXone files (land + ocean) for the granule search
spex_all_files = sorted(glob.glob(os.path.join(SPEX_DIR, 'PACE_SPEXONE.*.nc')))

print('Selecting best granules for Monterey Bay pixel-grid plot...')
pace_fp_mb,  pace_lat_mb,  pace_lon_mb,  pace_aot_mb,  pace_ang_mb  = best_granule_for_mb(pace_files,     read_pace_granule)
viirs_fp_mb, viirs_lat_mb, viirs_lon_mb, viirs_aot_mb, viirs_ang_mb = best_granule_for_mb(viirs_files,    read_viirs_granule)
spex_fp_mb,  spex_lat_mb,  spex_lon_mb,  spex_aot_mb,  spex_ang_mb  = best_granule_for_mb(spex_all_files, read_spexone_granule)

for lbl, fp, lat_g, lon_g, aot_g in [
    ('PACE',    pace_fp_mb,  pace_lat_mb,  pace_lon_mb,  pace_aot_mb),
    ('VIIRS',   viirs_fp_mb, viirs_lat_mb, viirs_lon_mb, viirs_aot_mb),
    ('SPEXone', spex_fp_mb,  spex_lat_mb,  spex_lon_mb,  spex_aot_mb),
]:
    if fp:
        n = np.sum(in_mb(lat_g, lon_g) & np.isfinite(aot_g)) if lat_g is not None else 0
        print(f'  {lbl}: {os.path.basename(fp)}  ({n:,} valid MB pixels)')


# In[ ]:


# --- Plot pixel grids for PACE, VIIRS, and SPEXone ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6),
                         subplot_kw={'projection': ccrs.PlateCarree()})

aot_norm = mcolors.Normalize(vmin=0.0, vmax=0.4)
aot_cmap = plt.cm.viridis_r

granule_info = [
    (pace_lat_mb,  pace_lon_mb,  pace_aot_mb,  'PACE',
     f'PACE OCI — AOT 865 nm\n{os.path.basename(pace_fp_mb)  if pace_fp_mb  else "(no data)"}',
     4,  'o'),
    (viirs_lat_mb, viirs_lon_mb, viirs_aot_mb, 'VIIRS',
     f'VIIRS DB — AOT 550 nm\n{os.path.basename(viirs_fp_mb) if viirs_fp_mb else "(no data)"}',
     12, 's'),
    (spex_lat_mb,  spex_lon_mb,  spex_aot_mb,  'SPEXone',
     f'SPEXone RemoTAP — AOT 550 nm\n{os.path.basename(spex_fp_mb) if spex_fp_mb else "(no data)"}',
     20, '^'),
]

sc = None
for ax, (lat_g, lon_g, aot_g, key, title, pt_sz, marker) in zip(axes, granule_info):
    ax.set_extent([MB_LON_MIN, MB_LON_MAX, MB_LAT_MIN, MB_LAT_MAX],
                  crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN.with_scale('10m'), color='#d0e8f5', zorder=0)
    ax.add_feature(cfeature.LAND.with_scale('10m'),  color='#e8e0d5', zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', lw=1, zorder=3)
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='gray', lw=0.6, zorder=2)

    if lat_g is not None:
        m = in_mb(lat_g, lon_g)
        if m.any():
            sc = ax.scatter(lon_g[m], lat_g[m], c=aot_g[m],
                            cmap=aot_cmap, norm=aot_norm, s=pt_sz, marker=marker,
                            alpha=0.75, transform=ccrs.PlateCarree(), zorder=4,
                            rasterized=True)

    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title(title, fontsize=9)

if sc is not None:
    cbar = fig.colorbar(sc, ax=axes, orientation='horizontal', fraction=0.03, pad=0.05,
                        shrink=0.55)
    cbar.set_label('AOT', fontsize=11)

fig.suptitle('Sample Satellite Pixel Grid — Monterey Bay', fontsize=13, y=1.01)
plt.savefig('coastal_CA_pixel_grid_monterey_bay.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_pixel_grid_monterey_bay.png')


# ## Figure 4 – Histogram: Distance to Coast
# 
# Distribution of **valid AOD retrievals and aircraft observation points** as a function
# of their signed distance to the nearest coastline.  
# **Positive = over water; negative = over land.**  
# Satellite histograms are normalised by their own total to allow shape comparison across
# the very different pixel counts.  Aircraft histograms are also normalised.

# In[ ]:


cirpas_distwb,cirpas_dist


# In[ ]:


DIST_BINS = np.arange(-10, 10, 0.2)  # km, 20 km bins
bin_ctr   = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])

datasets = [
    ('PACE',        pace_dist,    pace_distwb,   '-',  2.0, COLORS['PACE']),
    ('VIIRS',       viirs_dist,   viirs_distwb,  '-',  2.0, COLORS['VIIRS']),
    ('SPEXone',     spex_dist,    spex_distwb,   '-',   2.0, COLORS['SPEXone']),
    ('ER2',         er2_dist,     er2_distwb,      ':',   0.5, COLORS['ER2']),
    ('CIRPAS_TO',   cirpas_dist,  cirpas_distwb,   '-.',   0.5, COLORS['CIRPAS_TO']),
    ('AirSHARP_TO', airsharp_dist,airsharp_distwb, '--',   0.5, COLORS['AirSHARP_TO']),
]

fig, ax = plt.subplots(figsize=(13, 6))

for key, dist_arr, dist_arrwb, ls, lw, col in datasets:
    if len(dist_arr) == 0:
        continue
    valid_d = dist_arr[np.isfinite(dist_arr)]
    valid_d_wb = dist_arrwb
    if len(valid_d) == 0:
        continue
    counts, _ = np.histogram(valid_d, bins=DIST_BINS)
    countswb, _ = np.histogram(valid_d_wb, bins=DIST_BINS)
    norm_counts = counts  / countswb * 100.0
    if lw<1.0:
        ax.plot(bin_ctr, pd.Series(norm_counts).rolling(window=10).mean(), color=col, ls=ls, lw=lw*4,
            label=f"{LABELS[key]} (smoothed, N={len(valid_d):,})")
        ax.plot(bin_ctr, norm_counts, '.', color=col, markersize=1)#,
            #label=f"{LABELS[key]} (N={len(valid_d):,})")

    else:
        ax.plot(bin_ctr, norm_counts, color=col, ls=ls, lw=lw,
            label=f"{LABELS[key]} (N={len(valid_d):,})")

ax.axvline(0, color='black', lw=1.5, ls=':', label='Coastline (0 km)', zorder=5)
ax.axvspan(-300, 0, alpha=0.06, color='sienna', label='Over land')
ax.axvspan(0, 600, alpha=0.06, color='steelblue', label='Over ocean')

ax.set_xlabel('Distance to Nearest Coast (km)\n← over land | over ocean →', fontsize=12)
ax.set_ylabel('Fraction of Valid Observations (%)', fontsize=12)
ax.set_title('Distribution of Valid AOD Observations\nvs. Distance to California Coastline', fontsize=13)
ax.legend(fontsize=9, framealpha=0.85, ncol=2)
#ax.set_xlim(DIST_BINS[0], DIST_BINS[-1])
ax.set_xlim(-10,10)
ax.set_ylim(-0.2,102)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('coastal_CA_dist_to_coast_histogram.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_dist_to_coast_histogram.png')


# ## Figure 5 – Aerosol Properties vs. Distance to Coast
# 
# Binned statistics (median ± interquartile range) of:
# - **AOT** (PACE 865 nm, VIIRS 550 nm)
# - **Ångström Exponent** (443–865 nm for PACE; 550/865 nm for VIIRS)
# 
# as a function of distance to the California coastline.
# 
# Aircraft platforms do not carry a passive AOD column instrument in this dataset,
# so only satellite retrievals are plotted here.

# In[ ]:


star_aot.shape


# In[ ]:


from scipy.stats import binned_statistic

PROP_BINS = np.arange(-20, 20, 0.2)  # 25 km bins
prop_ctr  = 0.5 * (PROP_BINS[:-1] + PROP_BINS[1:])

def binned_median_iqr(x, y, bins):
    """Return bin centres, medians, 25th and 75th percentiles."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    med, _, _  = binned_statistic(x, y, statistic='median', bins=bins)
    p25, _, _  = binned_statistic(x, y, statistic=lambda v: np.percentile(v, 25), bins=bins)
    p75, _, _  = binned_statistic(x, y, statistic=lambda v: np.percentile(v, 75), bins=bins)
    cnt, _, _  = binned_statistic(x, y, statistic='count',  bins=bins)
    mask = cnt < 10
    med[mask] = np.nan
    p25[mask] = np.nan
    p75[mask] = np.nan
    return 0.5 * (bins[:-1] + bins[1:]), med, p25, p75, cnt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Aerosol Properties vs. Distance to California Coastline\n'
             'Median ± IQR (bins < 10 samples hidden)', fontsize=13)

aot_series = [
    (pace_dist,    pace_aot, 'PACE',       'OCI AOT$_{865}$','--'),
    (viirs_dist,   viirs_aot,'VIIRS',      'AOT$_{550}$',    '-.'),
    (spex_dist,    spex_aot, 'SPEXone',    'AOT$_{550}$',    ':' ),
    (er2_dist,     hsrl_aot, 'ER2',        'HSRL AOT$_{532}$',    ':'),
    (airsharp_dist,star_aot[192,:], 'AirSHARP_TO','4STAR AOT$_{550}$',    '--'),
    
]
ang_series = [
    (pace_dist,  pace_ang,  'PACE',    'OCI AE (443-865 nm)',   '--'),
    (viirs_dist, viirs_ang, 'VIIRS',   'AE (550-865 nm)',   '-.'),
    (spex_dist,  spex_ang,  'SPEXone', 'AE (440-870 nm)',   ':' ),
    (er2_dist,     hsrl_ang, 'ER2',        'HSRL AE (532-355 nm)',    ':'),
    (airsharp_dist,star_ang, 'AirSHARP_TO','4STAR AE (470-865 nm)',    '--'),
]

for ax, series, ylabel, title in [
    (axes[0], aot_series, 'Aerosol Optical Thickness', 'AOT vs. Distance to Coast'),
    (axes[1], ang_series, 'Angstrom Exponent',         'Angstrom Exponent vs. Distance to Coast'),
]:
    for dist_arr, prop_arr, key, lbl_sfx, ls in series:
        if len(dist_arr) == 0 or len(prop_arr) == 0:
            continue
        ctr, med, p25, p75, _ = binned_median_iqr(dist_arr, prop_arr, PROP_BINS)
        valid = np.isfinite(med)
        ax.plot(ctr[valid], med[valid], color=COLORS[key], ls=ls, lw=2,
                label=f'{key} {lbl_sfx} (median)')
        ax.fill_between(ctr[valid], p25[valid], p75[valid],
                        color=COLORS[key], alpha=0.2, label='_nolegend_')

    ax.axvline(0, color='black', lw=1.2, ls=':', zorder=5)
    ax.axvspan(-300, 0, alpha=0.06, color='sienna')
    ax.axvspan(0, 600, alpha=0.06, color='steelblue')
    ax.set_xlabel('Distance to Coast (km)\n← over land | over ocean →', fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, framealpha=0.85)
    ax.set_xlim(PROP_BINS[0], PROP_BINS[-1])
    ax.grid(True, alpha=0.3)

axes[0].set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('coastal_CA_aerosol_vs_dist_to_coast.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_aerosol_vs_dist_to_coast.png')


# ## Figure 6 – Aircraft Track Density vs. Distance to Coast
# 
# Normalised histogram of aircraft observation points as a function of distance to coast,
# showing how each platform samples the coastal gradient.

# In[ ]:


fig, ax = plt.subplots(figsize=(10, 5))

aircraft_datasets = [
    ('ER2',         er2_dist),
    ('CIRPAS_TO',   cirpas_dist),
    ('AirSHARP_TO', airsharp_dist),
]

for key, dist_arr in aircraft_datasets:
    valid_d = dist_arr[np.isfinite(dist_arr)]
    if len(valid_d) == 0:
        continue
    counts, _ = np.histogram(valid_d, bins=DIST_BINS)
    ax.plot(bin_ctr, counts / counts.sum() * 100.0,
            color=COLORS[key], lw=2,
            label=f"{LABELS[key]} (N={len(valid_d):,})")

ax.axvline(0, color='black', lw=1.5, ls=':', label='Coastline', zorder=5)
ax.axvspan(-300, 0, alpha=0.06, color='sienna', label='Over land')
ax.axvspan(0, 600, alpha=0.06, color='steelblue', label='Over ocean')

ax.set_xlabel('Distance to Nearest Coast (km)\n← over land | over ocean →', fontsize=12)
ax.set_ylabel('Fraction of Observations (%)', fontsize=12)
ax.set_title('Aircraft Observation Density vs. Distance to California Coastline', fontsize=13)
ax.legend(fontsize=9, framealpha=0.85)
ax.set_xlim(DIST_BINS[0], DIST_BINS[-1])
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('coastal_CA_aircraft_dist_histogram.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_aircraft_dist_histogram.png')


# # Summary Statistics

# In[ ]:


print('=== Dataset Summary ===')
print()

rows = [
    ('PACE OCI AOT_865',            pace_dist,  pace_aot,  'AOT'),
    ('VIIRS DB AOT_550',            viirs_dist, viirs_aot, 'AOT'),
    ('SPEXone RemoTAP AOT_550',     spex_dist,  spex_aot,  'AOT'),
    ('ER-2 nav',                    er2_dist,   None,      None),
    ('CIRPAS Twin Otter nav',       cirpas_dist,None,      None),
    ('AirSHARP Twin Otter nav',     airsharp_dist,None,    None),
]

for label, dist_arr, aot_arr, aot_lbl in rows:
    if len(dist_arr) == 0:
        print(f'{label}: no data\n')
        continue
    d = dist_arr[np.isfinite(dist_arr)]
    n_ocean = (d > 0).sum()
    n_land  = (d < 0).sum()
    n_total = len(d)
    print(f'{label}:')
    print(f'  Total valid: {n_total:,}  |  ocean: {n_ocean:,} ({100*n_ocean/max(n_total,1):.1f}%)  '
          f'|  land: {n_land:,} ({100*n_land/max(n_total,1):.1f}%)')
    if n_total > 0:
        print(f'  Distance: median={np.median(d):.1f} km, '
              f'IQR=[{np.percentile(d,25):.1f}, {np.percentile(d,75):.1f}] km')
    if aot_arr is not None and len(aot_arr) > 0:
        a = aot_arr[np.isfinite(dist_arr) & np.isfinite(aot_arr)]
        if len(a) > 0:
            print(f'  {aot_lbl}: median={np.median(a):.3f}, '
                  f'mean={np.mean(a):.3f}, '
                  f'IQR=[{np.percentile(a,25):.3f}, {np.percentile(a,75):.3f}]')
    print()


# In[ ]:




