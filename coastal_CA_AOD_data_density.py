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
# - `/data2/Proposals/Atmosphere_2025/data/PACE_PAX/*.ict` — PACE-PAX ER2, CIRPAS Twin Otter, AirSHARP Twin Otter MetNav ICARTT files  
# 
# **Output:** Maps and histograms displayed inline and saved as PNG files  
# 
# **Dependencies:** `numpy`, `scipy`, `matplotlib`, `cartopy`, `shapely`, `geopandas`, `netCDF4`, `pandas`  
# 
# **History:**  
# Written: Samuel LeBlanc, Santa Cruz, CA, 2026-09-01

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
matplotlib.rc_file('file.rc')


# # Paths and Configuration

# In[3]:


BASE_DIR    = '/data2/Proposals/Atmosphere_2025/data'
PACE_DIR    = os.path.join(BASE_DIR, 'PACE')
VIIRS_DIR   = os.path.join(BASE_DIR, 'VIIRS')
NAV_DIR     = os.path.join(BASE_DIR, 'PACE_PAX')

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

# Colours and styles per dataset
COLORS = {
    'PACE':       '#0077BB',
    'VIIRS':      '#EE7733',
    'ER2':        '#AA3377',
    'CIRPAS_TO':  '#009988',
    'AirSHARP_TO':'#CC3311',
}
LABELS = {
    'PACE':       'PACE OCI AOT$_{865}$',
    'VIIRS':      'VIIRS DB AOT$_{550}$',
    'ER2':        'ER-2 tracks',
    'CIRPAS_TO':  'CIRPAS Twin Otter tracks',
    'AirSHARP_TO':'AirSHARP Twin Otter tracks',
}


# # ICARTT (NASA Ames 1001) File Reader
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


# # Load Aircraft Navigation Data
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
df_airsharp= load_aircraft_group(os.path.join(NAV_DIR, 'AirSHARP-MetNav-1Hz_AirSHARP-TO_*.ict'), 'AirSHARP Twin Otter')


# # Load PACE OCI L2 Satellite Data
# 
# Variables used:
# - `navigation_data/latitude`, `navigation_data/longitude`
# - `geophysical_data/aot_865` — AOT at 865 nm (masked array, scale already applied by netCDF4)
# - `geophysical_data/l2_flags` — bitmask; apply `PACE_BAD_FLAGS` to exclude poor-quality pixels
# 
# Only ocean pixels outside bad-flag conditions are retained.  
# Spatial subset: California bounding box.

# In[6]:


pace_files = sorted(glob.glob(os.path.join(PACE_DIR, 'PACE_OCI.*.L2.OC_AOP.V3_2.nc')))
print(f'Found {len(pace_files)} PACE OCI granules')

pace_lat_all  = []
pace_lon_all  = []
pace_aot_all  = []
pace_land_all = []   # True if pixel is over land (for sign of distance)
pace_ang_all  = []   # Angstrom exponent 443–865 nm

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

pace_lat  = np.concatenate(pace_lat_all)  if pace_lat_all  else np.array([])
pace_lon  = np.concatenate(pace_lon_all)  if pace_lon_all  else np.array([])
pace_aot  = np.concatenate(pace_aot_all)  if pace_aot_all  else np.array([])
pace_ang  = np.concatenate(pace_ang_all)  if pace_ang_all  else np.array([])
pace_land = np.concatenate(pace_land_all) if pace_land_all else np.array([], dtype=bool)

print(f'PACE valid pixels in CA domain: {len(pace_lat):,}')
print(f'  AOT_865 range: {np.nanmin(pace_aot):.3f} – {np.nanmax(pace_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(pace_ang):.3f} – {np.nanmax(pace_ang):.3f}')


# # Load VIIRS NOAA-20 Deep Blue L2 Satellite Data
# 
# Variables used:
# - `Latitude`, `Longitude`
# - `Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate` — already QA-filtered (QA ≥ 2 land, QA = 3 ocean)
# - `Angstrom_Exponent_Land_Ocean_Best_Estimate`
# - `Aerosol_Optical_Thickness_550_Ocean` — used to determine sign (valid → over water)
# 
# Fill value: −999.  Spatial subset: California bounding box.

# In[7]:


viirs_files = sorted(glob.glob(os.path.join(VIIRS_DIR, 'AERDB_L2_VIIRS_NOAA20.*.nc')))
print(f'Found {len(viirs_files)} VIIRS granules')

viirs_lat_all   = []
viirs_lon_all   = []
viirs_aot_all   = []
viirs_ang_all   = []
viirs_ocean_all = []  # True if pixel is over ocean (sign determination)

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

viirs_lat   = np.concatenate(viirs_lat_all)   if viirs_lat_all   else np.array([])
viirs_lon   = np.concatenate(viirs_lon_all)   if viirs_lon_all   else np.array([])
viirs_aot   = np.concatenate(viirs_aot_all)   if viirs_aot_all   else np.array([])
viirs_ang   = np.concatenate(viirs_ang_all)   if viirs_ang_all   else np.array([])
viirs_ocean = np.concatenate(viirs_ocean_all) if viirs_ocean_all else np.array([], dtype=bool)

print(f'VIIRS valid pixels in CA domain: {len(viirs_lat):,}')
print(f'  AOT_550 range: {np.nanmin(viirs_aot):.3f} – {np.nanmax(viirs_aot):.3f}')
print(f'  Angstrom range: {np.nanmin(viirs_ang):.3f} – {np.nanmax(viirs_ang):.3f}')


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

# In[8]:


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


# In[10]:


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

# In[11]:


print('Computing distance to coast...')

# PACE: use LAND flag bit for sign (consistent with own retrieval surface class)
pace_dist = dist_to_coast_km(pace_lat, pace_lon, sign_from_flag=pace_land)
print(f'  PACE:  {np.sum(np.isfinite(pace_dist)):,} valid, '
      f'range {np.nanmin(pace_dist):.1f} – {np.nanmax(pace_dist):.1f} km')

# VIIRS: ocean pixels → positive, land pixels → negative
viirs_sign_land = ~viirs_ocean  # land if ocean retrieval is not valid
viirs_dist = dist_to_coast_km(viirs_lat, viirs_lon, sign_from_flag=viirs_sign_land)
print(f'  VIIRS: {np.sum(np.isfinite(viirs_dist)):,} valid, '
      f'range {np.nanmin(viirs_dist):.1f} – {np.nanmax(viirs_dist):.1f} km')

# Aircraft: use land-mask raster (no instrument flag available for surface type)
er2_dist     = dist_to_coast_km(df_er2['Latitude'],      df_er2['Longitude'])
cirpas_dist  = dist_to_coast_km(df_cirpas['Latitude'],   df_cirpas['Longitude'])
airsharp_dist= dist_to_coast_km(df_airsharp['Latitude'], df_airsharp['Longitude'])
print(f'  ER-2:        {np.sum(np.isfinite(er2_dist)):,} valid')
print(f'  CIRPAS TO:   {np.sum(np.isfinite(cirpas_dist)):,} valid')
print(f'  AirSHARP TO: {np.sum(np.isfinite(airsharp_dist)):,} valid')

print('Done.')


# # Figure 1 – California Coast Overview Map
# 
# Satellite valid-AOD pixel centroids (subsampled for clarity) and all aircraft tracks
# overlaid on a map of the California coast.

# In[13]:


# Subsample satellite pixels for plotting (every N-th point)
PACE_STRIDE  = max(1, len(pace_lat)  // 50000)
VIIRS_STRIDE = max(1, len(viirs_lat) // 50000)

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
               label=LABELS['PACE'], zorder=4, rasterized=True)

# VIIRS pixels
if len(viirs_lat) > 0:
    ax.scatter(viirs_lon[::VIIRS_STRIDE], viirs_lat[::VIIRS_STRIDE],
               c=COLORS['VIIRS'], s=3, alpha=0.25, marker='s',
               transform=ccrs.PlateCarree(), label=LABELS['VIIRS'], zorder=4, rasterized=True)

# Aircraft tracks
if len(df_er2) > 0:
    ax.plot(df_er2['Longitude'].values, df_er2['Latitude'].values,
            color=COLORS['ER2'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), label=LABELS['ER2'], zorder=5)
if len(df_cirpas) > 0:
    ax.plot(df_cirpas['Longitude'].values, df_cirpas['Latitude'].values,
            color=COLORS['CIRPAS_TO'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), label=LABELS['CIRPAS_TO'], zorder=5)
if len(df_airsharp) > 0:
    ax.plot(df_airsharp['Longitude'].values, df_airsharp['Latitude'].values,
            color=COLORS['AirSHARP_TO'], lw=0.5, alpha=0.7,
            transform=ccrs.PlateCarree(), label=LABELS['AirSHARP_TO'], zorder=5)

# Grid lines
gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
gl.top_labels = False
gl.right_labels = False

# Legend
handles = []
if len(pace_lat)  > 0: handles.append(Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['PACE'],    markersize=6, label=LABELS['PACE']))
if len(viirs_lat) > 0: handles.append(Line2D([0],[0], marker='s', color='w', markerfacecolor=COLORS['VIIRS'],   markersize=6, label=LABELS['VIIRS']))
if len(df_er2)    > 0: handles.append(Line2D([0],[0], color=COLORS['ER2'],        lw=1.5, label=LABELS['ER2']))
if len(df_cirpas) > 0: handles.append(Line2D([0],[0], color=COLORS['CIRPAS_TO'],  lw=1.5, label=LABELS['CIRPAS_TO']))
if len(df_airsharp)>0: handles.append(Line2D([0],[0], color=COLORS['AirSHARP_TO'],lw=1.5, label=LABELS['AirSHARP_TO']))
ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.85)

ax.set_title('Satellite AOD Pixel Coverage and Aircraft Tracks\nCalifornia Coast', fontsize=13)
#plt.tight_layout()
plt.savefig('coastal_CA_map_overview.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_map_overview.png')


# # Figure 2 – Monterey Bay Zoom Map

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

# PACE pixels (no subsampling in the small zoom box)
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

# Aircraft tracks
for df, key in [(df_er2, 'ER2'), (df_cirpas, 'CIRPAS_TO'), (df_airsharp, 'AirSHARP_TO')]:
    if len(df) == 0:
        continue
    m = in_mb(df['Latitude'].values, df['Longitude'].values)
    if m.any():
        ax.plot(df['Longitude'].values[m], df['Latitude'].values[m],
                color=COLORS[key], lw=1.0, alpha=0.85,
                transform=ccrs.PlateCarree(), label=LABELS[key], zorder=5)

gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
gl.top_labels = False
gl.right_labels = False

handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=COLORS['PACE'],     markersize=7, label=LABELS['PACE']),
    Line2D([0],[0], marker='s', color='w', markerfacecolor=COLORS['VIIRS'],    markersize=7, label=LABELS['VIIRS']),
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


# # Figure 3 – Sample Satellite Pixel Grid
# 
# One representative PACE granule and one VIIRS granule showing the pixel centroid
# grid pattern over the Monterey Bay area, colour-coded by AOT.

# In[ ]:


# Pick a single PACE granule nearest to Monterey Bay
def best_granule_for_mb(files, lat_var, lon_var, reader_fn):
    """Return the filepath and data for the granule with most pixels in the MB box."""
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
                       np.nan,
                       np.array(ds['Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate'][:]))
        ang = np.where(np.array(ds['Angstrom_Exponent_Land_Ocean_Best_Estimate'][:]) <= VIIRS_FILL,
                       np.nan,
                       np.array(ds['Angstrom_Exponent_Land_Ocean_Best_Estimate'][:]))
    return lat.ravel(), lon.ravel(), aot.ravel(), ang.ravel()

print('Selecting best granule for Monterey Bay pixel-grid plot...')
pace_fp_mb,   pace_lat_mb,  pace_lon_mb,  pace_aot_mb,  pace_ang_mb  = best_granule_for_mb(pace_files,  None, None, read_pace_granule)
viirs_fp_mb, viirs_lat_mb, viirs_lon_mb, viirs_aot_mb, viirs_ang_mb  = best_granule_for_mb(viirs_files, None, None, read_viirs_granule)

if pace_fp_mb:
    print(f'  PACE granule: {os.path.basename(pace_fp_mb)}')
    print(f'    {np.sum(in_mb(pace_lat_mb, pace_lon_mb) & np.isfinite(pace_aot_mb)):,} valid pixels in MB box')
if viirs_fp_mb:
    print(f'  VIIRS granule: {os.path.basename(viirs_fp_mb)}')
    print(f'    {np.sum(in_mb(viirs_lat_mb, viirs_lon_mb) & np.isfinite(viirs_aot_mb)):,} valid pixels in MB box')


# In[ ]:


# --- Plot the pixel grids ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7),
                         subplot_kw={'projection': ccrs.PlateCarree()})

aot_norm = mcolors.Normalize(vmin=0.0, vmax=0.4)
aot_cmap = plt.cm.viridis_r

for ax, (lat_g, lon_g, aot_g, title, src_lbl) in zip(
    axes,
    [
        (pace_lat_mb,  pace_lon_mb,  pace_aot_mb,
         f'PACE OCI — AOT 865 nm\n{os.path.basename(pace_fp_mb) if pace_fp_mb else "(no data)"}',
         'PACE'),
        (viirs_lat_mb, viirs_lon_mb, viirs_aot_mb,
         f'VIIRS DB — AOT 550 nm\n{os.path.basename(viirs_fp_mb) if viirs_fp_mb else "(no data)"}',
         'VIIRS'),
    ]
):
    ax.set_extent([MB_LON_MIN, MB_LON_MAX, MB_LAT_MIN, MB_LAT_MAX],
                  crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN.with_scale('10m'), color='#d0e8f5', zorder=0)
    ax.add_feature(cfeature.LAND.with_scale('10m'),  color='#e8e0d5', zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', lw=1, zorder=3)
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='gray', lw=0.6, zorder=2)

    if lat_g is not None:
        m = in_mb(lat_g, lon_g)
        sc = ax.scatter(lon_g[m], lat_g[m], c=aot_g[m],
                        cmap=aot_cmap, norm=aot_norm,
                        s=4 if src_lbl == 'PACE' else 12,
                        alpha=0.7, transform=ccrs.PlateCarree(), zorder=4,
                        rasterized=True)

    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.5, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    ax.set_title(title, fontsize=10)

cbar = fig.colorbar(sc, ax=axes, orientation='horizontal', fraction=0.03, pad=0.05,
                    shrink=0.6)
cbar.set_label('AOT', fontsize=11)

fig.suptitle('Sample Satellite Pixel Grid — Monterey Bay', fontsize=13, y=1.01)
plt.savefig('coastal_CA_pixel_grid_monterey_bay.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_pixel_grid_monterey_bay.png')


# # Figure 4 – Histogram: Distance to Coast
# 
# Distribution of **valid AOD retrievals and aircraft observation points** as a function
# of their signed distance to the nearest coastline.  
# **Positive = over water; negative = over land.**  
# Satellite histograms are normalised by their own total to allow shape comparison across
# the very different pixel counts.  Aircraft histograms are also normalised.

# In[ ]:


DIST_BINS = np.arange(-300, 601, 20)  # km, 20 km bins from -300 to 600
bin_ctr   = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])

datasets = [
    ('PACE',        pace_dist,      '--',  1.5, COLORS['PACE']),
    ('VIIRS',       viirs_dist,     '-.',  1.5, COLORS['VIIRS']),
    ('ER2',         er2_dist,       '-',   2.0, COLORS['ER2']),
    ('CIRPAS_TO',   cirpas_dist,    '-',   2.0, COLORS['CIRPAS_TO']),
    ('AirSHARP_TO', airsharp_dist,  '-',   2.0, COLORS['AirSHARP_TO']),
]

fig, ax = plt.subplots(figsize=(12, 6))

for key, dist_arr, ls, lw, col in datasets:
    valid_d = dist_arr[np.isfinite(dist_arr)]
    if len(valid_d) == 0:
        continue
    counts, _ = np.histogram(valid_d, bins=DIST_BINS)
    norm_counts = counts / counts.sum() * 100.0  # percent of total
    ax.plot(bin_ctr, norm_counts, color=col, ls=ls, lw=lw,
            label=f"{LABELS[key]} (N={len(valid_d):,})")

ax.axvline(0, color='black', lw=1.5, ls=':', label='Coastline (0 km)', zorder=5)
ax.axvspan(-300, 0, alpha=0.06, color='sienna', label='Over land')
ax.axvspan(0, 600,  alpha=0.06, color='steelblue', label='Over ocean')

ax.set_xlabel('Distance to Nearest Coast (km)\n← over land | over ocean →', fontsize=12)
ax.set_ylabel('Fraction of Valid Observations (%)', fontsize=12)
ax.set_title('Distribution of Valid AOD Observations\nvs. Distance to California Coastline', fontsize=13)
ax.legend(fontsize=9, framealpha=0.85, ncol=2)
ax.set_xlim(DIST_BINS[0], DIST_BINS[-1])
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('coastal_CA_dist_to_coast_histogram.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_dist_to_coast_histogram.png')


# # Figure 5 – Aerosol Properties vs. Distance to Coast
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


from scipy.stats import binned_statistic

PROP_BINS = np.arange(-300, 601, 25)  # 25 km bins
prop_ctr  = 0.5 * (PROP_BINS[:-1] + PROP_BINS[1:])

def binned_median_iqr(x, y, bins):
    """Return bin centres, medians, 25th and 75th percentiles."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    med, _, _  = binned_statistic(x, y, statistic='median', bins=bins)
    p25, _, _  = binned_statistic(x, y, statistic=lambda v: np.percentile(v, 25), bins=bins)
    p75, _, _  = binned_statistic(x, y, statistic=lambda v: np.percentile(v, 75), bins=bins)
    cnt, _, _  = binned_statistic(x, y, statistic='count',  bins=bins)
    # Mask bins with fewer than 10 samples
    mask = cnt < 10
    med[mask] = np.nan
    p25[mask] = np.nan
    p75[mask] = np.nan
    ctr = 0.5 * (bins[:-1] + bins[1:])
    return ctr, med, p25, p75, cnt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Aerosol Properties vs. Distance to California Coastline\n'
             'Median ± IQR (bins < 10 samples hidden)', fontsize=13)

# --- Panel 1: AOT ---
ax = axes[0]
for dist_arr, aot_arr, key, label_suffix, ls in [
    (pace_dist,  pace_aot,  'PACE',  'AOT$_{865}$', '--'),
    (viirs_dist, viirs_aot, 'VIIRS', 'AOT$_{550}$', '-.')
]:
    if len(dist_arr) == 0 or len(aot_arr) == 0:
        continue
    ctr, med, p25, p75, cnt = binned_median_iqr(dist_arr, aot_arr, PROP_BINS)
    valid = np.isfinite(med)
    ax.plot(ctr[valid], med[valid],
            color=COLORS[key], ls=ls, lw=2,
            label=f'{key} {label_suffix} (median)')
    ax.fill_between(ctr[valid], p25[valid], p75[valid],
                    color=COLORS[key], alpha=0.2, label='_nolegend_')

ax.axvline(0, color='black', lw=1.2, ls=':', zorder=5)
ax.axvspan(-300, 0, alpha=0.06, color='sienna')
ax.axvspan(0, 600, alpha=0.06, color='steelblue')
ax.set_xlabel('Distance to Coast (km)\n← over land | over ocean →', fontsize=11)
ax.set_ylabel('Aerosol Optical Thickness', fontsize=11)
ax.set_title('AOT vs. Distance to Coast', fontsize=11)
ax.legend(fontsize=9, framealpha=0.85)
ax.set_xlim(PROP_BINS[0], PROP_BINS[-1])
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)

# --- Panel 2: Ångström Exponent ---
ax = axes[1]
for dist_arr, ang_arr, key, label_suffix, ls in [
    (pace_dist,  pace_ang,  'PACE',  'AE (443–865 nm)', '--'),
    (viirs_dist, viirs_ang, 'VIIRS', 'AE (550/865 nm)', '-.')
]:
    if len(dist_arr) == 0 or len(ang_arr) == 0:
        continue
    ctr, med, p25, p75, cnt = binned_median_iqr(dist_arr, ang_arr, PROP_BINS)
    valid = np.isfinite(med)
    ax.plot(ctr[valid], med[valid],
            color=COLORS[key], ls=ls, lw=2,
            label=f'{key} {label_suffix} (median)')
    ax.fill_between(ctr[valid], p25[valid], p75[valid],
                    color=COLORS[key], alpha=0.2, label='_nolegend_')

ax.axvline(0, color='black', lw=1.2, ls=':', zorder=5)
ax.axvspan(-300, 0, alpha=0.06, color='sienna')
ax.axvspan(0, 600, alpha=0.06, color='steelblue')
ax.set_xlabel('Distance to Coast (km)\n← over land | over ocean →', fontsize=11)
ax.set_ylabel('Ångström Exponent', fontsize=11)
ax.set_title('Ångström Exponent vs. Distance to Coast', fontsize=11)
ax.legend(fontsize=9, framealpha=0.85)
ax.set_xlim(PROP_BINS[0], PROP_BINS[-1])
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('coastal_CA_aerosol_vs_dist_to_coast.png', dpi=150, bbox_inches='tight')
plt.show()
print('Saved coastal_CA_aerosol_vs_dist_to_coast.png')


# # Figure 6 – Aircraft Track Density vs. Distance to Coast
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
    ('PACE OCI AOT_865',          pace_dist,      pace_aot,  'AOT'),
    ('VIIRS DB AOT_550',          viirs_dist,     viirs_aot, 'AOT'),
    ('ER-2 nav',                  er2_dist,       None,      None),
    ('CIRPAS Twin Otter nav',     cirpas_dist,    None,      None),
    ('AirSHARP Twin Otter nav',   airsharp_dist,  None,      None),
]

for label, dist_arr, aot_arr, aot_lbl in rows:
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
    if aot_arr is not None:
        a = aot_arr[np.isfinite(dist_arr) & np.isfinite(aot_arr)]
        print(f'  {aot_lbl}: median={np.median(a):.3f}, '
              f'mean={np.mean(a):.3f}, '
              f'IQR=[{np.percentile(a,25):.3f}, {np.percentile(a,75):.3f}]')
    print()

