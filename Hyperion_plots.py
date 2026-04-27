#!/usr/bin/env python
# coding: utf-8

# # Intro
# 
# Name:  
# 
#     Hyperion_plots
# 
# Purpose:  
# 
#     Make some plot for the hyperion isofit reprocessing to view the 'before and after'
# 
# Input:
# 
#     none at command line
# 
# Output:
# 
#     plots
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - numpy
#     - Pyephem
#     - pandas
# 
# Needed Files:
# 
#   - hyperion L0, L1r, L2 files
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, 2026-03-23
#     Modified: 

# # Load the required modules

# In[684]:


import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
import os
import netCDF4 as nc
from datetime import datetime
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import spectral
import spectral.io.envi as envi
import json
import pandas as pd
from sklearn.decomposition import PCA
import cartopy.crs as ccrs


# In[5]:


from hyperion_isofitv3 import test, hyperion, hyperion2isofit


# In[6]:


if __name__ == "__main__":
    import map_utils as mu
    import matplotlib.colors as colors
    from path_utils import getpath


# In[7]:


if __name__ == "__main__":
    fp = getpath('Hyperion')#,path=r"/data2/SBG/Hyperion_data",make_path=True)
    fp


# # Define some basic functions

# In[8]:


def norm_rgb(rgb_raw,rgb_max=None):
    rgb_min = np.nanmin(rgb_raw)
    if rgb_max is None:
        rgb_max = np.nanmax(rgb_raw)
    return (rgb_raw - rgb_min) / (rgb_max - rgb_min)


# In[133]:


def plot_lon_lats(lons,lats,ax=plt.gca(),color='orange'):
    'quick plot of the latitudde and longitude lines'
    lon_lines = ax.contour(lons, levels=10, colors=color, linewidths=0.5, alpha=0.7)
    lat_lines = ax.contour(lats, levels=10, colors=color, linewidths=0.5, alpha=0.7)
    ax.clabel(lon_lines, inline=True, fontsize=8, 
              fmt=lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}")
    ax.clabel(lat_lines, inline=True, fontsize=8, 
              fmt=lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}")


# In[509]:


def select_pixels(rgb_array,verbose=False):
    #optimize to only select in steps of 10
    step = 10
    if type(rgb_array) is type(np.array([0])):
        subsampled = rgb_array[::step, ::step, :]
    else:
        subsampled = rgb_array.isel(**{'Along Track': slice(0, None, step), 'Cross Track': slice(0, None, step)}).values
    h_sub, w_sub, _ = subsampled.shape
    
    # 1. Flatten the image to a list of pixels: (9000*250, 3)
    pixels = subsampled.reshape(-1, 3)
    
    # 2. Calculate Intensity (Brightness) for each pixel
    # Simple average of R, G, B channels
    intensity = np.nanmean(pixels, axis=1)
    
    # --- 5th and 95th Percentile Brightness ---
    dim_val = np.nanpercentile(intensity, 5)
    bright_val = np.nanpercentile(intensity, 95)
    
    # Find the index of the pixel closest to these values
    idx_dim = np.nanargmin(np.abs(intensity - dim_val))
    idx_bright = np.nanargmin(np.abs(intensity - bright_val))
    
    # --- The "Greenest" Pixel ---
    # Defined as max(Green - Red) + (Green - Blue)
    greeness = (pixels[:, 1] - pixels[:, 0]) + (pixels[:, 1] - pixels[:, 2])
    idx_green = np.nanargmax(greeness)
    
    # --- 2 Random Pixels ---
    idx_random = np.random.choice(len(pixels), 2, replace=False)
    
    # 3. Consolidate Indices
    #target_indices = [idx_dim, idx_bright, idx_green, idx_random[0], idx_random[1]]
    target_sub_indices = [idx_dim, idx_bright, idx_green, idx_random[0], idx_random[1]]
    
    # 4. Map back to (y, x) coordinates
    coords = []#[(idx // 250, idx % 250) for idx in target_indices]
    for idx in target_sub_indices:
        y_sub = idx // w_sub
        x_sub = idx % w_sub
        coords.append((y_sub * step, x_sub * step)) # Scale back up
    
    labels = ["5th % Dimmest", "95th % Brightest", "Greenest", "Random 1", "Random 2"]
    if verbose:
        for label, (y, x) in zip(labels, coords):
            print(f"{label}: Pixel at y={y}, x={x}")
    return coords,labels


# In[134]:


def check_if_radcalnet(min_lat,max_lat,min_lon,max_lon,date):
    'function to verify if there is a radcalnet site'
    sites = {'Railroad Valley Playa':{'lat':38.497,'lon':-115.69,'fname':'RVUS_archive.nc','date_start':[2013,4,1],'date_end':[2026,3,22]},  # lat,lon
             'La Crau':{'lat':43.558889,'lon':4.864167,'fname':'LCFR_archive.nc','date_start':[2015,1,7],'date_end':[2026,3,12]},
             'Gobabeb':{'lat':-23.6002,'lon':15.11956,'fname':'GONA_archive.nc','date_start':[2017,7,19],'date_end':[2026,3,16]},
             'Baotou':{'lat':40.8514,'lon':109.6291,'fname':'BTCN_archive.nc','date_start':[2017,4,6],'date_end':[2025,12,10]},
             'Baotou Sand':{'lat':40.86587,'lon':109.6155,'fname':'BSCN_archive.nc','date_start':[2016,6,26],'date_end':[2025,12,10]}
            }
    for s in sites:
        if (sites[s]['lat'] >= min_lat) & (sites[s]['lat'] <= max_lat) & (sites[s]['lon'] >= min_lon) & (sites[s]['lon'] <= max_lon):
            if (datetime(*date) >= datetime(*sites[s]['date_start'])) & (datetime(*date) <= datetime(*sites[s]['date_end'])):
                return s,sites[s]
    return None,None


# In[156]:


def load_rad_cal_net(fname,path,date):
    'load the rad cal net reflectances'
    ff = nc.Dataset(os.path.join(path,fname),'r')
    datetimes = ((ff['Time']['Year'][:] - 1970).astype('datetime64[Y]') + 
                (ff['Time']['DOY'][:] - 1).astype('timedelta64[D]') + 
                ff['Time']['Hour_UTC'][:].astype('timedelta64[h]') + 
                ff['Time']['Minutes_UTC'][:].astype('timedelta64[m]'))
    idx = np.abs(datetimes -  np.datetime64(date)).argmin()
    time_diff = np.abs(datetimes[idx] - np.datetime64(date))
    is_within_limit = time_diff <= np.timedelta64(60, 'm')
    if is_within_limit:
        wvl_rad = ff['Reflectance']['Wavelength'][:]
        rfl_boa = ff['Reflectance']['Reflectance_BOA']['Reflectance_BOA_data'][:,idx]
    else:
        wvl_rad = [None]
        rfl_boa = [None]

    return wvl_rad,rfl_boa
    


# In[463]:


def plot_squares_on_select(coords,colors=['red', 'blue', 'green', 'teal', 'yellow'],square_size=20,ax=plt.gca(),labels=[]):
    'Plotting some squares on the image for the selected pixels'

    for i, (y, x) in enumerate(coords):
        # Rectangle( (x_start, y_start), width, height )
        # We subtract half the size to center the square on the pixel
        rect = patches.Rectangle(
            (x - square_size//2, y - square_size//2), 
            square_size, square_size*4, 
            linewidth=1, 
            edgecolor=colors[i], 
            facecolor='none'
        )
        
        rect.set_path_effects([pe.withStroke(linewidth=2, foreground='black')])
        
        # Add to plot
        ax.add_patch(rect)
    
        # Optional: Add a small label above each square
        if any(labels):
            ax.text(x, y - square_size, labels[i], color=colors[i], 
                    fontsize=9, ha='center')
            


# In[423]:


def print_to_json(dict_in,fullpath):
    # Load existing data or start a new list
    data = json.load(open(fullpath)) if os.path.exists(fullpath) else []
    dict_in['Processed_time'] = datetime.now().isoformat()
    # Append new entry (converting masked arrays to lists)
    data.append({k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in dict_in.items()})
    
    # Overwrite with the updated list
    print(f'..Saving json comparisons file to {fullpath}')
    json.dump(data, open(fullpath, 'w'), indent=4)
    


# In[361]:


def check_if_aeronet(min_lat,max_lat,min_lon,max_lon,date,fp):
    'function to verify if there is a radcalnet site'
    aero = pd.read_csv(os.path.join(fp,'aeronet_locations_extended_v3.txt'), skiprows=1, parse_dates=[6, 7], dayfirst=True)
    ids = np.where((aero['Latitude(decimal_degrees)'] >= min_lat) & (aero['Latitude(decimal_degrees)'] <= max_lat) & (aero['Longitude(decimal_degrees)'] >= min_lon) & (aero['Longitude(decimal_degrees)'] <= max_lon))[0]
    index = []
    if any(ids):
        if type(date) is list: 
            dt = datetime(*date) 
        else:
            dt = date
        ii = (dt >= aero['Data_Start_date(dd-mm-yyyy)'][ids]) & (dt <= aero['Data_End_Date(dd-mm-yyyy)'][ids])
        if any(ii):
            index = ids[ii]
    aeronet_name,info = [],[]
    
    for i in index:
        aeronet_name.append(aero['Name'][i])
        info.append({'lat':aero['Latitude(decimal_degrees)'][i],'lon':aero['Longitude(decimal_degrees)'][i],'date':dt})

    return aeronet_name,info


# In[376]:


def load_aeronet(a_name,date,cache_path='~/aeronet_cache',verbose=True):
    'simple function to load an aeronet aod from internet with cache'
    import requests_cache
    import urllib3
    from io import StringIO
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    url = f'https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={a_name}&year={date.year}&month={date.month}&day={date.day}&year2={date.year}&month2={date.month}&day2={date.day}&AOD20=1&AVG=10&if_no_html=1'
    bad = False
    
    if verbose: print('..getting the aeronet site: '+url)
    session = requests_cache.CachedSession('~/aeronet_cache', expire_after=-1,stale_if_error=True)
    response = session.get(url, verify=False)

    #check if there is a response issue
    try:
        response.raise_for_status()
    except:
        bad = True

    if len(response.text) < 200:
        bad = True
    if "<html" in response.text.lower():
        bad = True

    if bad:
        if verbose: print('..AOD not found')
        return None, None, None

    #now read and pull out nearest aod spectra
    if verbose: print('... reading aeronet')
    df = pd.read_csv(StringIO(response.text), skiprows=5,na_values=['-999.0', '-999'])
                     #parse_dates={'datetime': [1, 2]}, index_col='datetime')
    if verbose: print('.... converting to datetime')
    df['datetime'] = pd.to_datetime(df['Date(dd:mm:yyyy)'] + ' ' + df['Time(hh:mm:ss)'], format='%d:%m:%Y %H:%M:%S')
    df.set_index('datetime', inplace=True)
    df.drop(columns=['Date(dd:mm:yyyy)', 'Time(hh:mm:ss)'], inplace=True)
    aod_cols = [c for c in df.columns if c.startswith('AOD') and c.endswith('nm')]
    sorted_aod = sorted(aod_cols, key=lambda x: int(x.split('_')[1].replace('nm', '')))
    wvls = [int(s.split('_')[1].replace('nm', '')) for s in sorted_aod]
    target_time = pd.to_datetime(date)
    if verbose: print('..... pulling out nearest')
    df = df.groupby(level=0).mean(numeric_only=True).sort_index()
    nearest_idx = df.index.get_indexer([target_time], method='nearest')[0]
    nearest_aod = df.iloc[nearest_idx][sorted_aod].values
    pwv = df.iloc[nearest_idx]['Precipitable_Water(cm)']

    return nearest_aod, wvls, pwv


# In[590]:


def get_pretty_ticks(data):
    if any(data==0):
        data[data==0] = np.nan
    d_min, d_max = np.nanmin(data), np.nanmax(data)
    # Generate 3 values across the range, then round to 1 decimal
    raw_ticks = np.linspace(d_min, d_max, 5)[1:-1] # 5 points, drop ends to stay inside
    return np.round(raw_ticks, 1)


# In[601]:


def plot_pca_line(val, lat_vals,lon_vals, pca, mode='lat', color='white',ax=None):
    # 1. Create a dense set of points along the line
    # If it's a Lat line, Lon varies; if it's a Lon line, Lat varies
    num_points = 20
    if mode == 'lat':
        ind = np.nanargmin(abs(lat_vals-val))
        loni_start = ind-200 if ind>200 else 0
        loni_end = ind+200 if ind+200<len(lon_vals)-1 else len(lon_vals)-1
        lons = np.linspace(np.nanmin(lon_vals[loni_start:loni_end]), np.nanmax(lon_vals[loni_start:loni_end]), num_points)
        lats = np.full(num_points, val)
        direction = 'N' if val>0 else 'S'
            
    else:
        lons = np.full(num_points, val)
        ind = np.nanargmin(abs(lon_vals-val))
        lati_start = ind-2000 if ind>2000 else 0
        lati_end = ind+2000 if ind+2000<len(lat_vals)-1 else len(lat_vals)-1
        lats = np.linspace(np.nanmin(lat_vals[lati_start:lati_end]), np.nanmax(lat_vals[lati_start:lati_end]), num_points)
        #lats = np.linspace(np.nanmin(lat_vals), np.nanmax(lat_vals), num_points)
        direction = 'E' if val>0 else 'W'
    
    # 2. Transform this new line into PCA space
    line_coords = np.column_stack([lons, lats])
    line_pca = pca.transform(line_coords)

    ax.plot(line_pca[:, 0], line_pca[:, 1], color=color, alpha=0.7, lw=0.5, ls='--')
    
    ax.text(line_pca[int(num_points/2), 0], line_pca[int(num_points/2), 1], f"{val:.1f} °{direction}", 
                  color=color, fontsize=8, fontweight='bold')
    #print(line_coords,line_pca,line_pca[int(num_points/2), 0], line_pca[int(num_points/2), 1])
    #import ipdb; ipdb.set_trace()


# In[696]:


def compute_rotated_pole(lats, lons):
    """Compute rotated pole via brute-force angle sweep to minimize tilt."""
    n = len(lats)
    mid = n // 2
    lat0, lon0 = np.radians(lats[mid]), np.radians(lons[mid])
    i1, i2 = max(mid - n // 10, 0), min(mid + n // 10, n - 1)
    lat1, lon1 = np.radians(lats[i1]), np.radians(lons[i1])
    lat2, lon2 = np.radians(lats[i2]), np.radians(lons[i2])
    bearing = np.arctan2(np.sin(lon2-lon1)*np.cos(lat2), np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1))

    best_pole, best_score = None, np.inf
    for az in bearing + np.linspace(-np.pi, np.pi, 360):
        plat = np.arcsin(np.sin(lat0)*np.cos(np.pi/2) + np.cos(lat0)*np.sin(np.pi/2)*np.cos(az))
        plon = lon0 + np.arctan2(np.sin(az)*np.sin(np.pi/2)*np.cos(lat0), np.cos(np.pi/2) - np.sin(lat0)*np.sin(plat))
        trial = ccrs.RotatedPole(pole_longitude=np.degrees(plon), pole_latitude=np.degrees(plat))
        pts = trial.transform_points(ccrs.PlateCarree(), lons.astype(float), lats.astype(float))
        dx, dy = np.abs(pts[-1,0]-pts[0,0]), np.abs(pts[-1,1]-pts[0,1])
        score = dx / (dy + 1e-12)
        if score < best_score: best_pole, best_score = (np.degrees(plat), np.degrees(plon)), score

    # rescale so aspect ratio ~1 in rotated space
    trial = ccrs.RotatedPole(pole_longitude=best_pole[1], pole_latitude=best_pole[0])
    pts = trial.transform_points(ccrs.PlateCarree(), lons.astype(float), lats.astype(float))
    xspan, yspan = pts[:,0].ptp(), pts[:,1].ptp()
    return best_pole, (xspan / yspan if yspan > 0 else 1.0)


# In[685]:


def transform_coords(lons, lats, rotated_crs, geo_crs=ccrs.PlateCarree()):
    """Transform geographic coords into rotated projection space (x, y in meters)."""
    pts = rotated_crs.transform_points(geo_crs, np.asarray(lons, dtype=float), np.asarray(lats, dtype=float))
    return pts[:, 0], pts[:, 1]


# In[686]:


def get_pretty_ticks(data):
    d = np.where(data == 0, np.nan, data)
    return np.round(np.linspace(np.nanmin(d), np.nanmax(d), 5)[1:-1], 1)


# In[706]:


def plot_gridline(val, lats, lons, rotated_crs, mode='lon', color='cyan', ax=None, npts=300):
    """Draw a gridline clipped to data extent with a text label."""
    if mode == 'lon':
        ll = np.linspace(np.nanmin(lats), np.nanmax(lats), npts)
        lx, ly = np.full(npts, val), ll
        label = f'{val:.1f}°{"W" if val < 0 else "E"}'
    else:
        ll = np.linspace(np.nanmin(lons), np.nanmax(lons), npts)
        lx, ly = ll, np.full(npts, val)
        label = f'{val:.1f}°{"S" if val < 0 else "N"}'
    x, y = transform_coords(lx, ly, rotated_crs)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    keep = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    if not keep.any(): return
    ax.plot(x[keep], y[keep], color=color, linewidth=0.6, alpha=0.7)
    # label at midpoint of visible segment, rotated to follow the line
    idx = np.where(keep)[0]
    mid = idx[len(idx)//2]
    dx, dy = x[min(mid+1,len(x)-1)] - x[max(mid-1,0)], y[min(mid+1,len(y)-1)] - y[max(mid-1,0)]
    angle = np.degrees(np.arctan2(dy, dx))
    ax.text(x[mid], y[mid], label, color=color, fontsize=6, ha='center', va='bottom', rotation=angle, rotation_mode='anchor', 
            bbox=dict(facecolor='black', alpha=0.4, pad=0.5, edgecolor='none'))


# In[711]:


def overlay_geography(ax, rot_crs, rx, ry):
    """Overlay coastlines, country borders, and US states. Fails gracefully."""
    import cartopy.feature as cfeature
    from shapely.geometry import MultiLineString, MultiPolygon
    geo = ccrs.PlateCarree()
    inv = geo.transform_points(rot_crs, np.array([rx.min(), rx.max()]), np.array([ry.min(), ry.max()]))
    lon_min, lon_max = inv[:,0].min() - 1, inv[:,0].max() + 1
    lat_min, lat_max = inv[:,1].min() - 1, inv[:,1].max() + 1

    for feat, kwargs in [
        (cfeature.COASTLINE, dict(color='white', linewidth=0.8)),
        (cfeature.BORDERS, dict(color='white', linewidth=0.5, linestyle=':')),
        (cfeature.STATES, dict(color='gray', linewidth=0.3, linestyle='--')),
    ]:
        try:
            for geom in feat.geometries():
                b = geom.bounds
                if b[2] < lon_min or b[0] > lon_max or b[3] < lat_min or b[1] > lat_max: continue
                lines = []
                if hasattr(geom, 'exterior'):
                    lines.append(np.array(geom.exterior.coords))
                elif isinstance(geom, (MultiLineString, MultiPolygon)):
                    for part in geom.geoms:
                        lines.append(np.array(part.exterior.coords) if hasattr(part, 'exterior') else np.array(part.coords))
                else:
                    lines.append(np.array(geom.coords))
                for coords in lines:
                    tx, ty = transform_coords(coords[:,0], coords[:,1], rot_crs)
                    ax.plot(tx, ty, **kwargs)
        except Exception as e:
            print(f"Warning: could not overlay {feat}: {e}")
            continue


# # Define a class for combining all these plotting utils

# In[708]:


class make_hyperion_plots():   
    def __init__(self,scene):
        self.paths = hyperion2isofit.readHyperionL2(scene)
        paths = self.paths
        self.i_rgb = (20,10,4) 
        self.scene = scene
        self.date = datetime.strptime(self.paths['attrs']['Scene_Star'][:17],'%Y:%j:%H:%M:%S')
        if not 'working_directory_move' in paths:
            self.paths['working_directory_move'] = paths['L1dir']
            paths['working_directory_move'] = self.paths['working_directory_move']
        self.paths['radiance'] = self.paths['radiance'].replace(paths['L1dir'],paths['working_directory_move'])
        self.paths['obs'] = self.paths['obs'].replace(paths['L1dir'],paths['working_directory_move'])
        if not 'rfl' in self.paths:
            self.paths['rfl'] = os.path.join(paths['working_directory_move'],'output',scene+'_rfl.hdr')

    def prep_reflectance_calcs(self):
        obs = spectral.open_image(self.paths['obs'])
        self.sza = obs[:,:,4]

        self.paths['lut_full'] = os.path.join(self.paths['working_directory_move'],'lut_full','lut.nc')
        lut_fp = nc.Dataset(self.paths['lut_full'],'r')
        self.solar_irr = lut_fp['solar_irr'][:]

    def rad2rfl(self,rad):
        return np.pi*rad/self.solar_irr*np.cos(self.sza*np.pi/180.0)

    def make_masked_rfl(self,img):
        'function to make a flag for bad reflectance output'
        img_rgb = img.read_bands(self.i_rgb).astype(float)
        bad_pixel_mask = (img_rgb < 0).any(axis=2)
        img_rgb[bad_pixel_mask] = np.nan

        return img_rgb

    def plot_in_out_reflectance(self):
        i_rgb = self.i_rgb
        img = spectral.open_image(self.paths['rfl'])
        img_in = spectral.open_image(self.paths['radiance'])
        
        layout = [['I','O',  'S0'], ['I','O',  'S1'],['I','O',  'S2'],['I','O',  'S3'],['I','O',  'S4']]
        fig, axd = plt.subplot_mosaic(layout, figsize=(9, 6), 
                                      gridspec_kw={'width_ratios': [1, 1, 4]})
        
        axd['O'].set_title('Reflectance output')
        axd['I'].set_title('Reflectance input')
        
        # convert radiance to reflectance
        rfl_in = self.rad2rfl(img_in.load())
        masked_img_rgb = self.make_masked_rfl(img)
        
        view = axd['I'].imshow(norm_rgb(rfl_in[:,:,i_rgb],rgb_max=1.0),aspect=0.2)#,ax=axd['I']
        viewo = axd['O'].imshow(norm_rgb(masked_img_rgb,rgb_max=1.0),aspect=0.2)#,ax=axd['I'])
        axd['I'].set_xlim(1,img.shape[1])
        axd['O'].set_xlim(1,img.shape[1])
        
        plot_lon_lats(self.paths['lonlat'][0],self.paths['lonlat'][1],ax=axd['I'])
        plot_lon_lats(self.paths['lonlat'][0],self.paths['lonlat'][1],ax=axd['O'])
        
        try:
            add_coast(paths['lonlat'][0],paths['lonlat'][1],ax=axd['I'],xlen=img.shape[1],ylen=img.shape[0])
        except:
            print('-- Error adding coast, ignoring **')
        
        coords,labels = select_pixels(masked_img_rgb)
        colors = ['red', 'blue', 'green', 'teal', 'yellow']

        #radcalnet checks
        date = self.date
        site_label,site_info = check_if_radcalnet(self.paths['lonlat'][1].min(),self.paths['lonlat'][1].max(),
                                                  self.paths['lonlat'][0].min(),self.paths['lonlat'][0].max(),[date.year,date.month,date.day])
        dict_compare = {}
        wvl_rad = []
        if site_label:
            wvl_rad,rfl_boa = load_rad_cal_net(site_info['fname'],os.path.dirname(self.paths['working_directory_move']),date)
            if any(wvl_rad):
                labels[-1] = site_label
                dist_sq = (paths['lonlat'][0] - site_info['lon'])**2 + (paths['lonlat'][1] - site_info['lat'])**2
                row, col = np.unravel_index(dist_sq.argmin(), dist_sq.shape)
                coords[-1] = [row,col]
        
        plot_squares_on_select(coords,labels=labels,ax=axd['I'],colors=colors)
        plot_squares_on_select(coords,labels=labels,ax=axd['O'],colors=colors)
        
        #now plot the spectra
        wavelengths = self.paths['w']*1000.0
        self.wl = wavelengths
        for i, (name, (y, x)) in enumerate(zip(labels, coords)):
            ax_s = axd[f'S{i}']         

            try:
                spectrum_o = img.read_pixel(y,x)
                spectrum_i = self.rad2rfl(img_in.read_pixel(y,x))[y,x]
            except:
                import traceback, sys; traceback.print_exception(*sys.exc_info())
                import ipdb; ipdb.set_trace()
            
            ax_s.plot(wavelengths, spectrum_i, color=colors[i], linewidth=2,alpha=0.5,label='input')
            ax_s.plot(wavelengths, spectrum_o, color=colors[i], linewidth=2,label='output')
            if (i==4) & any(wvl_rad):
                ax_s.plot(wvl_rad, rfl_boa, color='k', linewidth=2.5,label='RadCalNet')
                #make a save file for this comparison for a combined comparison afterwards
                dict_compare = {'time':date.isoformat(),'site':site_label,'wvl_rad':wvl_rad,'rfl_rad':np.where(rfl_boa > 9000, None, rfl_boa),
                                'output_rfl':spectrum_o,'input_rfl':spectrum_i,'wvl':wavelengths}
                print_to_json(dict_compare,os.path.join(os.path.dirname(self.paths['working_directory_move']),'radcalnet_compare.json'))
                
            ax_s.legend(frameon=False,prop={'size': 7}, labelspacing=0.1, handletextpad=0.3, handlelength=1.0)
            
            ax_s.set_title(f"Spectrum: {name}", color=colors[i], fontsize=10, loc='left')
            ax_s.set_ylim(-0.1,1.1)
            ax_s.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
            ax_s.grid(True, alpha=0.2)
            axd[f'S{i}'].sharex(axd['S0'])
            if i < 4:
                ax_s.tick_params(labelbottom=False)
            else:
                ax_s.set_xlabel("Wavelength [nm]")
                ax_s.tick_params(labelbottom=True)
                ax_s.set_xticks([500,750,1000,1250,1500,1750,2000,2250])
        
        fig.tight_layout()
        plt.subplots_adjust(top=0.92) 
        start = datetime.strptime(self.paths['attrs']['Scene_Star'][:17],'%Y:%j:%H:%M:%S')
        stop = datetime.strptime(self.paths['attrs']['Scene_Stop'][:17],'%Y:%j:%H:%M:%S')
        mlon,mlat = np.nanmean(self.paths['lonlat'][0]),np.nanmean(self.paths['lonlat'][1])
        lon_fmt = lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}"
        lat_fmt = lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}"

        fig.suptitle(f'{self.scene} - {start:%Y/%m/%d - %H:%M:%S} to {stop:%H:%M:%S} (SZA:{np.nanmean(self.sza):2.0f}°)'+
                     ' at {lon_fmt(mlon)},{lat_fmt(mlat)}')
        return fig     


# In[753]:


def load_AOD(self):
    'check and load the aeronet AOD and PWV'
    fp = os.path.dirname(self.paths['working_directory_move'])
    a_names,a_infos = check_if_aeronet(self.paths['lonlat'][1].min(),self.paths['lonlat'][1].max(),
                                              self.paths['lonlat'][0].min(),self.paths['lonlat'][0].max(),self.date,fp)
    if a_names:
        a_aods,a_wvls,a_pwvs,a_latlons,a_coords = [],[],[],[],[]
        subs_loc = envi.open(paths['subs_loc'])[:,:,:]
        subs_lon = subs_loc[:,0,0]
        subs_lat = subs_loc[:,0,1]
        for i,a in enumerate(a_names):
            aod,wvl,pwv = load_aeronet(a[0],date)
            if aod is None: 
                continue
            a_aods.append(aod)
            a_wvls.append(wvl)
            a_pwvs.append(pwv)
            a_latlons.append([a_infos[i]['lat'],a_infos[i]['lon']])
            dist_sq = (subs_loc[:,0,0] - a_infos[i]['lon'])**2 + (subs_loc[:,0,1] - a_infos[i]['lat'])**2
            a_coords.append(dist_sq.argmin())
        return a_aods,a_wvls,a_pwvs,a_latlons,a_coords,a_names
    else:
        return None, None, None, None, None, None

def plot_AOD(self):
    'check if there is aeronet to compare to in the scene and plot comparison'
    a_aods,a_wlvs,a_pwvs,a_latlons,a_coords, a_names = self.load_AOD()
    if a_aods is None:
        print('** No AOD **')
        return None
        
    substate = envi.open(self.paths['subs_state'])
    aods = substate[:,0,-2]
    pwvs = substate[:,0,-1]

    subs_loc = envi.open(paths['subs_loc'])[:,:,:]
    
    layout = [['A','P',  'AOD'], ['A','P',  'PWV']]
    fig, axd = plt.subplot_mosaic(layout, figsize=(9, 6), gridspec_kw={'width_ratios': [1, 1, 4]})

    axd['A'].set_title('AOD'); axd['P'].set_title('PWV')
    axd['AOD'].set_title('AOD - Aerosol Optical Depth'); axd['PWV'].set_title('PWV - Precipitable Water Vapor [cm]')

    if any(subs_loc[:,0,0]==0): subs_loc[subs_loc[:,0,0]==0,0,0] = np.nan
    if any(subs_loc[:,0,1]==0): subs_loc[subs_loc[:,0,1]==0,0,1] = np.nan

    lons_raw, lats_raw = subs_loc[:, 0, 0], subs_loc[:, 0, 1]
    mask = ~(np.isnan(lons_raw) | np.isnan(lats_raw))

    (pole_lat, pole_lon), aspect = compute_rotated_pole(lats_raw[mask], lons_raw[mask])
    rot_crs = ccrs.RotatedPole(pole_longitude=pole_lon, pole_latitude=pole_lat)
    rx, ry = transform_coords(lons_raw[mask], lats_raw[mask], rot_crs)

    #coords_aligned, mask, pca = self.realign_coords(subs_loc)
    vmax = 0.2 if aods[mask].max() >= 0.2 else 0.4
    viewa = axd['A'].scatter(rx, ry, 4, aods[mask], vmin=0, vmax=vmax, cmap='plasma')
    viewp = axd['P'].scatter(rx, ry, 4, pwvs[mask], vmin=0, vmax=2)

    # lock axis limits to data extent + small pad before drawing gridlines
    pad_x, pad_y = 0.02 * rx.ptp(), 0.02 * ry.ptp()
    for a in [axd['A'], axd['P']]:
        a.set_xlim(rx.min()-pad_x, rx.max()+pad_x); a.set_ylim(ry.min()-pad_y, ry.max()+pad_y)
        overlay_geography(a, rot_crs, rx, ry)

    axd['P'].set_yticks([]); axd['P'].set_xticks([]); axd['A'].set_yticks([]); axd['A'].set_xticks([])
    axd['P'].set_xlabel('cross track'); axd['A'].set_ylabel('along track')
    
    cbara, cbarp = plt.colorbar(viewa,ax=axd['A'],extend='max'), plt.colorbar(viewp,ax=axd['P'],extend='max')

    #plot some lat/lon contours
    lon_levels, lat_levels = get_pretty_ticks(lons_raw), get_pretty_ticks(lats_raw)
    for lon in lon_levels:
        for a in [axd['A'], axd['P']]: plot_gridline(lon, lats_raw, lons_raw, rot_crs, mode='lon', color='cyan', ax=a)
    for lat in lat_levels:
        for a in [axd['A'], axd['P']]: plot_gridline(lat, lats_raw, lons_raw, rot_crs, mode='lat', color='yellow', ax=a)

    colors = ['red', 'blue', 'green', 'teal', 'yellow']
    #latlon_transformed = pca.transform(np.array([[a[1],a[0]] for a in a_latlons]))
    ax, ay = transform_coords(np.array([ll[1] for ll in a_latlons]), np.array([ll[0] for ll in a_latlons]), rot_crs)
    for i,name in enumerate(a_names):
        axd['A'].plot(ax[i], ay[i], 's', color=colors[i])
        axd['A'].text(ax[i], ay[i], name, color=colors[i])

    #import ipdb; ipdb.set_trace()
    self.make_subplot_AOD_PWV(a_names,a_coords,a_wlvs,a_aods,a_pwvs,aods,pwvs,axd,colors)
        
    fig.tight_layout()
    plt.subplots_adjust(top=0.92) 
    start = datetime.strptime(self.paths['attrs']['Scene_Star'][:17],'%Y:%j:%H:%M:%S')
    stop = datetime.strptime(self.paths['attrs']['Scene_Stop'][:17],'%Y:%j:%H:%M:%S')
    mlon,mlat = np.nanmean(self.paths['lonlat'][0]),np.nanmean(self.paths['lonlat'][1])
    lon_fmt,lat_fmt = lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}", lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}"

    fig.suptitle(f'{self.scene} - {start:%Y/%m/%d - %H:%M:%S} to {stop:%H:%M:%S} (SZA:{np.nanmean(self.sza):2.0f}°) at {lon_fmt(mlon)},{lat_fmt(mlat)}')
    return fig

def make_subplot_AOD_PWV(self,a_names,a_coords,a_wlvs,a_aods,a_pwvs,aods,pwvs,axd,colors):
    'Plot the spectra plots for AOD and the PWV'
    for i, (name, y) in enumerate(zip(a_names, a_coords)):      
        try:
            aod_val = aods[y]
            pwv_val = pwvs[y]
        except:
            import traceback, sys; traceback.print_exception(*sys.exc_info())
        afl = np.isfinite(a_aods[i])
        axd['AOD'].plot(np.array(a_wlvs[i])[afl], np.array(a_aods[i])[afl],'x-', color=colors[i], linewidth=2,alpha=0.5,label='AERONET: '+a_names[i])
        axd['AOD'].plot([550], aod_val[0],'o', color=colors[i],label='Isofit')
        axd['PWV'].scatter(a_pwvs[i], pwv_val[0],color=colors[i], alpha=0.5,label='AERONET: '+a_names[i])
        
        dict_compare = {'time':self.date.isoformat(),'site':a_names[i],'wvl':a_wlvs[i],'aeronet_AOD':a_aods[i],
                            'isofit_AOD':aod_val,'isofit_wvl':550,'aeronet_pwv':a_pwvs[i],'isofit_pwv':pwv_val}
        print_to_json(dict_compare,os.path.join(os.path.dirname(self.paths['working_directory_move']),'aeronet_compare.json'))

    
    axd['AOD'].legend(frameon=False,prop={'size': 7}, labelspacing=0.1, handletextpad=0.3, handlelength=1.0)
    #axd['AOD'].set_title(f"AOD comparisons", color=colors[i], fontsize=10, loc='left')
    if aod_val[0]<0.2:
        axd['AOD'].set_ylim(-0.02,0.2); axd['AOD'].set_yticks([0.0,0.05,0.1,0.15,0.2])
    elif aod_val[0]<0.6:
        axd['AOD'].set_ylim(-0.02,0.6); axd['AOD'].set_yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6])
    else:
        axd['AOD'].set_ylim(-0.02,1.0); axd['AOD'].set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
    axd['AOD'].grid(True, alpha=0.2)
    axd['AOD'].set_xlabel("Wavelength [nm]");  axd['AOD'].tick_params(labelbottom=True)
    axd['AOD'].set_xlim(340,1050); axd['AOD'].set_xticks([350,450,500,650,750,850,1000])


    axd['PWV'].plot([0,1,2],[0,1,2],':k')
    axd['PWV'].grid(True,alpha=0.2)
    axd['PWV'].set_xlabel("AERONET PWV [cm]"); axd['PWV'].set_ylabel("Isofit PWV [cm]")
    axd['PWV'].legend(frameon=False,prop={'size': 7}, labelspacing=0.1, handletextpad=0.3, handlelength=1.0)
    

make_hyperion_plots.load_AOD = load_AOD
make_hyperion_plots.make_subplot_AOD_PWV = make_subplot_AOD_PWV
make_hyperion_plots.plot_AOD = plot_AOD


# # Main processing scripts

# In[185]:


def make_error_figure(e,scene=''):
    'make a simple figure with just the error'
    error_fig, ax_err = plt.subplots(figsize=(4, 2))
    ax_err.axis('off') # Hide the plot box/axes
    
    # Format the error message text
    error_text = f"Error making the figure for scene: {scene}\n\nError Detail:\n{str(e)}"
    
    ax_err.text(0.5, 0.5, error_text, 
                color='red', fontweight='bold', ha='center', va='center',
                wrap=True, fontsize=12, family='monospace')
    return error_fig


# In[347]:


def make_figs_for_one_scene(scene,path=''):
    'Main run for building the figures for one scene'
    try:
        hyp_plot = make_hyperion_plots(scene)
        hyp_plot.prep_reflectance_calcs()
        fig = hyp_plot.plot_in_out_reflectance()
        fig_name = scene+'_Reflectance.png'
        if not path:
            path = hyp_plot.paths['L1Rdir']
        rts = 0
    except Exception as e:
        print(f"***Error: {scene} - figure was not made. {str(e)}")
        fig = make_error_figure(e,scene)
        fig_name = scene+'_Reflectance_ERROR.png'
        if not path:
            path = '.'
        rts = -1
    try:
        fig.savefig(os.path.join(path,fig_name),dpi=600,transparent=True)
    except Exception as e:
        print(f"***Error: {scene} - Saving difficulty: {str(e)}")
        rts = -2
    try:
        plt.close(fig)
    except:
        pass

    # now check for aeronet
    try:
        figa = hyp_plot.plot_AOD()
        figa_name = scene+'_AOD_PWV.png'
        rts = 0
    except Exception as e:
        print(f"***Error: {scene} - PWV AOD figure was not made. {str(e)}")
        figa = make_error_figure(e,scene)
        figa_name = scene+'_AOD_PWV_ERROR.png'
        rts = -1
    try:
        figa.savefig(os.path.join(path,figa_name),dpi=600,transparent=True)
    except Exception as e:
        print(f"***Error: {scene} - Saving difficulty: {str(e)}")
        rts = -2
    try:
        plt.close(figa)
    except:
        pass
    return rts        


# In[187]:


def main_figs_looper(in_paths='.'):
    'looping over the paths to select the scenes to plot'
    scenes = [l for l in os.listdir(in_paths) if os.path.isdir(os.path.join(in_paths,l)) and l.startswith('EO1')]
    n = len(scenes)
    for i,scene in enumerate(scenes):
        print(f' {i} of {n}: on scene {scene}')
        try:
            rts = make_figs_for_one_scene(scene,path=in_paths)
        except Exception as e:
            print(f'*** ERROR with scene: {scene}, error: {str(e)} ***' )
        
    print('!Done!')


# In[203]:


if __name__ == "__main__":
    main_figs_looper(in_paths=fp)


# # Select and load a scene

# ## Testing out the scripts

# In[15]:


if __name__ == "__main__":
    i_rgb = (22,11,5) #near rgb
    #i_rgb = (150,11,5) # not real


#  Load the radiance L2 values and plot em

# In[16]:


if __name__ == "__main__":
    img = spectral.open_image(paths['radiance'])
    view = spectral.imshow(img, bands=i_rgb,figsize=(9,9),aspect=0.1)
    ax = plt.gca() # Get current axes
    
    # 2. Extract Lon and Lat arrays
    lons = paths['lonlat'][0]
    lats = paths['lonlat'][1]
    
    # 3. Overlay the grid lines using contour
    # 'levels' defines how many grid lines you want
    lon_lines = ax.contour(lons, levels=10, colors='orange', linewidths=0.5, alpha=0.7)
    lat_lines = ax.contour(lats, levels=10, colors='orange', linewidths=0.5, alpha=0.7)
    
    # 4. Optional: Add labels to the grid lines
    ax.clabel(lon_lines, inline=True, fontsize=11, 
              fmt=lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}")
    
    # Latitude: Positive = N, Negative = S
    ax.clabel(lat_lines, inline=True, fontsize=11, 
              fmt=lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}")
    plt.title('Reflectance')
    plt.show()


# Load the L1R and plot it

# In[65]:


if __name__ == "__main__":
    view = spectral.imshow(image_L1r.values,figsize=(9,9),aspect=0.1)
    ax = plt.gca() # Get current axes
    
    # 2. Extract Lon and Lat arrays
    lons = paths['lonlat'][0]
    lats = paths['lonlat'][1]
    
    # 3. Overlay the grid lines using contour
    # 'levels' defines how many grid lines you want
    lon_lines = ax.contour(lons, levels=10, colors='orange', linewidths=0.5, alpha=0.7)
    lat_lines = ax.contour(lats, levels=10, colors='orange', linewidths=0.5, alpha=0.7)
    
    # 4. Optional: Add labels to the grid lines
    ax.clabel(lon_lines, inline=True, fontsize=11, 
              fmt=lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}")
    
    # Latitude: Positive = N, Negative = S
    ax.clabel(lat_lines, inline=True, fontsize=11, 
              fmt=lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}")
    plt.title('Reflectance')
    
    coords,labels = select_pixels(image_L1r)
    plot_squares_on_select(coords,labels=labels,ax=ax)
    
    plt.show()


# In[82]:


if __name__ == "__main__":
    layout = [['I', 'S0'], ['I', 'S1'],['I', 'S2'],['I', 'S3'],['I', 'S4']]
    fig, axd = plt.subplot_mosaic(layout, figsize=(9, 6), 
                                  gridspec_kw={'width_ratios': [1, 3]})
    
    # plot the main figure
    
    view = axd['I'].imshow(norm_rgb(image_L1r.values), aspect=0.1)#, ax=axd['I'])
    axd['I'].set_title('L1R Reflectance Scene', fontsize=14)
    
    plot_lon_lats(paths['lonlat'][0],paths['lonlat'][1],ax=axd['I'])
    
    coords,labels = select_pixels(image_L1r)
    colors = ['red', 'blue', 'green', 'teal', 'yellow']
    
    plot_squares_on_select(coords,labels=labels,ax=axd['I'],colors=colors)
    
    #now plot the spectra
    wavelengths = wl.values
    for i, (name, (y, x)) in enumerate(zip(labels, coords)):
        # Select the specific subplot for this pixel (S0, S1, etc.)
        ax_s = axd[f'S{i}']
        spectrum = image.isel(**{'Along Track': y, 'Cross Track': x}).values
        ax_s.plot(wavelengths, spectrum, color=colors[i], linewidth=2)
        
        ax_s.set_title(f"Spectrum: {name}", color=colors[i], fontsize=10, loc='left')
        ax_s.set_ylim(-2,55)
        ax_s.grid(True, alpha=0.2)
        axd[f'S{i}'].sharex(axd['S0'])
        if i < 4:
            #ax_s.set_xticklabels([''])
            ax_s.tick_params(labelbottom=False)
        else:
            ax_s.set_xlabel("Wavelength [nm]")
            ax_s.tick_params(labelbottom=True)
            ax_s.set_xticks([500,750,1000,1250,1500,1750,2000,2250])
    
    plt.tight_layout()
    plt.show()


# In[214]:


if __name__ == "__main__":
    img = spectral.open_image(paths['rfl'])
    img_in = spectral.open_image(paths['radiance'])
    layout = [['I','O',  'S0'], ['I','O',  'S1'],['I','O',  'S2'],['I','O',  'S3'],['I','O',  'S4']]
    fig, axd = plt.subplot_mosaic(layout, figsize=(9, 6), 
                                  gridspec_kw={'width_ratios': [1, 1, 4]})
    
    axd['O'].set_title('Reflectance output')
    axd['I'].set_title('Radiance input')
    
    #view = spectral.imshow(img, bands=i_rgb,aspect=0.1)
    view = axd['I'].imshow(norm_rgb(img_in.read_bands(i_rgb)),aspect=0.2)#,ax=axd['I']
    viewo = axd['O'].imshow(norm_rgb(img.read_bands(i_rgb)),aspect=0.2)#,ax=axd['I'])
    
    #ax = plt.gca() # Get current axes
    plot_lon_lats(paths['lonlat'][0],paths['lonlat'][1],ax=axd['I'])
    plot_lon_lats(paths['lonlat'][0],paths['lonlat'][1],ax=axd['O'])
    
    #plt.title('Reflectance output')
    
    colors = ['red', 'blue', 'green', 'teal', 'yellow']
    
    plot_squares_on_select(coords,labels=labels,ax=axd['I'],colors=colors)
    plot_squares_on_select(coords,labels=labels,ax=axd['O'],colors=colors)
    
    #now plot the spectra
    wavelengths = wl.values
    for i, (name, (y, x)) in enumerate(zip(labels, coords)):
        # Select the specific subplot for this pixel (S0, S1, etc.)
        ax_s = axd[f'S{i}']
       # spectrum = img.isel(**{'Along Track': y, 'Cross Track': x}).values
        spectrum_o = img.read_pixel(y,x)
        spectrum_i = img_in.read_pixel(y,x)
        
        ax_s.plot(wavelengths, spectrum_i, color=colors[i], linewidth=2,alpha=0.5,label='input')
        ax_s.plot(wavelengths, spectrum_o, color=colors[i], linewidth=2,label='output')
        ax_s.legend(frameon=False)
        
        ax_s.set_title(f"Spectrum: {name}", color=colors[i], fontsize=10, loc='left')
        ax_s.set_ylim(-0.2,1.2)
        ax_s.grid(True, alpha=0.2)
        axd[f'S{i}'].sharex(axd['S0'])
        if i < 4:
            #ax_s.set_xticklabels([''])
            ax_s.tick_params(labelbottom=False)
        else:
            ax_s.set_xlabel("Wavelength [nm]")
            ax_s.tick_params(labelbottom=True)
            ax_s.set_xticks([500,750,1000,1250,1500,1750,2000,2250])
    
    fig.tight_layout()
    plt.show()


# In[231]:


if __name__ == "__main__":
    img = spectral.open_image(paths['rfl'])
    img_in = spectral.open_image(paths['radiance'])
    layout = [['I','O',  'S0'], ['I','O',  'S1'],['I','O',  'S2'],['I','O',  'S3'],['I','O',  'S4']]
    fig, axd = plt.subplot_mosaic(layout, figsize=(9, 6), 
                                  gridspec_kw={'width_ratios': [1, 1, 4]})
    
    axd['O'].set_title('Reflectance output')
    axd['I'].set_title('Reflectance input')
    
    # convert radiance to reflectance
    rfl_in = rad2rfl(img_in.load())
    
    view = axd['I'].imshow(norm_rgb(rfl_in[:,:,i_rgb],rgb_max=1.0),aspect=0.2)#,ax=axd['I']
    viewo = axd['O'].imshow(norm_rgb(img.read_bands(i_rgb),rgb_max=1.0),aspect=0.2)#,ax=axd['I'])
    axd['I'].set_xlim(1,img.shape[1])
    axd['O'].set_xlim(1,img.shape[1])
    
    plot_lon_lats(paths['lonlat'][0],paths['lonlat'][1],ax=axd['I'])
    plot_lon_lats(paths['lonlat'][0],paths['lonlat'][1],ax=axd['O'])
    
    add_coast(paths['lonlat'][0],paths['lonlat'][1],ax=axd['I'],xlen=img.shape[1],ylen=img.shape[0])
    
    coords,labels = select_pixels(img.read_bands(i_rgb))
    colors = ['red', 'blue', 'green', 'teal', 'yellow']
    
    plot_squares_on_select(coords,labels=labels,ax=axd['I'],colors=colors)
    plot_squares_on_select(coords,labels=labels,ax=axd['O'],colors=colors)
    
    #now plot the spectra
    wavelengths = wl.values
    for i, (name, (y, x)) in enumerate(zip(labels, coords)):
        # Select the specific subplot for this pixel (S0, S1, etc.)
        ax_s = axd[f'S{i}']
       # spectrum = img.isel(**{'Along Track': y, 'Cross Track': x}).values
        spectrum_o = img.read_pixel(y,x)
        spectrum_i = rad2rfl(img_in.read_pixel(y,x))[y,x]
        
        ax_s.plot(wavelengths, spectrum_i, color=colors[i], linewidth=2,alpha=0.5,label='input')
        ax_s.plot(wavelengths, spectrum_o, color=colors[i], linewidth=2,label='output')
        ax_s.legend(frameon=False,prop={'size': 7}, labelspacing=0.1, handletextpad=0.3, handlelength=1.0)
        
        ax_s.set_title(f"Spectrum: {name}", color=colors[i], fontsize=10, loc='left')
        ax_s.set_ylim(-0.2,1.2)
        ax_s.grid(True, alpha=0.2)
        axd[f'S{i}'].sharex(axd['S0'])
        if i < 4:
            #ax_s.set_xticklabels([''])
            ax_s.tick_params(labelbottom=False)
        else:
            ax_s.set_xlabel("Wavelength [nm]")
            ax_s.tick_params(labelbottom=True)
            ax_s.set_xticks([500,750,1000,1250,1500,1750,2000,2250])
    
    fig.tight_layout()
    plt.subplots_adjust(top=0.92) 
    start = datetime.strptime(paths['attrs']['Scene_Star'][:17],'%Y:%j:%H:%M:%S')
    stop = datetime.strptime(paths['attrs']['Scene_Stop'][:17],'%Y:%j:%H:%M:%S')
    mlon,mlat = np.nanmean(paths['lonlat'][0]),np.nanmean(paths['lonlat'][1])
    lon_fmt = lambda x: f"{abs(x):.2f}{'E' if x >= 0 else 'W'}"
    lat_fmt = lambda x: f"{abs(x):.2f}{'N' if x >= 0 else 'S'}"
    fig.suptitle(f'{scene} - {start:%Y/%m/%d - %H:%M:%S} to {stop:%H:%M:%S} at {lon_fmt(mlon)},{lat_fmt(mlat)}')
    plt.show()


# In[234]:


if __name__ == "__main__":
    fig.savefig(os.path.join(paths['L1dir'],scene+'_Reflectance.png'),dpi=600,transparent=True)


# # Tester area

# In[241]:


if __name__ == "__main__":
    hyp_plot = make_hyperion_plots(scene)
    hyp_plot.prep_reflectance_calcs()
    fig = hyp_plot.plot_in_out_reflectance()


# In[519]:


if __name__ == "__main__":
    scene1 = 'EO1H0400332013312110P1'
    hyp_plot = make_hyperion_plots(scene1)
    hyp_plot.prep_reflectance_calcs()
    fig = hyp_plot.plot_in_out_reflectance()
    print(hyp_plot.wl[hyp_plot.i_rgb[0]],hyp_plot.wl[hyp_plot.i_rgb[1]],hyp_plot.wl[hyp_plot.i_rgb[2]])


# In[755]:


if __name__ == "__main__":
    scene1 = 'EO1H0400332013312110P1'
    hyp_plot = make_hyperion_plots(scene1)
    hyp_plot.prep_reflectance_calcs()
    figa = hyp_plot.plot_AOD()
    figa_name = scene1+'_AOD_PWV.png'


# In[745]:


if __name__ == "__main__":
    scene1 = 'EO1H0890822006354110PX'
    hyp_plot = make_hyperion_plots(scene1)
    hyp_plot.prep_reflectance_calcs()
    figa = hyp_plot.plot_AOD()
    figa_name = scene1+'_AOD_PWV.png'


# # end

# In[ ]:




