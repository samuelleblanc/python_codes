# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%config InlineBackend.rc = {}
import matplotlib 
matplotlib.rc_file('C:\\Users\\sleblan2\\Research\\python_codes\\file.rc')
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
import scipy.io as sio
from mpl_toolkits.basemap import Basemap

# <codecell>

fp = 'C:\Users\sleblan2\Research\Calipso\proposal\\'

# <headingcell level=2>

# Load CERES data

# <codecell>

ceres = sio.netcdf_file(fp+'CERES\CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Subset_200708-200710.nc','r')

# <codecell>

var = ceres.variables

# <codecell>

var

# <codecell>

ceres.comment

# <codecell>

ceres.title

# <codecell>

ceres.filename

# <codecell>

ceres.dimensions

# <codecell>

var['time'].data

# <codecell>

ceres.variables['lon'].data.max()

# <codecell>

ceres.variables['lon'].data.min()

# <codecell>

ceres.variables['lat'].data.max()

# <codecell>

ceres.variables['cldarea_low_mon'].data.shape

# <codecell>

ctr = plt.contourf(ceres.variables['lon'].data,ceres.variables['lat'].data,ceres.variables['cldarea_low_mon'].data[0,:,:],20)
plt.colorbar(ctr)
plt.title('CERES Low Cloud Area Fraction')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# <codecell>

def create_map(lower_left=[-20,-25],upper_right=[20,10],ax=plt.gca()):
    "Create a basemap for the region of South East Atlantic + Africa"
    m = Basemap(projection='stere',lon_0=(upper_right[0]+lower_left[0]),lat_0=(upper_right[1]+lower_left[1]),
            llcrnrlon=lower_left[0], llcrnrlat=lower_left[1],
            urcrnrlon=upper_right[0], urcrnrlat=upper_right[1],resolution='h',ax=ax)
    m.drawcoastlines()
    #m.fillcontinents(color='#AAAAAA')
    m.drawstates()
    m.drawcountries()
    m.drawmeridians(np.linspace(lower_left[0],upper_right[0],8).astype(int),labels=[0,0,0,1])
    m.drawparallels(np.linspace(lower_left[1],upper_right[1],8).astype(int),labels=[1,0,0,0])
    return m

# <codecell>

aero_darf = var['toa_comp_sw-up_naer_mon'][2,:,:]-var['toa_comp_sw-up_all_mon'][2,:,:]

# <codecell>

aero_darf_sep = var['toa_comp_sw-up_naer_mon'][1,:,:]-var['toa_comp_sw-up_all_mon'][1,:,:]

# <codecell>

cld_frac_sep = var['cldarea_low_mon'][1,:,:]

# <codecell>

cld_frac = var['cldarea_low_mon'][2,:,:]

# <codecell>

import cmaps

# <codecell>

cmaps.cmaps()

# <codecell>

fig,(ax0,ax) = plt.subplots(1,2,figsize=(12,4))
m = create_map(ax=ax)
clon,clat = np.meshgrid(var['lon'].data,var['lat'].data,sparse=False)
x,y = m(clon,clat)
ctr = m.contourf(x,y,aero_darf,40,cmap=plt.cm.RdBu_r)
cbr = plt.colorbar(ctr,ax=ax)
cbr.set_label('DARF [W/m$^{2}$]')
ctrl = m.contour(x,y,cld_frac,8,colors='g')
plt.clabel(ctrl, fontsize=9, inline=1,fmt='%2i \%%')
ax.set_title('CERES TOA Shortwave DARF - October 2007')


ax0.plot(cld_frac.flatten(),aero_darf.flatten(),'o')

ax0.text(160,-11.5,'Low Cloud Fraction',color='g')
ax0.set_ylabel('CERES Shortwave DARF [W/m$^{2}$]')
ax0.set_xlabel('Low Cloud Fraction [\%]')
ax0.set_title('Shortwave DARF as a function of Cloud fraction')



plt.savefig(fp+'plots/CERES_SW_DARF_cloud_fraction_Oct2007.png',dpi=600,transparent=True)

# <codecell>

fig,(ax0,ax) = plt.subplots(1,2,figsize=(12,4))
m = create_map(ax=ax)
clon,clat = np.meshgrid(var['lon'].data,var['lat'].data,sparse=False)
x,y = m(clon,clat)
ctr = m.contourf(x,y,aero_darf_sep,40,cmap=plt.cm.RdBu_r)
cbr = plt.colorbar(ctr,ax=ax)
cbr.set_label('DARF [W/m$^{2}$]')
ctrl = m.contour(x,y,cld_frac_sep,8,colors='g')
plt.clabel(ctrl, fontsize=9, inline=1,fmt='%2i \%%')
ax.set_title('CERES TOA Shortwave DARF - September 2007')

ax0.plot(cld_frac_sep.flatten(),aero_darf_sep.flatten(),'o')

ax0.text(160,-12,'Low Cloud Fraction',color='g')
ax0.set_ylabel('CERES Shortwave DARF [W/m$^{2}$]')
ax0.set_xlabel('Low Cloud Fraction [\%]')
ax0.set_title('Shortwave DARF as a function of Cloud fraction')


plt.savefig(fp+'plots/CERES_SW_DARF_cloud_fraction_Sept2007.png',dpi=600,transparent=True)

# <codecell>

aero_darf_aug = var['toa_comp_sw-up_naer_mon'][0,:,:]-var['toa_comp_sw-up_all_mon'][0,:,:]
cld_frac_aug = var['cldarea_low_mon'][0,:,:]

# <codecell>

fig,(ax0,ax) = plt.subplots(1,2,figsize=(12,4))
m = create_map(ax=ax)
clon,clat = np.meshgrid(var['lon'].data,var['lat'].data,sparse=False)
x,y = m(clon,clat)
ctr = m.contourf(x,y,aero_darf_aug,40,cmap=plt.cm.RdBu_r)
cbr = plt.colorbar(ctr,ax=ax)
cbr.set_label('DARF [W/m$^{2}$]')
ctrl = m.contour(x,y,cld_frac_aug,8,colors='g')
plt.clabel(ctrl, fontsize=9, inline=1,fmt='%2i \%%')
ax.set_title('CERES TOA Shortwave DARF - August 2007')

ax0.plot(cld_frac_aug.flatten(),aero_darf_aug.flatten(),'o')

ax0.text(160,-12,'Low Cloud Fraction',color='g')
ax0.set_ylabel('CERES Shortwave DARF [W/m$^{2}$]')
ax0.set_xlabel('Low Cloud Fraction [\%]')
ax0.set_title('Shortwave DARF as a function of Cloud fraction')


plt.savefig(fp+'plots/CERES_SW_DARF_cloud_fraction_Aug2007.png',dpi=600,transparent=True)

# <codecell>


