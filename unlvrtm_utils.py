#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     The methods and classes to use with unlvrtm radiative transfer calculations
# 
# Input:
# 
#     none at command line
# 
# Output:
# 
#     input files for unlvrtm
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#    - ...
# 
# Needed Files:
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2022-04-14
#     Modified:
# 

# # Define classes

# In[ ]:


class unl_input
    """
    Classe to define the unl inputs needed to build the namelist.ini files
    
    """
    from datetime import datetime
    import os
    default_control = {
        'freq_as_input':False,
        'spectra_set':[755.0,775.0],
        'interval':1.3,
        'FWHM': 2.6,
        'atmos_type':2,
        'num_layers': -1,
        'surf_alt': -1,
        'surf_press': 1013.0,
        'run_dir': './',
        'data_dir': '/home/sam/unl-vrtm/data/',
        'out_dir': './',
        'format':"""%%% CONTROL MENU %%%    :
Frequency as input?     : {freq_as_input} 
Spectra set (nm or 1/cm): {spectra_set[0]} {spectra_set[1]} 
  - interval(nm or 1/cm): {interval} 
  - FWHM    (nm or 1/cm): {FWHM} 
Atmos. type[1-6,-1 or-2]: {atmos_type} 
  - number of layers    : {num_layers} 
  - surf altitude (m)   : {surf_alt} 
  - surf pressure (hPa) : {surf_press:4.2f} 
Run    directory        : {run_dir} 
Data   directory        : {data_dir} 
Output directory        : {out_dir}  """}
    
    default_geo = {
        'sza': [30],
        'saz': [0],
        'solar_pos_by_time': False,
        'yyymmmdd': '20120101',
        'hhmmss': '103000',
        'lat': 10.0,
        'lon': -20.0,
        'vza': [10 20 30 40 50 60 70 80 90],
        'vaz': [180],
        'format': """%%% GEOMETRY MENU %%%   :
Solar zenith angle      : 30
Solar azimuthal angle   : 0
Solar position by time? : F
  - yyyymmdd  hhmmss    : 20120101 103000
  - latitude  longitude : 10 -20
View zenith angle(s)    : 10 20 30 40 50 60 70 80 90
View azimuthal angle(s) : 180"""}
    
    default_rtm = {'vlidort_on':True,
                   'num_stokes':3,
                   'num_streams':15,
                   'num_legendre':180,
                   'receptor_levels':33,
                   'receptor_direction':1,
                   'pseudo_spherical':True,
                   'SS_correction_scenario':2,
                   'solar_source':True,
                   'flux_factor':1.0,
                   'atmos_emission':False,
                   'surface_emission':False,
                   'surface_T':300,
        'format':"""%%% RTM MENU %%%        :
Turn on RTM (VLIDORT)?  : {vlidort_on}
# of Stokes components  : {num_stokes}
# of dsicrete streams   : {num_streams}
# of P Legendre terms   : {num_legendre}
Receptor level(s)       : {receptor_levels}
Receptor direction      : {receptor_direction}
Do pseudo-spherical?    : {pseudo_spherical}
SS correction scenario  : {SS_correction_scenario}
Do Solar source         : {solar_source}
  - flux factor         : {flux_factor}
Do atmos emission?      : {atmos_emission}
Do surface emission?    : {surface_emission}
  - surface T (K)       : {surface_T}"""}
    
    default_jacobian = {'jacobian_profile':False,
                        'jacobian_column':False,
                        'wrt_gas':False,
                        'wrt_aod':False,
                        'wrt_ssa':False,
                        'wrt_aero_vol':False,
                        'wrt_mode_fraction':False,
                        'wrt_refrac':False,
                        'wrt_shape_factor':False,
                        'wrt_size_dist':False,
                        'wrt_aero_prof':False,
                        'non_varying_vol':False,
                        'non_varying_AOT':False,
                        'jacobian_surface':False,
                        'surf_wrt_brdf':False,
                        'surf_wrt_brdf_param':False,
                        'FD_verification':False,
        'format':"""%%% JACOBIAN MENU %%%   :
Turn on Profile Jacob.? : {jacobian_profile}
Turn on Column Jacob.?  : {jacobian_column}
  - wrt Gas?            : {wrt_gas}
  - wrt AOD?            : {wrt_aod}
  - wrt SSA?            : {wrt_ssa}
  - wrt aerosol volume? : {wrt_aero_vol}
  - wrt mode fraction?  : {wrt_mode_fraction}
  - wrt refractivity?   : {wrt_refrac}
  - wrt shape factor?   : {wrt_shape_factor}
  - wrt size dist?      : {wrt_size_dist}
  - wrt aerosol profile?: {wrt_aero_prof}
Non-varying volume?     : {non_varying_vol}
Non-varying AOT?        : {non_varying_AOT}
Turn on surface Jacob.? : {jacobian_surface}
  - wrt BRDF factor?    : {surf_wrt_brdf}
  - wrt BRDF parameter? : {surf_wrt_brdf_param}
Do FD verication?       : {FD_verification}"""}
    
    default_surface = {'format':"""%%% SURFACE MENU %%%    :
Do Lambertian surface?  : T
  - surface reflecetance: 0.53
Do BRDF/BPDF surface?   : F
  - # of kernels        : 2
  - kernel entries      : Name      Index  Factor  #PARS  PAR(1)  PAR(2)  PAR(3)
           ==> kernel#1 : Lambertian  1    0.02    0      0.0     0.0     0.0
           ==> kernel#2 : Cox-Munk    9    1.0     2      0.0286  1.7796  0.0"""}

    default_aerosol = {'format':"""%%% AEROSOL MENU %%%    :
Turn on aerosol?        : T
Number of aerosol modes : 2
Columnar loading(s)     : 0.3        0.00
  - is AOD?             : T
  - is Vol (um3/um2)?   : F
Mie/T-Mat/User (1/2/-1) : 1           1
Mode #1 Properties      : ....................Mode#1...................
  - refractive index    : 1.44      0.001
  - shape factor        : -1         1.4
  - monodisperse?       : F          1.0
  - size range [um]     : 0.01       10
  - size distribution   : Index      PAR(1)    PAR(2)    PAR(3)
            ==> Entries : 4          0.20      1.6       0
  - vertical range[km]  : 5          6
  - vertical profile    : Index      PAR(1)    PAR(2)
            ==> Entries : 2          8.0       2.0
Mode #2 Properties      : ....................Mode#2...................
  - refractive index    : 1.558     0.001
  - shape factor        : -1         1.0
  - monodisperse?       : F          1.0
  - size range [um]     : 0.05       15
  - size distribution   : Index      PAR(1)    PAR(2)    PAR(3)
            ==> Entries : 5          1.90      0.41      0
  - vertical range[km]  : 2          7.0
  - vertical profile    : Index      PAR(1)    PAR(2)
            ==> Entries : 3          2.0       2.0
Reference aerosol mode  : ..................Mode#REF...................
  - reference AOD?      : F
  - wavelength [nm]     : 550
  - refractive index    : 1.46       0.01
  - monodisperse?       : F          1.0
  - size range [um]     : 0.01       10
  - size distribution   : Index      PAR(1)    PAR(2)    PAR(3)
            ==> Entries : 4          0.20      1.6       0"""}
                       
    default_tracegas = {'format':"""%%% TRACE-GAS MENU %%%  :
Turn on trace gases?    : T
Use pre-calculated data?: T
  - data filename       : oxa_gas/OxA_height_test.unlvrtm.nc
Include HITRAN lines?   : T
Include SAO UV-Vis xsec?: T
Include continuum?      : T
Num. of trace gases     : 22
Trace gas Entries   ==> : TR#   Name    M-Weight   Include?  Scaling
Trace gas #1            :   1   H2O     18.015       T       1.0
Trace gas #2            :   2   CO2     44.010       T       1.0
Trace gas #3            :   3   O3      47.998       T       1.0
Trace gas #4            :   4   N2O     44.010       T       1.0
Trace gas #5            :   5   CO      28.011       T       1.0
Trace gas #6            :   6   CH4     16.043       T       1.0
Trace gas #7            :   7   O2      31.999       T       1.0
Trace gas #8            :   8   NO      30.010       T       1.0
Trace gas #9            :   9   SO2     64.060       T       1.0
Trace gas #10           :  10   NO2     46.010       T       1.0
Trace gas #11           :  11   NH3     17.030       F       1.0
Trace gas #12           :  12   HNO3    63.010       F       1.0
Trace gas #13           :  13   OH      17.000       F       1.0
Trace gas #14           :  14   HF      20.010       F       1.0
Trace gas #15           :  15   HCL     36.460       F       1.0
Trace gas #16           :  16   HBR     80.920       F       1.0
Trace gas #17           :  17   HI      127.91       F       1.0
Trace gas #18           :  18   CLO     51.450       F       1.0
Trace gas #19           :  19   OCS     60.080       F       1.0
Trace gas #20           :  20   H2CO    30.030       F       1.0
Trace gas #21           :  21   HOCL    52.460       F       1.0
Trace gas #22           :  22   N2      28.014       F       1.0"""}
                        
    default_rayleigh = {'format':"""%%% RAYLEIGH MENU %%%   :
Turn on Rayleigh?       : T
Turn on anisotropy?     : T"""
    }
                        
    default_diagnostic = {'format':"""%%% DIAGNOSTIC MENU %%% :
Turn on DIAGNOSTIC?     : T
Output NC file prefix   : OxA_height_5km_PP_sza30
DIAG01: Model inputs    : T
DIAG02: Atmos profiles  : T
DIAG03: Linearized Mie  : T
DIAG04: Optic profiles  : T
DIAG05: Surface property: T
DIAG06: VLIDORT IOP     : T
DIAG07: Radiances       : T
DIAG08: Jacobians       : T"""
    }
    
    default_debug = {'format':"""%%% DEBUG MENU %%%      :
Write VLIDORT inputs?   : T
Turn on screen print?   : F
  - print aerosol calc? : F
  - print Mie calc?     : F
  - print Rayleigh calc?: F
  - print gas calc?     : F
  - print surface calc? : F
  - print RTM calc?     : F"""}
    
    header = """UNL-VRTM v2.1 namelist file {daystr} {user}
************************+******************************************************
"""
    
    split_line = """
------------------------+------------------------------------------------------
"""
    def __init__(self,input_dict,fp='/home/sam/unl-vrtm/unl-vrtm-2.1/run/',user=None):
        'initialize the unl_input class'

        self.fp = fp
        self.control['run_dir'] = fp
        self.fp_namelist = fp+'namelist.ini'
        if user:
            self.user = user 
        else:
            self.user = os.path.expanduser('~').split(os.path.sep)[-1]
        
    def print_namelist_file(self,fp_name=None):
        'Method to print to a namelist file as input (path:fp)'
        if fp_name:
            ff = fp_name
        else:
            ff = self.fp_namelist
        while open(ff,'w') as f:
            f.write(header.format(daystr=datetime.now(),user=self.user))
            f.write(self.control.['format'].format(**self.control))
         
    def __getitem__(self,i):
        'Method to call only the variables in the unl_input class like a dict'
        return self.__dict__.get(i)
    
    def keys(self):
        'Method to wrap the dict call to the unl_input class object'
        return self.__dict__.keys()


# In[10]:


import os


# In[19]:


os.path.expanduser('~').split(os.path.sep)[-1]


# # Plot out data

# In[ ]:




