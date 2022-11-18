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

# In[124]:


class unl_input:
    """
    Classe to define the unl inputs needed to build the namelist.ini files
    
    """
    from datetime import datetime
    import os
    default_control = {
        'freq_as_input':False,
        'spectra_set':[750.0,785.0],
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
        'sza_arr': [30],        
        'saz_arr': [0],
        'solar_pos_by_time': False,
        'yyyymmdd': '20120101',
        'hhmmss': '103000',
        'lat': 10.0,
        'lon': -20.0,
        'vza_arr': [10, 20, 30, 40, 50, 60, 70, 80, 90],
        'vaz_arr': [180],
        'format': """%%% GEOMETRY MENU %%%   :
Solar zenith angle      : {sza_str}
Solar azimuthal angle   : {saz_str}
Solar position by time? : {solar_pos_by_time}
  - yyyymmdd  hhmmss    : {yyyymmdd} {hhmmss}
  - latitude  longitude : {lat} {lon}
View zenith angle(s)    : {vza_str}
View azimuthal angle(s) : {vaz_str}"""}
    default_geo['sza_str'] = ' '.join(map(str,default_geo['sza_arr']))
    default_geo['saz_str'] = ' '.join(map(str,default_geo['saz_arr']))
    default_geo['vza_str'] = ' '.join(map(str,default_geo['vza_arr']))
    default_geo['vaz_str'] = ' '.join(map(str,default_geo['vaz_arr']))
    
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
    
    default_surface = {'lambertian_surface':True,
                       'surface_reflectance':0.3,
                       'do_BRDF_BPDF':False,
                       'BRDF_num_kernels':2,
                       'BRDF_kernel1_name': 'Lambertian',
                       'BRDF_kernel1_index': 1,
                       'BRDF_kernel1_factor': 0.02,
                       'BRDF_kernel1_parsnum': 0,
                       'BRDF_kernel1_pars': [0.0, 0.0, 0.0],
                       'BRDF_kernel2_name': 'Cox-Munk',
                       'BRDF_kernel2_index': 9,
                       'BRDF_kernel2_factor': 1.00,
                       'BRDF_kernel2_parsnum': 2,
                       'BRDF_kernel2_pars': [0.0286, 1.7786, 0.0],                       
                       'format':"""%%% SURFACE MENU %%%    :
Do Lambertian surface?  : {lambertian_surface}
  - surface reflecetance: {surface_reflectance}
Do BRDF/BPDF surface?   : {do_BRDF_BPDF}
  - # of kernels        : {BRDF_num_kernels}
  - kernel entries      : Name      Index  Factor  #PARS  PAR(1)  PAR(2)  PAR(3)
           ==> kernel#1 : {BRDF_kernel1_name: <12}  {BRDF_kernel1_index: <2.2f}    {BRDF_kernel1_factor: <2.2f}  {BRDF_kernel1_parsnum:2.0f}   {BRDF_kernel1_pars[0]}   {BRDF_kernel1_pars[1]}   {BRDF_kernel1_pars[2]}
           ==> kernel#2 : {BRDF_kernel2_name: <12}  {BRDF_kernel2_index: <2.2f}    {BRDF_kernel2_factor: <2.2f}  {BRDF_kernel2_parsnum:2.0f}   {BRDF_kernel2_pars[0]}   {BRDF_kernel2_pars[1]}   {BRDF_kernel2_pars[2]}"""}

    default_aerosol = {'turn_on_aerosol':True,
                       'column_load_is_AOD':True,
                       'column_load_is_vol':False,
                       'column_load_arr': [0.3, 0.0],
                       'mie_tmat_user_select_arr': [1,1],
                       'REF_is_AOD': False,
                       'REF_AOD_wvl':550.0,
                       'REF_ref_index_real':1.46,
                       'REF_ref_index_imag':0.01,
                       'REF_is_monodisperse':False,
                       'REF_monodisperse_factor':1.0,
                       'REF_start_size_um':0.01,
                       'REF_end_size_um':10,
                       'REF_index_size':4,
                       'REF_size_par1':0.20,
                       'REF_size_par2':1.6,
                       'REF_size_par3':0.0,  
                       'mode_str':"""Mode #{mode_num} Properties      : ....................Mode#{mode_num}...................
  - refractive index    : {ref_real_index}     {ref_imag_index}
  - shape factor        : {shape_factor_real}         {shape_factor_imag}
  - monodisperse?       : {is_monodisperse}          {monodisperse_factor}
  - size range [um]     : {start_size_um}       {end_size_um}
  - size distribution   : Index      PAR(1)    PAR(2)    PAR(3)
            ==> Entries : {size_index}          {size_par1}      {size_par2}       {size_par3}
  - vertical range[km]  : {start_alt_km}          {end_alt_km}
  - vertical profile    : Index      PAR(1)    PAR(2)
            ==> Entries : {alt_index}          {alt_par1}      {alt_par2}
""",
                       'modes':[{'mode_num':1,'ref_real_index':1.44,'ref_imag_index':0.001,'shape_factor_real':-1,'shape_factor_imag':1.4,
                                 'is_monodisperse':False,'monodisperse_factor':1.0,'start_size_um':0.01,'end_size_um':10.0,
                                 'size_index':4,'size_par1':0.20,'size_par2':1.6,'size_par3':0.0,'start_alt_km':5.0,'end_alt_km':6.0,
                                 'alt_index':2,'alt_par1':8.0,'alt_par2':2.0},
                                {'mode_num':2,'ref_real_index':1.558,'ref_imag_index':0.001,'shape_factor_real':-1,'shape_factor_imag':1.0,
                                 'is_monodisperse':False,'monodisperse_factor':1.0,'start_size_um':0.05,'end_size_um':15.0,
                                 'size_index':5,'size_par1':1.90,'size_par2':0.41,'size_par3':0.0,'start_alt_km':2.0,'end_alt_km':7.0,
                                 'alt_index':3,'alt_par1':2.0,'alt_par2':2.0}],
        'format':"""%%% AEROSOL MENU %%%    :
Turn on aerosol?        : {turn_on_aerosol}
Number of aerosol modes : {num_aerosol_modes}
Columnar loading(s)     : {column_load_str}
  - is AOD?             : {column_load_is_AOD}
  - is Vol (um3/um2)?   : {column_load_is_vol}
Mie/T-Mat/User (1/2/-1) : {mie_tmat_user_select_str}
{aero_modes}
Reference aerosol mode  : ..................Mode#REF...................
  - reference AOD?      : {REF_is_AOD}
  - wavelength [nm]     : {REF_AOD_wvl}
  - refractive index    : {REF_ref_index_real}       {REF_ref_index_imag}
  - monodisperse?       : {REF_is_monodisperse}          {REF_monodisperse_factor}
  - size range [um]     : {REF_start_size_um}       {REF_end_size_um}
  - size distribution   : Index      PAR(1)    PAR(2)    PAR(3)
            ==> Entries : {REF_index_size}          {REF_size_par1}      {REF_size_par2}       {REF_size_par3}"""}
    default_aerosol['aero_modes'] = ''
    for m in default_aerosol['modes']:
        default_aerosol['aero_modes'] = default_aerosol['aero_modes']+default_aerosol['mode_str'].format(**m)
    default_aerosol['aero_modes'] = default_aerosol['aero_modes'].rstrip()
    
    default_aerosol['num_aerosol_modes'] = len(default_aerosol['modes'])
    default_aerosol['column_load_str'] = ' '.join(map(str,default_aerosol['column_load_arr']))
    default_aerosol['mie_tmat_user_select_str'] = ' '.join(map(str,default_aerosol['mie_tmat_user_select_arr']))
                       
    default_trace_gas = {'turn_on_gas':True,
                        'use_precalc_gases':False,
                        'precalc_filename':'oxa_gas/OxA_height_test.unlvrtm.nc',
                        'include_HITRAN':True,
                        'include_SAO_UV_Vis':True,
                        'include_continuum':True,
                        'trace_gas_line_fmt':"""Trace gas #{num: <2}           :   {num: <2}   {name: <4}     {M_weight:02.3f}       {include}       {scale: >3}
""",
                        'gases':[{'num':1,'name':'H2O','M_weight':18.015,'include':True,'scale':1.0},
                                 {'num':2,'name':'CO2','M_weight':44.010,'include':True,'scale':1.0},
                                 {'num':3,'name':'O3','M_weight':47.998,'include':True,'scale':1.0},
                                 {'num':4,'name':'N2O','M_weight':44.010,'include':True,'scale':1.0},
                                 {'num':5,'name':'CO','M_weight':28.011,'include':True,'scale':1.0},
                                 {'num':6,'name':'CH4','M_weight':16.043,'include':True,'scale':1.0},
                                 {'num':7,'name':'O2','M_weight':31.999,'include':True,'scale':1.0},
                                 {'num':8,'name':'NO','M_weight':30.010,'include':True,'scale':1.0},
                                 {'num':9,'name':'SO2','M_weight':64.060,'include':True,'scale':1.0},
                                 {'num':10,'name':'NO2','M_weight':46.010,'include':True,'scale':1.0},
                                 {'num':11,'name':'NH3','M_weight':17.030,'include':False,'scale':1.0},
                                 {'num':12,'name':'HNO3','M_weight':63.010,'include':False,'scale':1.0},
                                 {'num':13,'name':'OH','M_weight':17.000,'include':False,'scale':1.0},
                                 {'num':14,'name':'HF','M_weight':20.010,'include':False,'scale':1.0},
                                 {'num':15,'name':'HCL','M_weight':36.460,'include':False,'scale':1.0},
                                 {'num':16,'name':'HBR','M_weight':80.920,'include':False,'scale':1.0},
                                 {'num':17,'name':'HI','M_weight':127.91,'include':False,'scale':1.0},
                                 {'num':18,'name':'CLO','M_weight':51.450,'include':False,'scale':1.0},
                                 {'num':19,'name':'OCS','M_weight':60.080,'include':False,'scale':1.0},
                                 {'num':20,'name':'H2CO','M_weight':30.030,'include':False,'scale':1.0},
                                 {'num':21,'name':'HOCL','M_weight':52.460,'include':False,'scale':1.0},
                                 {'num':22,'name':'N2','M_weight':28.014,'include':False,'scale':1.0}],
        'format':"""%%% TRACE-GAS MENU %%%  :
Turn on trace gases?    : {turn_on_gas}
Use pre-calculated data?: {use_precalc_gases}
  - data filename       : {precalc_filename}
Include HITRAN lines?   : {include_HITRAN}
Include SAO UV-Vis xsec?: {include_SAO_UV_Vis}
Include continuum?      : {include_continuum}
Num. of trace gases     : {num_trace_gas}
Trace gas Entries   ==> : TR#   Name    M-Weight   Include?  Scaling
{trace_gas_strs}"""}
    default_trace_gas['num_trace_gas'] = len(default_trace_gas['gases'])
    default_trace_gas['trace_gas_strs'] = ''
    for m in default_trace_gas['gases']:
        default_trace_gas['trace_gas_strs'] = default_trace_gas['trace_gas_strs']+default_trace_gas['trace_gas_line_fmt'].format(**m)
    default_trace_gas['trace_gas_strs'] = default_trace_gas['trace_gas_strs'].rstrip()
                        
    default_rayleigh = {'turn_on_rayleigh':True,
                        'turn_on_anisotropy':True,
        'format':"""%%% RAYLEIGH MENU %%%   :
Turn on Rayleigh?       : {turn_on_rayleigh}
Turn on anisotropy?     : {turn_on_anisotropy}"""
    }
                        
    default_diagnostic = {'turn_on_diagnostic':True,
                          'output_NC_prefix':'OxA_height_5km_PP_sza30',
                          'model_input':True,
                          'atmos_profile':True,
                          'linearized_mie':True,
                          'optic_profiles':True,
                          'surface_property':True,
                          'VLIDORT_IOP': True,
                          'radiances': True,
                          'jacobians': True,
        'format':"""%%% DIAGNOSTIC MENU %%% :
Turn on DIAGNOSTIC?     : {turn_on_diagnostic}
Output NC file prefix   : {output_NC_prefix}
DIAG01: Model inputs    : {model_input}
DIAG02: Atmos profiles  : {atmos_profile}
DIAG03: Linearized Mie  : {linearized_mie}
DIAG04: Optic profiles  : {optic_profiles}
DIAG05: Surface property: {surface_property}
DIAG06: VLIDORT IOP     : {VLIDORT_IOP}
DIAG07: Radiances       : {radiances}
DIAG08: Jacobians       : {jacobians}"""
    }
    
    default_debug = {'write_VLIDORT_inputs':True,
                     'print_to_screen':False,
                     'print_aerosol_calc':False,
                     'print_mie_calc':False,
                     'print_rayleigh_calc':False,
                     'print_gas_calc':False,
                     'print_surface_calc':False,
                     'print_RTM_calc':False,
                     'format':"""%%% DEBUG MENU %%%      :
Write VLIDORT inputs?   : {write_VLIDORT_inputs}
Turn on screen print?   : {print_to_screen}
  - print aerosol calc? : {print_aerosol_calc}
  - print Mie calc?     : {print_mie_calc}
  - print Rayleigh calc?: {print_rayleigh_calc}
  - print gas calc?     : {print_gas_calc}
  - print surface calc? : {print_surface_calc}
  - print RTM calc?     : {print_RTM_calc}"""}
    
    header = """UNL-VRTM v2.1 namelist file {daystr} {user}
************************+******************************************************
"""
    ender = """END OF FILE             :
***********************+******************************************************"""
    
    split_line = """
------------------------+------------------------------------------------------
"""
    defaults = {'control':default_control,'geo':default_geo,'rtm':default_rtm,'jacobian':default_jacobian,'surface':default_surface,
                'aerosol':default_aerosol,'trace_gas':default_trace_gas,'rayleigh':default_rayleigh,'diagnostic':default_diagnostic,'debug':default_debug}
    
    def __init__(self,input_dict,fp='/home/sam/unl-vrtm/unl-vrtm-2.1/run/',user=None):
        'initialize the unl_input class'
        import os

        self.fp = fp        
        self.fp_namelist = fp+'namelist.ini'
        if user:
            self.user = user 
        else:
            self.user = os.path.expanduser('~').split(os.path.sep)[-1]
            
        for defs in self.defaults:
            print(defs)
            self[defs] = self.defaults[defs]
            if defs in input_dict:
                if 'geo' in defs:
                    if 'sza_arr' in input_dict['geo']: input_dict['geo']['sza_str'] = ' '.join(map(str,input_dict['geo']['sza_arr']))
                    if 'saz_arr' in input_dict['geo']: input_dict['geo']['saz_str'] = ' '.join(map(str,input_dict['geo']['saz_arr']))
                    if 'vza_arr' in input_dict['geo']: input_dict['geo']['vza_str'] = ' '.join(map(str,input_dict['geo']['vza_arr']))
                    if 'vaz_arr' in input_dict['geo']: input_dict['geo']['vaz_str'] = ' '.join(map(str,input_dict['geo']['vaz_arr']))
                if 'aerosol' in defs:
                        if 'modes' in input_dict['aerosol']:
                            input_dict['aerosol']['aero_modes'] = ''
                            for m in input_dict['aerosol']['modes']:
                                input_dict['aerosol']['aero_modes'] = input_dict['aerosol']['aero_modes']+self.default_aerosol['mode_str'].format(**m)
                            input_dict['aerosol']['aero_modes'] = input_dict['aerosol']['aero_modes'].rstrip()
                            input_dict['aerosol']['num_aerosol_modes'] = len(input_dict['aerosol']['modes'])
                        if 'column_load_arr' in input_dict['aerosol']: input_dict['aerosol']['column_load_str'] = ' '.join(map(str,input_dict['aerosol']['column_load_arr']))
                        if 'mie_tmat_user_select_arr' in input_dict['aerosol']: input_dict['aerosol']['mie_tmat_user_select_str'] = ' '.join(map(str,input_dict['aerosol']['mie_tmat_user_select_arr']))
                if 'trace_gas' in defs:
                    if 'gases' in input_dict[defs]:
                        input_dict[defs]['num_trace_gas'] = len(input_dict[defs]['gases'])
                        input_dict[defs]['trace_gas_strs'] = ''
                        for m in input_dict[defs]['gases']:
                            input_dict[defs]['trace_gas_strs'] = input_dict[defs]['trace_gas_strs']+self.default_trace_gas['trace_gas_line_fmt'].format(**m)
                        input_dict[defs]['trace_gas_strs'] = input_dict[defs]['trace_gas_strs'].rstrip()
                for dd in input_dict[defs]:
                    self[defs][dd] = input_dict[defs][dd]
                    
        self.control['run_dir'] = fp
        
    def print_namelist_file(self,fp_name=None):
        'Method to print to a namelist file as input (path:fp)'
        from datetime import datetime
        if fp_name:
            ff = fp_name
        else:
            ff = self.fp_namelist
        with open(ff,'w') as f:
            f.write(self.header.format(daystr=datetime.now(),user=self.user))
            for dicts in [self.control,self.geo,self.rtm,self.jacobian,self.surface,self.aerosol,self.trace_gas,self.rayleigh,self.diagnostic,self.debug]:
                f.write(dicts['format'].format(**dicts))
                f.write(self.split_line)
            f.write(self.ender)
         
    def __getitem__(self,i):
        'Method to call only the variables in the unl_input class like a dict'
        return self.__dict__.get(i)
    
    def __setitem__(self, key, value):
        'Method to set items like a dict'
        setattr(self, key, value)
    
    def keys(self):
        'Method to wrap the dict call to the unl_input class object'
        return self.__dict__.keys()


# # Test

# In[ ]:




