#!/usr/bin/env python
# coding: utf-8

# # Intro
# Name:  
# 
#     LeBlanc_2015_paper_file_restore
# 
# Purpose:  
# 
#     Laod the cloud property retrieval data from the paper, for saving into a shareable format with Kokhanovsky
#   
# Input:
# 
#     none at command line
#   
# Output:
# 
#     ict file save
#   
# Keywords:
# 
#     none
#   
# Dependencies:
# 
#     - Sp_parameters.py : for Sp class definition, and for defining the functions used to build parameters
#     - matplotlib
#     - mpltools
#     - numpy
#     - scipy : for saving and reading
#     - plotting_utils (user defined plotting routines)
#     - hdf5storage
#   
# Needed Files:
# 
#   - file.rc : for consistent creation of look of matplotlib figures
#   - retrieved_kisq_20120525_v9.out  retrieved_kisq_20120806_v9.out  retrieved_kisq_20130110_v9.out
#   
#  Modification History:
#  
#      Written: by Samuel LeBlanc, Santa Cruz, 2019-11-29

# # Import modules

# In[73]:


import numpy as np
import hdf5storage as hs
import os
import write_utils as wu
import scipy.io as sio
from path_utils import getpath
import matplotlib.pyplot as plt


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


from Sp_parameters import smooth
from load_utils import load_from_json, mat2py_time,toutc
import write_utils as wu


# In[4]:


import Sun_utils as su
from datetime import datetime
import pytz
import map_utils as mu


# In[5]:


name = 'SSFR3'


# In[6]:


vv = 'v9'


# In[7]:


fp = getpath('SSFR3')
fp_dat = getpath('SSFR3')+'data/'


# # Load the files

# In[8]:


days = ['20120525','20120806','20130110']
labels = ['Liquid','Mixed-phase','Ice']


# In[9]:


times = [[15.0,16.0],[22.0,23.0],[17.5,19.5]]


# ## Load the 15 parameters retrievals

# In[10]:


f_liq = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[0],v=vv))
f_mix = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[1],v=vv))
f_ice = sio.idl.readsav(fp_dat+'retrieved_kisq_{d}_{v}.out'.format(d=days[2],v=vv))


# In[11]:


f_liq.keys()


# ## Load the slope retrievals

# In[12]:


sl_liq = sio.idl.readsav(fp_dat+'{d}_cld_parms3_1600nm.out'.format(d=days[0]))
sl_mix = sio.idl.readsav(fp_dat+'{d}_cld_parms3_1600nm.out'.format(d=days[1]))
sl_ice = sio.idl.readsav(fp_dat+'{d}_cld_parms3_ic_1600nm.out'.format(d=days[2]))


# In[13]:


sl_liq.keys()


# ## Load the 2wvl retrievals

# In[14]:


twv_liq = sio.idl.readsav(fp_dat+'{d}_cld_parms3_2wvl_1600nm.out'.format(d=days[0]))
twv_mix = sio.idl.readsav(fp_dat+'{d}_cld_parms3_2wvl_1600nm.out'.format(d=days[1]))
twv_ice = sio.idl.readsav(fp_dat+'{d}_cld_parms3_ic_2wvl_1600nm.out'.format(d=days[2]))


# In[15]:


twv_liq.keys()


# ## Load the spectra

# In[16]:


sp_liq = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[0],v=vv))
sp_mix = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[1],v=vv))
sp_ice = sio.idl.readsav(fp_dat+'{d}_calibspcs.out'.format(d=days[2],v=vv))


# In[17]:


sp_liq.keys()


# In[18]:


sp_liq['zspectra'].shape


# In[19]:


sp_liq['tmhrs'].shape


# In[20]:


sp_liq['status']


# ## Load the Kurudz solar downwelling

# In[23]:


fso = getpath('uvspec_dat')


# In[26]:


ku = np.genfromtxt(fso+'solar_flux/kurudz_1.0nm.dat')  # mW / (m2 nm), with wavelength in nm


# In[27]:


ku


# In[45]:


ku[980-250,0]


# In[29]:


vis_slit = np.genfromtxt(fp+'pro/rtm2/vis_1nm.dat')


# In[30]:


nir_slit = np.genfromtxt(fp+'pro/rtm2/nir_1nm.dat')


# Calculate the downwelling solar irradiance at TOA by convolving the slit functions and the kurudz

# In[36]:


wvl_vis = sp_liq['zenlambda'][0:194]


# In[38]:


wvl_nir = sp_liq['zenlambda'][194:]


# In[43]:


vis_slit[:,1] = vis_slit[:,1]/np.sum(vis_slit[:,1])


# In[44]:


nir_slit[:,1] = nir_slit[:,1]/np.sum(nir_slit[:,1])


# In[64]:


ku_vis = np.array([np.sum([ku[i+7+int(j),1]*v for j,v in vis_slit]) for i,k in enumerate(ku[7:980-250,1])])


# In[65]:


ku_vis


# In[71]:


ku_nir = np.array([np.sum([ku[i+980-250+int(j),1]*v for j,v in nir_slit]) for i,k in enumerate(ku[980-250:-15,1])])


# In[72]:


plt.figure()
plt.plot(ku[:,0],ku[:,1],'-k',label='orig')
plt.plot(ku[7:980-250,0],ku_vis,'-g',label='slit_vis')
plt.plot(ku[980-250:-15,0],ku_nir,'-r',label='slit_nir')
plt.legend()


# In[ ]:





# # Verify the results and write out

# In[160]:


hdict = {'PI':'Samuel LeBlanc',
     'Institution':'NASA Ames Research Center',
     'Instrument':'Solar Spectral Flux Radiometer - 3 (SSFR3)',
     'campaign':'University of Colorado Skywatch Observatory',
     'special_comments':'',
     'PI_contact':'Samuel.leblanc@nasa.gov',
     'platform':'Roof top University of Colorado',
     'location':'University of Colorado, Duane rooftop, Boulder, Colorado, USA, Lat: 40.01 N, Lon: 105.25 W, Alt: 1660 m',
     'instrument_info':'Derived product from SSFR with zenith narrow field of view radiance light collector',
     'data_info':'Using the cloud property retrieval method based on spectral transmitted light measurements described by LeBlanc, Pileskie, Schmidt, and Coddington (2015), AMT, https://doi.org/10.5194/amt-8-1361-2015',
     'uncertainty':'See included variables.',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'N/A',
     'stipulations':'',
     'rev_comments':'R0: Data from LeBlanc et al., 2015 publication, released to A. Kohkanovsky in November, 2019.'
    }
order = ['COD','COD_err_low','COD_err_up','REF','REF_err_low','REF_err_up','Phase','Ki_square']


# ## 15 parameters retrieval

# In[34]:


f_liq.keys()


# In[45]:


plt.figure()
plt.plot(f_liq['tmhrs'],f_liq['tau_rtm'],'-+')
plt.plot(f_liq['tmhrs'],smooth(f_liq['tau_rtm'],4,old=True),'-x')


# In[46]:


plt.figure()
plt.plot(f_liq['ki_rtm'])


# In[48]:


f_liq['ref_err'].shape


# In[77]:


f = [f_liq,f_mix,f_ice]


# In[86]:


def fo(gu,nn):

    gu['tau_rtm'][nn] = np.nan
    gu['ref_rtm'][nn] = np.nan
    gu['tau_err'][0,nn] = np.nan
    gu['ref_err'][0,nn] = np.nan
    gu['tau_err'][1,nn] = np.nan
    gu['ref_err'][1,nn] = np.nan
    gu['ki_rtm'][nn] = np.nan
    gu['wp_rtm'][nn] = np.nan
    return gu


# In[91]:


rtr = []
for i,g in enumerate(f):
    ff = {}
    for k in g.keys(): ff[k] = np.array(g[k])
    ff['tau_rtm'] = smooth(ff['tau_rtm'],4,old=True)
    ff['ref_rtm'] = smooth(ff['ref_rtm'],4,old=True)
    tr = ff['tau_rtm']==0
    if any(tr): ff = fo(ff,tr)
    trm = ff['tau_rtm']>99.0
    if any(trm): ff = fo(ff,trm)
    rr = ff['ref_rtm']==0
    if any(rr): ff = fo(ff,rr)
    rrm = ff['ref_rtm']>99.0
    if any(rrm): ff = fo(ff,rrm)
    kr = ff['ki_rtm']>0.69
    if any(kr): ff = fo(ff,kr)
    rem = ff['ref_err'][0,:]>3.0
    if any(rem): ff['ref_err'][0,rem] = 3.0
    if i==0:
        wr = ff['wp_rtm']==1
        if any(wr): ff = fo(ff,wr)
    elif i==2:
        wr = ff['wp_rtm']==0
        if any(wr): ff = fo(ff,wr)
    rtr.append(ff)


# In[93]:


rtr[0].keys()


# ### save liquid cloud

# In[161]:


dict_fliq =  {'Start_UTC':{'data':rtr[0]['tmhrs']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtr[0]['tau_rtm'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtr[0]['ref_rtm'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err_low':{'data':rtr[0]['tau_err'][0,:],'unit':'None',
                     'long_description':'Lower value of retrieval uncertainty of Cloud Optical Depth'},
      'COD_err_up':{'data':rtr[0]['tau_err'][1,:],'unit':'None',
                    'long_description':'Upper value of retrieval uncertainty of Cloud Optical Depth'},
      'REF_err_low':{'data':rtr[0]['ref_err'][0,:],'unit':'micrometer',
                     'long_description':'Lower value of retrieval uncertainty of Cloud effective radius.'},
      'REF_err_up':{'data':rtr[0]['ref_err'][1,:],'unit':'micrometer',
                    'long_description':'Upper value of retrieval uncertainty of Cloud effective radius.'},
      'Phase':{'data':rtr[0]['wp_rtm'],'unit':'None',
               'long_description':'Thermodynamic phase, 0 for liquid cloud, 1 for ice cloud'},
      'Ki_square':{'data':rtr[0]['ki_rtm'],'unit':'None',
                   'long_description':'Ki square fit parameter. It is the remainder of the ki square fit, values higher than 0.69 are considered to be failed retrievals.'}}


# In[162]:


hdict_fliq = hdict
hdict_fliq['special_comments'] = 'Liquid cloud case from LeBlanc et al., 2015'


# In[163]:


wu.write_ict(hdict_fliq,dict_fliq,filepath=fp,
              data_id='SSFR_15params_CLD',loc_id='Boulder',date=days[0],rev='R0',order=order)    


# ### Save mix phase cloud

# In[164]:


dict_fmix =  {'Start_UTC':{'data':rtr[1]['tmhrs']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtr[1]['tau_rtm'],'unit':'None',
             'long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtr[1]['ref_rtm'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err_low':{'data':rtr[1]['tau_err'][0,:],'unit':'None',
                     'long_description':'Lower value of retrieval uncertainty of Cloud Optical Depth'},
      'COD_err_up':{'data':rtr[1]['tau_err'][1,:],'unit':'None',
                    'long_description':'Upper value of retrieval uncertainty of Cloud Optical Depth'},
      'REF_err_low':{'data':rtr[1]['ref_err'][0,:],'unit':'micrometer',
                     'long_description':'Lower value of retrieval uncertainty of Cloud effective radius.'},
      'REF_err_up':{'data':rtr[1]['ref_err'][1,:],'unit':'micrometer',
                    'long_description':'Upper value of retrieval uncertainty of Cloud effective radius.'},
      'Phase':{'data':rtr[1]['wp_rtm'],'unit':'None',
               'long_description':'Thermodynamic phase, 0 for liquid cloud, 1 for ice cloud'},
      'Ki_square':{'data':rtr[1]['ki_rtm'],'unit':'None',
                   'long_description':'Ki square fit parameter. It is the remainder of the ki square fit, values higher than 0.69 are considered to be failed retrievals.'}}


# In[165]:


hdict_mix = hdict
hdict_mix['special_comments'] = 'Mixed-phase cloud case from LeBlanc et al., 2015'


# In[166]:


wu.write_ict(hdict_mix,dict_fmix,filepath=fp,
              data_id='SSFR_15params_CLD',loc_id='Boulder',date=days[1],rev='R0',order=order) 


# ### Save ice cloud

# In[167]:


dict_fice =  {'Start_UTC':{'data':rtr[2]['tmhrs']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtr[2]['tau_rtm'],'unit':'None',
             'long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtr[2]['ref_rtm'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err_low':{'data':rtr[2]['tau_err'][0,:],'unit':'None',
                     'long_description':'Lower value of retrieval uncertainty of Cloud Optical Depth'},
      'COD_err_up':{'data':rtr[2]['tau_err'][1,:],'unit':'None',
                    'long_description':'Upper value of retrieval uncertainty of Cloud Optical Depth'},
      'REF_err_low':{'data':rtr[2]['ref_err'][0,:],'unit':'micrometer',
                     'long_description':'Lower value of retrieval uncertainty of Cloud effective radius.'},
      'REF_err_up':{'data':rtr[2]['ref_err'][1,:],'unit':'micrometer',
                    'long_description':'Upper value of retrieval uncertainty of Cloud effective radius.'},
      'Phase':{'data':rtr[2]['wp_rtm'],'unit':'None',
               'long_description':'Thermodynamic phase, 0 for liquid cloud, 1 for ice cloud'},
      'Ki_square':{'data':rtr[2]['ki_rtm'],'unit':'None',
                   'long_description':'Ki square fit parameter. It is the remainder of the ki square fit, values higher than 0.69 are considered to be failed retrievals.'}}


# In[168]:


hdict_ice = hdict
hdict_ice['special_comments'] = 'Ice cloud case from LeBlanc et al., 2015'


# In[169]:


wu.write_ict(hdict_ice,dict_fice,filepath=fp,
              data_id='SSFR_15params_CLD',loc_id='Boulder',date=days[2],rev='R0',order=order) 


# ## slope retrieval

# In[114]:


sl_liq


# In[115]:


sl = [sl_liq,sl_mix,sl_ice]


# In[133]:


rsl = []
for i,g in enumerate(sl):
    fsl = {}
    tr = g['tau']!=0
    fsl['COD'] = g['tau'][tr]
    fsl['REF'] = g['ref'][tr]
    fsl['REF_ERR'] = g['eref'][tr]
    fsl['COD_ERR'] = g['etau'][tr]
    fsl['utc'] = g['tmhrs'][tr]
    rsl.append(fsl)


# In[131]:


rsl[2]['utc']


# ### Save liquid cloud

# In[170]:


hdict_sl = {'PI':'Samuel LeBlanc',
     'Institution':'NASA Ames Research Center',
     'Instrument':'Solar Spectral Flux Radiometer - 3 (SSFR3)',
     'campaign':'University of Colorado Skywatch Observatory',
     'special_comments':'',
     'PI_contact':'Samuel.leblanc@nasa.gov',
     'platform':'Roof top University of Colorado',
     'location':'University of Colorado, Duane rooftop, Boulder, Colorado, USA, Lat: 40.01 N, Lon: 105.25 W, Alt: 1660 m',
     'instrument_info':'Derived product from SSFR with zenith narrow field of view radiance light collector',
     'data_info':'Using the cloud property retrieval method based on slope at 1600 nm from transmitted light measurements described by McBride et al., 2011, doi:10.5194/acp-11-7235-2011.',
     'uncertainty':'See included variables.',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'N/A',
     'stipulations':'',
     'rev_comments':'R0: Data from LeBlanc et al., 2015 publication, released to A. Kohkanovsky in December, 2019.'
    }
order_sl = ['COD','COD_err','REF','REF_err']


# In[171]:


dict_sl_liq =  {'Start_UTC':{'data':rsl[0]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rsl[0]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rsl[0]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rsl[0]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rsl[0]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[172]:


hdict_sl_liq = hdict_sl
hdict_sl_liq['special_comments'] = 'Liquid cloud case from LeBlanc et al., 2015'


# In[173]:


wu.write_ict(hdict_sl_liq,dict_sl_liq,filepath=fp,
              data_id='SSFR_slope_CLD',loc_id='Boulder',date=days[0],rev='R0',order=order_sl)    


# ### Save mix-phase cloud

# In[174]:


dict_sl_mix =  {'Start_UTC':{'data':rsl[1]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rsl[1]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rsl[1]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rsl[1]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rsl[1]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[175]:


hdict_sl_mix = hdict_sl
hdict_sl_mix['special_comments'] = 'Mixed-phase cloud case from LeBlanc et al., 2015'


# In[176]:


wu.write_ict(hdict_sl_mix,dict_sl_mix,filepath=fp,
              data_id='SSFR_slope_CLD',loc_id='Boulder',date=days[1],rev='R0',order=order_sl)    


# ### Save ice cloud

# In[177]:


dict_sl_ice =  {'Start_UTC':{'data':rsl[2]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rsl[2]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rsl[2]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rsl[2]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rsl[2]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[178]:


hdict_sl_ice = hdict_sl
hdict_sl_ice['special_comments'] = 'Ice cloud case from LeBlanc et al., 2015'


# In[179]:


wu.write_ict(hdict_sl_ice,dict_sl_ice,filepath=fp,
              data_id='SSFR_slope_CLD',loc_id='Boulder',date=days[2],rev='R0',order=order_sl)    


# ## 2wvl retrieval
# 

# In[143]:


twv = [twv_liq,twv_mix,twv_ice]


# In[144]:


twv_liq.keys()


# In[146]:


rtw = []
for i,g in enumerate(twv):
    ftw = {}
    tr = g['tau']!=0
    ftw['COD'] = g['tau'][tr]
    ftw['REF'] = g['ref'][tr]
    ftw['REF_ERR'] = g['eref'][tr]
    ftw['COD_ERR'] = g['etau'][tr]
    ftw['utc'] = g['tmhrs'][tr]
    rtw.append(ftw)


# In[148]:


rtw[0]['REF']


# In[180]:


hdict_twv = {'PI':'Samuel LeBlanc',
     'Institution':'NASA Ames Research Center',
     'Instrument':'Solar Spectral Flux Radiometer - 3 (SSFR3)',
     'campaign':'University of Colorado Skywatch Observatory',
     'special_comments':'',
     'PI_contact':'Samuel.leblanc@nasa.gov',
     'platform':'Roof top University of Colorado',
     'location':'University of Colorado, Duane rooftop, Boulder, Colorado, USA, Lat: 40.01 N, Lon: 105.25 W, Alt: 1660 m',
     'instrument_info':'Derived product from SSFR with zenith narrow field of view radiance light collector',
     'data_info':'Using the 2 wavelength cloud property retrieval method from transmitted light measurements.',
     'uncertainty':'See included variables.',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'N/A',
     'stipulations':'',
     'rev_comments':'R0: Data from LeBlanc et al., 2015 publication, released to A. Kohkanovsky in December, 2019.'
    }
order_twv = ['COD','COD_err','REF','REF_err']


# ### Save liquid cloud

# In[181]:


dict_tw_liq =  {'Start_UTC':{'data':rtw[0]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtw[0]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtw[0]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rtw[0]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rtw[0]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[182]:


hdict_tw_liq = hdict_twv
hdict_tw_liq['special_comments'] = 'Liquid cloud case from LeBlanc et al., 2015'


# In[183]:


wu.write_ict(hdict_tw_liq,dict_tw_liq,filepath=fp,
              data_id='SSFR_2wvl_CLD',loc_id='Boulder',date=days[0],rev='R0',order=order_twv)  


# ### Save mix phase cloud

# In[184]:


dict_tw_mix =  {'Start_UTC':{'data':rtw[1]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtw[1]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtw[1]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rtw[1]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rtw[1]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[185]:


hdict_tw_mix = hdict_twv
hdict_tw_mix['special_comments'] = 'Mixed-phase cloud case from LeBlanc et al., 2015'


# In[186]:


wu.write_ict(hdict_tw_mix,dict_tw_mix,filepath=fp,
              data_id='SSFR_2wvl_CLD',loc_id='Boulder',date=days[1],rev='R0',order=order_twv) 


# ### Save ice cloud

# In[263]:


dict_tw_ice =  {'Start_UTC':{'data':rtw[2]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'COD':{'data':rtw[2]['COD'],'unit':'None','long_description':'Cloud Optical Depth of overlying cloud'},
      'REF':{'data':rtw[2]['REF'],'unit':'micrometer',
             'long_description':'Cloud drop effective radius for liquid clouds'},
      'COD_err':{'data':rtw[2]['COD_ERR'],'unit':'None',
                     'long_description':'Retrieval uncertainty of Cloud Optical Depth'},
      'REF_err':{'data':rtw[2]['REF_ERR'],'unit':'micrometer',
                    'long_description':'Retrieval uncertainty of Cloud effective radius.'}}


# In[188]:


hdict_tw_ice = hdict_twv
hdict_tw_ice['special_comments'] = 'Ice cloud case from LeBlanc et al., 2015'


# In[189]:


wu.write_ict(hdict_tw_ice,dict_tw_ice,filepath=fp,
              data_id='SSFR_2wvl_CLD',loc_id='Boulder',date=days[2],rev='R0',order=order_twv) 


# ## Spectra

# In[216]:


lat = 40.01
lon = -105.25
alt = 1660.0


# In[193]:


sp_liq.keys()


# In[225]:


sp = [sp_liq,sp_mix,sp_ice]


# In[230]:


ut = [[15.0,16.0],[22.0,23.0],[17.5,19.5]]


# In[257]:


rsp = []
for i,g in enumerate(sp):
    fsp = {}
    dtu = np.array([datetime(int(days[i][0:4]),int(days[i][4:6]),int(days[i][6:8])+(int(u)/24+1),
                             int(u%24),int((u-int(u))*60),int(((u-int(u))*60)-int((u-int(u))*60))*60,
                             tzinfo=pytz.timezone('UTC')) for u in g['tmhrs']])
    sza, azi = mu.get_sza_azi(lat,lon,dtu,alt=alt)
    tp = (g['tmhrs']>=ut[i][0]) & (g['tmhrs']<=ut[i][1]) & (g['sat']==0) & (g['zspectra'][:,50]>0.0003)
    fsp['utc'] = g['tmhrs'][tp]
    fsp['sza'] = np.array(sza)[tp]
    fsp['wvl'] = g['zenlambda']
    fsp['rad'] = g['zspectra'][tp,:]
    rsp.append(fsp)


# In[254]:


hdict_sp = {'PI':'Samuel LeBlanc',
     'Institution':'NASA Ames Research Center',
     'Instrument':'Solar Spectral Flux Radiometer - 3 (SSFR3)',
     'campaign':'University of Colorado Skywatch Observatory',
     'special_comments':'',
     'PI_contact':'Samuel.leblanc@nasa.gov',
     'platform':'Roof top University of Colorado',
     'location':'University of Colorado, Duane rooftop, Boulder, Colorado, USA, Lat: 40.01 N, Lon: 105.25 W, Alt: 1660 m',
     'instrument_info':'SSFR spectral radiances with zenith narrow field of view radiance light collector',
     'data_info':'Radiance calibrated with lab measurement using NIST traceable lamp and spectralon panel. See details in LeBlanc et al., 2015, AMT, DOI:10.5194/amt-8-1361-2015 ',
     'uncertainty':'month to month variability of accuracy at 8% in radiance. spectral and day to day precision at higher than 0.1%',
     'DM_contact':'Samuel LeBlanc, samuel.leblanc@nasa.gov',
     'project_info':'N/A',
     'stipulations':'',
     'rev_comments':'R0: Data from LeBlanc et al., 2015 publication, released to A. Kohkanovsky in December, 2019.'
    }


# In[259]:


rsp[0]['wvl']


# ### Save the liquid cloud

# In[275]:


dict_sp_liq =  {'Start_UTC':{'data':rsp[0]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'SZA':{'data':rsp[0]['sza'],'unit':'degrees','long_description':'Solar Zenith Angle of observations'}}
order_sp = ['SZA']


# In[276]:


for i,n in enumerate(rsp[0]['wvl']):
    dict_sp_liq['rad_{:06.1f}'.format(n)]={'data':rsp[0]['rad'][:,i],'unit':'W/m^2/nm/sr','long_description':'Zenith radiance measured at {} nm'.format(n)}
    order_sp.append('rad_{:06.1f}'.format(n))


# In[277]:


hdict_sp_liq = hdict_sp
hdict_sp_liq['special_comments'] = 'Liquid cloud case from LeBlanc et al., 2015'


# In[278]:


wu.write_ict(hdict_sp_liq,dict_sp_liq,filepath=fp,
              data_id='SSFR_Zenith_radiance',loc_id='Boulder',date=days[0],rev='R0',order=order_sp) 


# ### Save the mix phase cloud

# In[279]:


dict_sp_mix =  {'Start_UTC':{'data':rsp[1]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'SZA':{'data':rsp[1]['sza'],'unit':'degrees','long_description':'Solar Zenith Angle of observations'}}


# In[280]:


for i,n in enumerate(rsp[1]['wvl']):
    dict_sp_mix['rad_{:06.1f}'.format(n)]={'data':rsp[1]['rad'][:,i],'unit':'W/m^2/nm/sr','long_description':'Zenith radiance measured at {} nm'.format(n)}


# In[281]:


hdict_sp_mix = hdict_sp
hdict_sp_mix['special_comments'] = 'Mixed-phase cloud case from LeBlanc et al., 2015'


# In[282]:


wu.write_ict(hdict_sp_mix,dict_sp_mix,filepath=fp,
              data_id='SSFR_Zenith_radiance',loc_id='Boulder',date=days[1],rev='R0',order=order_sp)


# ### Save the ice cloud

# In[287]:


dict_sp_ice =  {'Start_UTC':{'data':rsp[2]['utc']*3600.0,'unit':'seconds from midnight UTC',
                           'long_description':'time keeping, based on UTC midnight'},
      'SZA':{'data':rsp[2]['sza'],'unit':'degrees','long_description':'Solar Zenith Angle of observations'}}
order_spi = ['SZA']


# In[288]:


for i,n in enumerate(rsp[2]['wvl']):
    dict_sp_ice['rad_{:06.1f}'.format(n)]={'data':rsp[2]['rad'][:,i],'unit':'W/m^2/nm/sr','long_description':'Zenith radiance measured at {} nm'.format(n)}
    order_spi.append('rad_{:06.1f}'.format(n))


# In[289]:


hdict_sp_ice = hdict_sp
hdict_sp_ice['special_comments'] = 'Ice cloud case from LeBlanc et al., 2015'


# In[290]:


wu.write_ict(hdict_sp_ice,dict_sp_ice,filepath=fp,
              data_id='SSFR_Zenith_radiance',loc_id='Boulder',date=days[2],rev='R0',order=order_spi)


# # Update for specific format requested by Kokhanovsky

# From: Kokhanovsky Alexander <a.kokhanovsky@vitrocisetbelgium.com>
# Sent: Tuesday, December 3, 2019 1:46 AM
# 
# My ideal format is
# 
# Date, SZA, T at 3 wavelengths (440, 1020,1640nm), surface albedo at these 3 wavelengths, your results for COT, LWP,a_ef(spectral method)
# 
#  

# In[ ]:




