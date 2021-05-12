#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Info
def __init__():
    """
Purpose:

    Collections of functions to read in, analyze and plot 5STAR radiometer data

Input:

    See each function

Output:

    Figures and dict

Keywords:

    none

Dependencies:

    - load_utils.py
    - matplotlib
    - numpy
    - write_utils
    - path_utils
    - scipy

Needed Files:
  - 5STAR * RADIOMETERS.dat files

Modification History:

    Written: Samuel LeBlanc, Santa Cruz, CA, 2021-05-12, based on muSSTAR_rad_test_data.ipynb
    Modified:
    """
    pass


# In[5]:


def reader_RADIOMETER(file_in,lat = 36.970087, lon = -122.052537, alt = 74.0,islabjack=False):
    """
Purpose:
    Read in the 5STARG_RADIOMETERS*.dat files, and do preliminary calculations

Input:
    file_in: full path of the file to read

Output:
    Dict with headers as variables, 
        including UTC calculations as datetime objects, 
        sza, azi angles, 
        ref voltage values, 
        calculate voltages above reference

Keywords:
    lat: latitude of meaurements for sza calc (defaults to Conrad's house, but does not need to be exact)
    lon: longitude of measurements for sza calc  (defaults to Conrad's house, but does not need to be exact)
    alt: altitude of measurements for sza calc  (defaults to Conrad's house, but does not need to be exact)
    islabjack: (default False) If True, computes the voltages without the 5v reference, for use with labjack instead of ni

Dependencies:
    - pandas
    - datetime
    - numpy
    - map_utils

Needed Files:
  - file_in

Modification History:
    Written: Samuel LeBlanc, Santa Cruz, CA, 2020-11-18
    Modified: Samuel LeBlanc, Santa Cruz, CA, 2021-04,06, added handling for labjack
    Modified: Samuel LeBlanc, Santa Cruz, CA, 2021-05-12, adapted to python 3
    
    """
    import pandas as pd
    import map_utils as mu
    from parse_and_move_incoming_fx import get_filters_from_json, filetypes
    from pathlib2 import Path
    import dateutil.parser
    from datefinder import find_dates
    from datetime import date, datetime
    import json, re
    import numpy as np
    
    # read in
    s = pd.read_csv(file_in,header=5)
    
    # sanitize header names (remove spaces)
    for k in list(s.keys()):
        if k.strip().startswith(('CH','YY','MM','DD','J_')):
            s[k.strip()] = s.pop(k)
        else:
            s[k.strip().lower()] = s.pop(k)
    
    # get file metadata
    filters = get_filters_from_json('/data/sunsat/_incoming_gdrive/')
    ft = filetypes(file_in,dirlabel=None,dirdate=None,filters=filters)
    s.instname = ft.instname
    if s.instname in ['5STAR','5STARG']:
        s.label = ft.label.replace('_RADIOMETERS','')
    else:
        s.label = ft.label
    s.daystr = ft.daystr
    
    # calculate time and sun position
    s['UTC'] = pd.to_datetime(s['hh:mm:ss'],format='%H:%M:%S.%f')
    s['UTC'] = s['UTC'] + pd.offsets.DateOffset(years=s['YYYY'][0]-s['UTC'][0].year,
                                            months=s['MM'][0]-s['UTC'][0].month,
                                            days=s['DD'][0]-s['UTC'][0].day)
    s['sza'],s['azi'],s['sf'],s['dec'] = mu.get_sza_azi(lat,lon,s['UTC'],alt=alt,return_sunf_and_dec=True)
    
    # now get the reference voltages
    s['ref5v'] = s.get('-5v_0_66',np.zeros_like(s['sza']))*3.0
    s['ref12v_pos'] = s['12pos_0_33']*3.0
    s['ref12v_neg'] = s['12neg_0_33']*3.0
    s['5vtef_tempC'] = s['5vref_temp']/0.0021-273.15
    
    # run some calculations
    if '1020_outa' in list(s.keys()): s['nir_1020_teca'] = s['1020_outa']
    if '1020_outb' in list(s.keys()): s['nir_1020_tecb'] = s['1020_outb']
    if 'nir_1020_tecb' in list(s.keys()): s['nir_tec_V'] = s['nir_1020_tecb'] - s['nir_1020_teca']
    v2c = 100.0 #  "Temp = 10mv/DegC"
    if 'nir_block' in list(s.keys()):
        if (min(s['nir_block']) > 0.0) & (max(s['nir_block']) < 1.0):
            s['nir_block'] = s['nir_block']*v2c
            s['nir_board'] = s['nir_board']*v2c
            s['vis_block'] = s['vis_block']*v2c
            s['vis_board'] = s['vis_board']*v2c
            try:
                s['mezz'] = s['mezz']*v2c
            except:
                print('issue with the mezz need to multiply by {}'.format(v2c))
    if 'nir_block_temp' in list(s.keys()):
        if (min(s['nir_block_temp']) > 0.0) & (max(s['nir_block_temp']) < 1.0):
            s['nir_block'] = s['nir_block_temp']*v2c
            s['nir_board'] = s['nir_board_temp']*v2c
            s['vis_block'] = s['vis_block_temp']*v2c
            s['vis_board'] = s['vis_board_temp']*v2c
            try:
                s['mezz'] = s['mezz']*v2c
            except:
                print('issue with the mezz need to multiply by {}'.format(v2c))
                
    if 'J7_1' in list(s.keys()):
        s['5Vref_power'] = s['J7_1']
        s['cold_block_tempmon'] = s['J7_2']
        s['hot_block_tempmon'] = s['J7_3']
        s['cold_block_T'] = s['J7_4']
        s['nir_board_T'] = s['J7_5']
        s['hot_block_T'] = s['J7_6']
        s['air'] = s['J7_7']
        s['lidT'] = s['J7_8']
    
    print(len(s['UTC']),len(s['ref5v']))
    
    # now run through and calculate the voltages  compared to reference voltages
    V = np.zeros((9,3,len(s['ref5v'])))+np.nan
    for k in list(s.keys()):
        if k.startswith('CH') and k[-1].isdigit():
            if islabjack:
                s[k.replace('CH','V')] = -s[k]
            else:
                s[k.replace('CH','V')] = -s[k]-s['ref5v']
            # make array for easier use. channel, gain stage, then measurements
            icht,istaget = [int(i) for i in k.replace('CH','').split('_')]
            V[icht-1,int(istaget/2),:] = s[k.replace('CH','V')]
    define_long_names(s)
    #s['daystr'] = s['YYYY'][0]+s['MM'][0]+s['DD'][0]
    return s


# In[ ]:


def define_long_names(s):
    'function to give long names of common'    
    names = {
        'UTC': 'Datetime object of measurement time in UTC',
        'Year': 'Year of the measurement acquisistion',
        'YYYY': 'Year of the measurement acquisistion',
        'MM': 'Month of the measurement acquisistion',
        'DD': 'Day of the month of the measurement acquisistion',
        'hh:mm:ss': 'Hour minute seconds from the instrument time (in UTC) of the measurement acquisistion',
        'CH1_0':'340 nm 10x',
        'CH1_2':'340 nm 100x',
        'CH1_4':'340 nm 100x 100x',
        'CH2_0':'440 nm 1x',
        'CH2_2':'440 nm 100x',
        'CH2_4':'440 nm 100x 100x',
        'CH3_0':'675 nm 1x',
        'CH3_2':'675 nm 100x',
        'CH3_4':'675 nm 100x 100x',
        'CH4_0':'870 nm 1x',
        'CH4_2':'870 nm 100x',
        'CH4_4':'870 nm 100x 100x',
        'CH5_0':'1020 nm (Si) 1x',
        'CH5_2':'1020 nm (Si) 100x',
        'CH5_4':'1020 nm (Si) 100x 100x',
        'CH6_0':'1020 nm (InGaAs) 1x',
        'CH6_2':'1020 nm (InGaAs) 100x',
        'CH6_4':'1020 nm (InGaAs) 100x 100x',        
        'CH7_0':'1240 nm 1x',
        'CH7_2':'1240 nm 100x',
        'CH7_4':'1240 nm 100x 100x',
        'CH8_0':'1640 nm 1x',
        'CH8_2':'1640 nm 100x',
        'CH8_4':'1640 nm 100x 100x',
        'CH9_0':'2200 nm 1x',
        'CH9_2':'2200 nm 100x',
        'CH9_4':'2200 nm 100x 100x',
        'spare':'Spare - Empty [V]',
        '5vref_temp':'5Vref Temp, 2.1mV/DegC [V]',
        '5vtef_tempC':'Temp of 5Vref [\N{DEGREE SIGN}C]',
        '12pos_0_33':'1/3 of +12V rail [V/3]',
        '12neg_0_33':'1/3 of -12v rail [V/3]',
        '-5v_0_66':'1/3 of our -5v reference signal [V/3]',
        'ref5v':'calculated -5V reference signal [V]',
        'ref12v_pos':'calculated +12V reference signal [V]',
        'ref12v_neg':'calculated -12V reference signal [V]',
        'V':'Calculated voltage output [V], reference adjusted',
        'lat':'Latitude of measurement [\N{DEGREE SIGN}]',
        'lon':'Longitude of measurement [\N{DEGREE SIGN}]',
        'alt':'Altitude of measurement [m]',
        'sza':'Solar Zenith Angle [\N{DEGREE SIGN}], calculated',
        'azi':'Solar Azimuth Angle [\N{DEGREE SIGN}], calculated',
        'sf':'Solar-Earth factor due to distance between Earth and sun [unitless], calculated',
        'dec':'Solar declination per respect to earth [\N{DEGREE SIGN}],calculated',
        'nir_block':'NIR Block temp [\N{DEGREE SIGN}C]',
        'nir_board':'NIR Board temp [\N{DEGREE SIGN}C]',
        'vis_block':'VIS Block temp [\N{DEGREE SIGN}C]',
        'vis_board':'VIS Board temp [\N{DEGREE SIGN}C]',
        'mezz':'Mezzanine board temp [\N{DEGREE SIGN}C]',
        'nir_tec_V':'TEC difference voltage [V]',
        'air':'Ambient air temp [\N{DEGREE SIGN}C]',
        'J7_1' : '5v Power supply monitoring',
        'J7_2' : 'Cold block controller temp monitor',
        'J7_3' : 'Hot block controller temp monitor',
        'J7_4' : 'Cold block temp sensor',
        'J7_5' : 'NIR Board temp sensor',
        'J7_6' : 'Hot block temp sensor',
        'J7_7':'Ambient Temperature',
        'J7_8':'Lid Temperature'}
    for k in list(names.keys()):
        if k.startswith('CH'):
            names[k.replace('CH','V')] = names[k]+' [V]'
            names[k] = names[k] + ' [rawV]'
    for k in list(s.keys()):
        s[k].name = names.get(k,'Undefined')


# In[ ]:


def plot_housekeeping(s):
    'function to make the housekeeping plot'
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import numpy as np
    myFmt = mdates.DateFormatter('%H:%M:%S')
    
    
    fig = plt.figure(dpi=100)
    plt.plot(s['UTC'],s['spare'])
    plt.plot(s['UTC'],s['5vref_temp'])
    plt.plot(s['UTC'],s['ref5v'],label='(-5v_0_66)*3:\n {:02.5f}+/-{:02.5f}'.format(                                        np.nanmean(s['ref5v']),np.nanstd(s['ref5v'])))
    plt.plot(s['UTC'],s['12pos_0_33'])
    plt.plot(s['UTC'],s['12neg_0_33'])
    plt.legend()
    plt.xlabel('UTC [Hour]')
    plt.ylabel('Voltage [V]')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(7))
    try:
        fig.gca().set_title('{f.instname} - {f.daystr} {f.label}- Raw voltages'.format(f=s))
    except:
        pass

    fig1,ax = plt.subplots(4,1,sharex=True,dpi=100)
    ax[2].plot(s['UTC'],s['ref5v'],'.',label='5V ref: {:6.4f} +/- {:6.4f} V'.format(np.nanmean(s['ref5v']),np.nanstd(s['ref5v'])))
    ax[2].xaxis.set_major_formatter(myFmt)
    ax[2].grid()
    ax[2].legend()
    #ax[2].set_xlabel('UTC [Hour]')
    ax[2].set_ylabel('-5V ref \n[V]')

    ax[0].plot(s['UTC'],s['ref12v_pos'],'.',label='+12V ref: {:6.4f} +/- {:6.4f} V'.format(                                    np.nanmean(s['ref12v_pos']),np.nanstd(s['ref12v_pos'])))
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].grid()
    ax[0].legend()
    #ax[0].set_xlabel('UTC [Hour]')
    ax[0].set_ylabel('+12V ref \n[V]')
       
    ax[1].plot(s['UTC'],s['ref12v_neg'],'.',label='-12V ref: {:6.4f} +/- {:6.4f} V'.format(                                    np.nanmean(s['ref12v_neg']),np.nanstd(s['ref12v_neg'])))
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].grid()
    ax[1].legend()
    #ax[1].set_xlabel('UTC [Hour]')
    ax[1].set_ylabel('-12V ref \n[V]')
    
    ax[3].plot(s['UTC'],s['5vref_temp']/0.0021-273.15,'.',label='-5V Temperature ref: {:6.4f} +/- {:6.4f} $^{{\circ}}$C'.format(                                    np.nanmean(s['5vref_temp']/0.0021-273.15),                                    np.nanstd(s['5vref_temp']/0.0021-273.15)))
    ax[3].xaxis.set_major_formatter(myFmt)
    ax[3].grid()
    ax[3].legend()
    ax[3].set_xlabel('UTC [Hour]')
    ax[3].set_ylabel('Temperature\n[$^{{\circ}}$C]')
    fig1.tight_layout(rect=[0,0,1,0.96])
    ax[3].xaxis.set_major_locator(plt.MaxNLocator(7))
    try:
        ax[0].set_title('{f.instname} - {f.daystr} {f.label}- Calculated input voltages'.format(f=s))
    except:
        pass
    
    if 'vis_block' in list(s.keys()):
        fig2,ax2 = plt.subplots(1,1,dpi=100)
        plt.plot(s['UTC'],s['vis_block'])
        plt.plot(s['UTC'],s['nir_block'])
        plt.plot(s['UTC'],s['vis_board'])
        plt.plot(s['UTC'],s['nir_board'])
        if 'mezz' in list(s.keys()):
            plt.plot(s['UTC'],s['mezz'])
        if 'air' in list(s.keys()):
            plt.plot(s['UTC'],s['air'])
        plt.grid()
        plt.legend()
        plt.xlabel('UTC [Hour]')
        plt.ylabel('Temperature [\N{DEGREE SIGN}C]')
        fig2.tight_layout(rect=[0,0,1,0.96])
        ax2.xaxis.set_major_formatter(myFmt)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(7))
        try:
            ax2.set_title('{f.instname} - {f.daystr} {f.label}- Temperatures'.format(f=s))
        except:
            pass
    else:
        fig2 = None
        
    if 'nir_tec_V' in list(s.keys()):
        fig3,ax = plt.subplots(2,1,sharex=True,dpi=100)
        if 'nir_1020_teca' in list(s.keys()):
            ax[0].plot(s['UTC'],s['nir_1020_teca'],'.')
            ax[0].plot(s['UTC'],s['nir_1020_tecb'],'.')
        if 'CH8_TEC_B' in list(s.keys()):
            ax[0].plot(s['UTC'],s['CH8_TEC_A'],'.')
            ax[0].plot(s['UTC'],s['CH8_TEC_B'],'.')
        ax[1].plot(s['UTC'],s['nir_tec_V'],'-',label='tecb - teca',lw=0.2)
        ax[0].legend()
        ax[1].legend()
        ax[1].set_xlabel('UTC [Hour]')
        ax[0].set_ylabel('TEC Voltages [V]')
        ax[1].set_ylabel('TEC Difference [V]')
        ax[0].set_title('{s.instname} - {s.daystr} {s.label}- NIR 1020 nm TEC Voltages'.format(s=s))
        ax[1].xaxis.set_major_locator(plt.MaxNLocator(7))
        myFmt = mdates.DateFormatter('%H:%M:%S')
        ax[1].xaxis.set_major_formatter(myFmt)
        ax[1].grid()
        ax[0].grid()
    else:
        fig3 = None

    return [fig,fig1,fig2,fig3]


# In[ ]:


def plot_v(s,gain=0):
    'function to plot out the different gain channels'
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(2,1,sharex=True,dpi=100)
    
    ig = ['0','2','4']
    ig_str = ['1x','100x','100x-100x']
    for ut in ['V1_','V2_','V3_','V4_','V5_']:
        u = ut+ig[gain]
        ax[0].plot(s['UTC'],s[u],label=s[u].name)
    #ax[0].set_title('5STARG - test {} - VIS {} Channels'.format(s['daystr'],ig_str[gain]))
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_ylim(-0.5,10.5)

    for ut in ['V6_','V7_','V8_','V9_']:
        u = ut+ig[gain]
        ax[1].plot(s['UTC'],s[u],label=s[u].name)
    ax[1].set_title('NIR {} Channels'.format(ig_str[gain]))
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylim(-0.5,10.5)
    ax[1].xaxis.set_major_locator(plt.MaxNLocator(7))
    try:
        ax[0].set_title('{f.instname} - {f.daystr} {f.label}- VIS {} Channels'.format(ig_str[gain],f=s))
    except:
        pass
    
    return fig,ax


# In[ ]:


def plot_v_expect(s,gain=0):
    'function to plot out the different gain channels, gain refers to the base gain'
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(2,1,sharex=True,figsize=(9,6),dpi=100)
    
    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    
    ig = ['0','2','4']
    ig_str = ['1x','100x','100x-100x']
    for iu,ut in enumerate(['V1_','V2_','V3_','V4_','V5_']):
        u = ut+ig[gain]
        u1 = ut+ig[gain+1]
        ax[0].plot(s['UTC'],s[u]*100.0,label=s[u].name.split()[0]+' expect',c=cs[iu],alpha=0.4,ls=':')
        ax[0].plot(s['UTC'],s[u1],label=s[u1].name,c=cs[iu],ls='-')
    ax[0].set_title('VIS {} Channels'.format(ig_str[gain+1]))
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V]')
    ax[0].grid()
    ax[0].legend(bbox_to_anchor=(1.05, 1.0, 0.3, 0.2), loc='upper left')
    ax[0].set_ylim(-0.5,10.5)

    for iu,ut in enumerate(['V6_','V7_','V8_','V9_']):
        u = ut+ig[gain]
        u1 = ut+ig[gain+1]
        ax[1].plot(s['UTC'],s[u]*100.0,label=s[u].name.split()[0]+' expect',c=cs[iu],alpha=0.4,ls=':')
        ax[1].plot(s['UTC'],s[u1],label=s[u1].name,c=cs[iu],ls='-')
    ax[1].set_title('NIR {} Channels'.format(ig_str[gain+1]))
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend(bbox_to_anchor=(1.05, 1.0, 0.3, 0.2), loc='upper left')
    ax[1].set_ylim(-0.5,10.5)
    
    plt.tight_layout()
    ax[1].xaxis.set_major_locator(plt.MaxNLocator(7))

    try:
        ax[0].set_title('{f.instname} - {f.daystr} {f.label}- Radiometer expected voltages'.format(f=s))
    except:
        pass
    
    return fig,ax


# In[ ]:


def plot_channel_multigain(s,ch=1,ns=300):
    'Plotting of a single wavelength (channel) over multiple gain stages. s:5STAR pd, ch: channel number, ns: smoothing window'
    from Sp_parameters import smooth
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(7,7),dpi=100)

    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    ic = (ch-1)%5
    
    u,u1,u2 = 'V{}_0'.format(ch),'V{}_2'.format(ch),'V{}_4'.format(ch)
    ax[0].plot(s['UTC'],s[u],'.',label='raw, mean={m:2.5f}, median={e:2.5f}, std={s:2.5f}'.format(
               m=s[u].mean(),e=s[u].median(),s=s[u].std()),c=cs[ic],alpha=0.02)
    ax[0].plot(s['UTC'],smooth(s[u],ns,old=True),'-',label='smooth n={}'.format(ns),c='k')
    ax[0].legend()
    ax[1].plot(s['UTC'],s[u1],'.',label='raw, mean={m:2.5f}, median={e:2.5f}, std={s:2.5f}'.format(
               m=s[u1].mean(),e=s[u1].median(),s=s[u1].std()),c=cs[ic],alpha=0.02)
    ax[1].plot(s['UTC'],smooth(s[u1],ns,old=True),'-',c='k')
    ax[1].legend()
    ax[2].plot(s['UTC'],s[u2],'.',label='raw, mean={m:2.5f}, median={e:2.5f}, std={s:2.5f}'.format(
               m=s[u2].mean(),e=s[u2].median(),s=s[u2].std()),c=cs[ic],alpha=0.02)
    ax[2].plot(s['UTC'],smooth(s[u2],ns,old=True),'-',c='k')
    ax[2].legend()

    ax[0].set_title('{s.instname} - {s.daystr} - {u.name}'.format(u=s[u],s=s))
    ax[1].set_title('{u.name}'.format(u=s[u1]))
    ax[2].set_title('{u.name}'.format(u=s[u2]))
    ax[2].xaxis.set_major_formatter(myFmt)
    #ax[2].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V] (1x)')
    ax[1].set_ylabel('Irradiance [V] (100x)')
    ax[2].set_ylabel('Irradiance [V] (100x-100x)')
    ax[2].set_xlabel('UTC [H]')
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    plt.tight_layout()
    
    return fig, ax


# In[ ]:


def plot_gains(s,ch,ns=300):
    'Plotting the ratio between the different gain stages for one wavelength'
    from Sp_parameters import smooth
    import plotting_utils as pu
    import matplotlib.pyplot as plt

    fig,ax = plt.subplots(1,2,figsize=(7,3),dpi=100)
    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    ic = (ch-1)%5
    igood = (s['V{}_0'.format(ch)] < 10.3) & (s['V{}_2'.format(ch)] < 10.3) & (s['V{}_2'.format(ch)] > 0.0)
    ax[0].plot(smooth(s['V{}_0'.format(ch)][igood],ns,old=True),s['V{}_2'.format(ch)][igood],'.',
               alpha=0.02,label=s['V{}_0'.format(ch)].name,c=cs[ic])
    pu.plot_lin(smooth(s['V{}_0'.format(ch)][igood],ns,old=True),s['V{}_2'.format(ch)][igood],
                lblfmt='2.5f',ax=ax[0],color=cs[ic])
    ax[0].legend()
    ax[0].set_ylabel(s['V{}_2'.format(ch)].name)
    ax[0].set_xlabel('smooth ('+s['V{}_0'.format(ch)].name+')')
    ax[0].xaxis.set_major_locator(plt.MaxNLocator(5))

    igood2 = (s['V{}_2'.format(ch)] < 10.3) & (s['V{}_4'.format(ch)] < 10.3) & (s['V{}_4'.format(ch)] > 0.0)
    ax[1].plot(smooth(s['V{}_2'.format(ch)][igood2],ns,old=True),s['V{}_4'.format(ch)][igood2],'.',
               alpha=0.02,label=s['V{}_2'.format(ch)].name,c=cs[ic])
    pu.plot_lin(smooth(s['V{}_2'.format(ch)][igood2],ns,old=True),s['V{}_4'.format(ch)][igood2],
                lblfmt='2.5f',ax=ax[1],color=cs[ic])
    ax[1].legend()
    ax[1].set_ylabel(s['V{}_4'.format(ch)].name)
    ax[1].set_xlabel('smooth ('+s['V{}_2'.format(ch)].name+')')
    ax[1].xaxis.set_major_locator(plt.MaxNLocator(5))

    fig.tight_layout()
    fig.suptitle('V{}'.format(ch))
    return fig,ax


# In[ ]:


def plot_corr_temp(s,ig=0,tmp='nir_block'):
    'Plot some correlation plots of the voltages vs. temp'
    import plotting_utils as pu
    import matplotlib.pyplot as plt
    #tmp = ['nir_block','nir_board','vis_block','vis_board','mezz'][it]
    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    fig,ax = plt.subplots(3,3,sharex=True,dpi=100)
    for ia, x in enumerate(ax.ravel()):
        goods = (s['V{}_{}'.format(ia+1,ig)]<10.3) & (s['V{}_{}'.format(ia+1,ig)]>0.0) 
        r = np.corrcoef(s[tmp][goods],s['V{}_{}'.format(ia+1,ig)][goods])[0,1]
        x.plot(s[tmp][goods],s['V{}_{}'.format(ia+1,ig)][goods],'.',alpha=0.01,label='R$^2$={:1.4f}'.format(r*r),c=cs[ig])
        pu.plot_lin(s[tmp][goods],s['V{}_{}'.format(ia+1,ig)][goods],ax=x,lblfmt='1.4f',color=cs[ig])
        x.legend(frameon=False, prop={'size': 6})
        x.set_ylabel('V{}_{}, '.format(ia+1,ig)+s['V{}_{}'.format(ia+1,ig)].name.split(' ')[0]+' [V]')
        if (ia+1)>6: 
            x.set_xlabel(s[tmp].name)
    fig.suptitle('{s.instname} - {s.daystr} {s.label}- Correlation with temp'.format(s=s))
    fig.tight_layout(rect=(0,0,1,0.97))
    return fig


# In[ ]:


def plot_channels(s,fp,ns=300,dpi=600):
    'Run through each channel to plot the multigains'
    for i in range(9):
        fig,ax = plot_channel_multigain(s,ch=i+1,ns=ns)
        u = 'V{}_0'.format(i+1)
        fig.savefig(fp+'{s.instname}_{s.daystr}_rad_{wvl}_allch.png'.format(s=s,wvl=s[u].name.split()[0]),
                    dpi=dpi,transparent=True)


# In[ ]:


def plot_all_gainratios(s,fpp,ns=300,dpi=300):
    'Run through and plot each wavelenght and gain stage ratios'
    for i in range(9):
        try:
            fig,ax = plot_gains(s,i+1,ns=ns)
            u = 'V{}_0'.format(i+1)
            fig.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_{wvl}_gainratios.png'.format(s=s,wvl=s[u].name.split()[0]),
                        dpi=dpi,transparent=True)
        except:
            print('issue with gain ratios for wvl {wvl}'.format(wvl=s[u].name.split()[0]))
            pass


# In[ ]:


def clean_up_and_prep_for_plots(s,fpp):
    'to clean up the label and make the path for plotting, and other pre-plotting setups'
    import os
    if (len(s.label.split('_'))>1) & s.label.split('_')[0].isdigit():
        s.label = s.label.split('_')[-1]
    if not os.path.isdir(fpp+'{s.daystr}_{s.label}/'.format(s=s)):
        os.makedirs(fpp+'{s.daystr}_{s.label}/'.format(s=s))
    define_long_names(s)


# In[ ]:


def plt_hskp(s,fpp,dpi=300):
    'plot all the housekeeping files and save them'
    import numpy as np
    fig = plot_housekeeping(s)
    fig[0].savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=dpi,transparent=True)
    fig[1].savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=dpi,transparent=True)
    try:
        fig[2].savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        pass
    try:
        fig[3].savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_NIR_TECvoltages.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        pass
    #return fig


# In[ ]:


def plt_gains(s,fpp,dpi=300):
    'plot all the different channels seperated by the gains'
    try:
        fig1,ax1 = plot_v(s,gain=0)
        ax1[0].autoscale(axis='y')
        ax1[1].autoscale(axis='y')
        fig1.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        print('issue with gain 0, 1x')
        pass
        
    try:
        fig2,ax2 = plot_v(s,gain=1)
        ax2[0].autoscale(axis='y')
        ax2[1].autoscale(axis='y')
        fig2.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        print('issue with gain 0, 100x')
        pass
    
    try:
        fig3,ax3 = plot_v(s,gain=2)
        ax3[0].autoscale(axis='y')
        ax3[1].autoscale(axis='y')
        fig3.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        print('issue with gain 0, 100x - 100x')
        pass
    
    try:
        fig4,ax4 = plot_v_expect(s,gain=0)
        #plt.title('5STARG radiometers 1x 2020-11-17')
        ax4[0].autoscale(axis='y')
        ax4[1].autoscale(axis='y')
        fig4.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        print('issue with gain 0, 100x, expect')
        pass
        
    try:
        fig5,ax5 = plot_v_expect(s,gain=1)
        #plt.title('5STARG radiometers 1x 2020-11-17')
        ax5[0].autoscale(axis='y')
        ax5[1].autoscale(axis='y')
        fig5.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=dpi,transparent=True)
    except:
        print('issue with gain 0, 100x - 100x expect')
        pass
    
    #return [fig1,fig2,fig3,fig4,fig5],[ax1,ax2,ax3,ax4,ax5]


# In[ ]:


def plt_corrs(s,fpp,dpi=300):
    'Run through the different gains and temperatures'
    for ig in [0,2,4]:
        for it in ['nir_block','vis_block','nir_tec_V']:
            try:
                fig = plot_corr_temp(s,ig,it)
                fig.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_corr_gain{}_temp{}.png'.format(ig,it,s=s),
                             dpi=dpi,transparent=True)
            except:
                print('problem with corr gain {} temp {}'.format(ig,it))

