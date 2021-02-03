#!/usr/bin/env python
# coding: utf-8

# # Info
# Purpose:
# 
#     Analyse initial dataset from the radiometer sun measurements from Conrad of muSSTAR
# 
# Input:
# 
#     none
# 
# Output:
# 
#     Figures 
# 
# Keywords:
# 
#     none
# 
# Dependencies:
# 
#     - load_utils.py
#     - matplotlib
#     - numpy
#     - Sp_parameters
#     - write_utils
#     - path_utils
#     - hdf5storage
#     - scipy
# 
# Needed Files:
#   - file.rc : for consistent creation of look of matplotlib figures
#   - ...
# 
# Modification History:
# 
#     Written: Samuel LeBlanc, Santa Cruz, CA, 2020-10-28
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2021-01-15
#               Added new test for dark drift.
#     Modified: Samuel LeBlanc, Santa Cruz, CA, 2021-01-26
#               Added new dark test with temp sensors.
# 

# # Notes

# From Conrad Esch, Oct 16, 2020, 5:01 PM:
# Hi All,
# 
# I got 5STAR going today and managed to grab some data. The data is included in the text file. The text file is quite large so I'm open to alternative data storage methods or loop speeds.
# 
# The 0 to -10v signals were referenced to -5v so consider ~5 volts to be 0 and ~-5v to be -10v.
# 
# I ran the capture loop at 10 hz and it was started at 3:39 pm.
# 
# The instrument was manually placed so that it was looking directly at the sun and then was left alone as the sun's position changed. Also included is a photo.
# 
# Best,
# 
# Conrad Esch
# ![5STAR10.16.jpg](attachment:5STAR10.16.jpg)
# 

# Conrad,  can you tell me what the labels represent?  I see:
# 
# CH1_0, CH1_2, CH1_4: (Are these the same detector but with different gains and/or dynamic range?)
# 
# CH2_0, CH2_2, CH2_4,
# 
# CH3_0, CH3_2, CH3_4,
# 
# CH4_0, CH4_2, CH4_4,
# 
# CH5_0, CH5_2, CH5_4,
# 
# -5v_0_66:  (This signal has a mean value of -1.673 with standard deviation 0.0197 ~= 0.002.)
# 
# CH6_0, CH6_2, CH6_4,
# 
# CH7_0, CH7_2, CH7_4,
# 
# CH8_0, CH8_2, CH8_4,
# 
# CH9_0, CH9_2, CH9_4,
# 
# Spare,
# 
# 5vRef_Temp  (mean value of 0.554, standard deviation 0.0197)
# 
# 12POS_0_33 (mean value 0f 3.9382, standard deviation of 0.100)
# 
# 12NEG_0_33 (mean value -3.9281, standard dev 0.143)
# 
#  
# 
# So, I'd guess that these are sequences of the same detector (CHx_0, CHx_2, CHx_4) recorded with different gains and ranges?  Do you have a mapping of the channel number to wavelength?  And -5v_0_66, 5vRef_Temp, 12POS_0_33 and 12NEG_0_33 are reference voltages?
# 
#  
# 
# Was the 5STAR supposedly tracking the sun most of the time?  Some of these traces appear to track with each other and others not so well. 
# 
#  
# 
# Do we have any blocked, shuttered, or dark measurements?
# 
#  
# 
# Here is a quick plot of all 32 signals (including references) after removing the minimum value as a baseline.   Time is on the horizontal.
# 
# ![image001.png](attachment:image001.png)

# From Conrad Esch, 
# Oct 20, 2020, 12:57 PM
# 
# Hi Connor,
# 
# Sorry for confusing labelling, yes the ones you grouped together are the same channels but with different gains. 
# 
# CHX_Y is Channel Number X times 10^Y. So for example CH1_0 is Channel1x1 and CH1_4 is Channel1x100x100. I'll write out what I believe is the correct wavelengths and channel descriptions below.
# 
# Radiometer Channels - Visible Bank
#  - CH1_Y - 340-2 Channel
#  - CH2_Y - 440-10 Channel
#  - CH3_Y - 675-10 Channel
#  - CH4_Y - 870-10 Channel
#  - CH5_Y - 1020 Channel
# 
# Radiometer Channels - NIR Bank
#  - CH6_Y - 1020 Channel
#  - CH7_Y - 1240 Channel
#  - CH8_Y - 1640 Channel
#  - CH9_Y - 2200 Channel
# 
# Housekeeping Channels
#  - Spare - Should be floating and unused but the signal appears similar to some of the other channels so it needs to be investigated.
#  - -5v_0_66 - Should actually be -5v_0_33, it is 1/3 of our -5v reference signal, the channel name was supposed to represent -5v*0.333
#  - 5vRef_Temp - Attached to  a temperature output from 5v precision reference, conversion is 2.1mV/Degree Celsius, datasheet says typical value is 630mV
#  - 12POS_0_33 - 1/3 of +12V rail
#  - 12NEG_0_33 - 1/3 of -12v rail
# 
# In this data set the Radiometer channels were all referenced to -5v to try to provide a higher bit resolution. Normally the channel outputs would range from 0v to -10v where 0v is dark and -10v is light. With the reference to -5v, a signal of 5v is dark and -5v is light. I haven't tuned the -5v reference yet so it is more like -5.019v.
# 
# In terms of linearity there is an issue I'm trying to track down with the third gain stages for the channels (CHX_4). So at this time I would ignore those gain stages for every channel.
# 
# The instrument was not tracking for this test, I pointed the instrument by hand toward the sun and then let the sun drift out of line with the instrument over time. We are not quite set up for tracking yet.
# 
# I haven't done any dark measurements or drift/noise tests yet. Friday's test was more to show major system/integration/procedural flaws. Once I've worked out the issues with the third gain stage I'll proceed to other tests.
# 
# Let me know if you have any questions.
# 
# Best,
# 
# Conrad Esch

# ## New Test, 2020-11-04
# I still haven't figured out the issues with the 3rd gain stage but I thought we could evaluate other aspects of the instrument and test procedures while I'm sending emails back and forth with Analog Devices.
# 
# Before the tests occurred i rebiased all of the op-amps and adjusted the 5v offset.
# 
# I got 4 tests:
# 1. Sweep Test: Instrument was pointed at the sun at 12:38PM (43 minutes after solar noon). It was rotated off of the sun by hand so that the sun would pass across the sensors.
# 
# 2. MoveOffTest: Instrument was pointed at the sun at 1:20 PM and then the sun was allowed to move away from alignment with the instrument, same as the original test a few weeks ago.
# 
# 3. Dark Test outside: Tinfoil was placed over the gershun tubes and the instrument was pointed down immediately after the MoveOffTest ended. I noticed some readings that seemed high so i performed a second dark test indoors. I think there might be a light leakage problem coming from the interface between the photodiodes and the gershun tubes.
# 
# 4. Dark Test indoors: Same as test 3 but indoors to remove natural light sources.
# 
# Hope this data is helpful, I will talk more about it at the meeting.
# 
# Best,
# 
# Conrad Esch

# ## Test 2021-01-14
# 
# 
# I got a 2 hour drift/noise test done today. I managed to hook up the TEC drivers so the NIR photodies were cooled. This test was done with the Daburn/Polytron AC-DC power supply.
# 
# I biased the amplifiers prior to the test and every channel should have started with negative voltages across the 3 gains.
# 
# The blocks did heat up during the test, especially the VIS block which the TEC drivers were dumping their heat into. I had a hard time measuring the block temperature but It seemed warm to the touch at the end of the test. The board temperatures produced better results with the thermal temperature sensor. Sections of the board rose to a maximum of 96 degrees fahrenheit by the end of the test.
# 
# Watching the test, the values definitely changed. I suspect that it is due to temperature but it's possible that light leakage played a role. Early on in the test, the NIR x100x100 values railed beyond the ADC's measurement value. They should have been already cooled down by the time I biased them so I wonder if there is a board warmup time that should be considered.
# 
# If you would let me know what you think of the data and how I should re conduct the tests, I would appreciate it. Maybe the next step is a lot of temperature sensors.
# 
# Best,
# 
# Conrad Esch

# ## Test 2021-01-25
# I got the updated drift/noise test done last night with temperature sensors.
# 
# In this test the temperature of the boards and blocks actually went down, probably because the ambient temperature started dropping immediately after starting the test.
# 
# I think that part of the issue with the previous test was that a heater vent may have been blowing a small amount of warm air toward the system.
# 
# The next test I'll be running is the system on battery power with the TEC drivers turned off. After that I'll run a short test with a blow dryer pointed toward the system.
# 
# Best,
# 
# Conrad Esch
# 
# 
#  5STARG_20210125_154600_RADIOMETERS.dat

# ## Test 2021-01-26
# I got another test done this morning. It is also a 2 hour drift test.
# 
# This is probably a very low priority test for a few reasons:
# 
# This test was done with the battery powered PDN Board
# The batteries ran out of power part way through the test
# For this test I left the NIR TEC drivers off, so the NIR cans were uncooled (as per last meeting's discussion)
# 
# The next thing I will do is our normal test (like the one I sent you this morning) with a hair dryer pointed at the system. It should be a much shorter test.
# 
# If there are any other tests you think I should run, please let me know. The system is currently set up and idle so now is a good time. Setting up and taking down the system is a bit time consuming.
# 
# Thank you for doing all of the analysis. I know you're busy.
# 
# Best,
# 
# Conrad Esch
# 
#  5STARG_20210126_103658_RADIOMETERS.dat
#  
#  
#  
# ====
# I completed two additional tests today.
# 
# For Test 1, I used a hair dryer to warm the system by about 15-20 degrees.
# 
# For Test 2, I waved a group of incandescent lamps in front of the gershun tubes.
# 
# I'm guessing that I will have to redo some of these tests in terms of duration so let me know if there's any useful information in the tests.
# 
# Best,
# 
# Conrad
# 
# 5STARG_20210126_143946_RADIOMETERS_HeatTest.dat
# 
# 5STARG_20210126_154223_RADIOMETERS_IncandescentTest.dat
# 

# ## Test 2021-01-27
# I have completed the "Speaker Wire" test. Based on what I was watching on the screen, I don't think there is any useful information in this but I thought I would send it along anyways.
# 
# I played tones through the speakers starting with off and then on in 30 second increments.
# 
# The pattern would follow this.
# 
# off-100-off-200-off-300-off-400-off-500-off-600-off-700-off
# 
# Where 100 is a 100hz. tone.
# 
# Best,
# 
# Conrad Esch

# ## Test 2021-01-28
# Another test.
# 
# I hooked two adc pins of the labjack into the OUTA and OUTB pins of the 1020nm NIR TEC Driver. It has an OUTA and an OUTB because it can drive the TEC in either direction to warm or cool the photodiode.
# 
# To do the test I had to remove some of the temperature processing. The temp sensors were still there but the data in the file is in raw voltages instead of Deg C. The conversion is "Temp = 10mv/DegC" so to get the correct temps the raw voltages need to be multiplied by 100.
# 
# For the TEC driver voltages, I think the most useful way to look at the data will be to take the difference, "NIR_1020_TECB-NIR_1020_TECA", of those data columns.
# 
# Best,
# 
# Conrad Esch
# 
# 5STARG_20210128_155528_RADIOMETERS_TECDriverVoltages.dat
# 

# ## Test 2021-01-29

# Last friday I ran a test where I ran the system for 1 hour with the NIR detectors covered and then ran for another hour with the NIR detectors uncovered.
# 
# Included is the data.
# 
# I'm currently running a 6 hour drift test and will be digging into the spec sheets on our photodiodes looking for indicated differences between the NIR and VIS diodes.
# 
# Best,
# 
# Conrad Esch
# 
# 
# 5STARG_20210129_163011_RADIOMETERS_CoverONOFF.dat
# 

# ## Test 2021-02-01 6 hour drift test

# Got the 6 hour drift test done.
# 
# I think there is some very useful data in it for the NIR 100x channels.
# 
# It is a very large file, does this track with the file size for other systems like 4STAR?
# 
# Best,
# 
# Conrad Esch
# 
# 5STARG_20210201_120748_RADIOMETERS_6HOUR.dat
# 

# # Prepare python environment

# In[1]:


import numpy as np
import Sp_parameters as Sp
import load_utils as lu
import write_utils as wu
from path_utils import getpath
import hdf5storage as hs
import scipy.io as sio
import matplotlib.pyplot as plt
import os
import pandas as pd


# In[2]:


get_ipython().magic(u'matplotlib notebook')


# In[3]:


name = 'muSSTAR'
fp = getpath(name)


# # Make functions for faster analysis

# In[97]:


def reader_RADIOMETER(file_in,lat = 36.970087, lon = -122.052537, alt = 74.0):
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

Dependencies:
    - pandas
    - datetime
    - numpy
    - map_utils

Needed Files:
  - file_in

Modification History:
    Written: Samuel LeBlanc, Santa Cruz, CA, 2020-11-18
    Modified:
    
    """
    import pandas as pd
    import map_utils as mu
    from parse_and_move_incoming_fx import get_filters_from_json, filetypes
    from pathlib2 import Path
    import dateutil.parser
    from datefinder import find_dates
    from datetime import date, datetime
    import json, re
    
    # read in
    s = pd.read_csv(file_in,header=5)
    
    # sanitize header names (remove spaces)
    for k in s.keys():
        if k.strip().startswith(('CH','YY','MM','DD')):
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
    s['ref5v'] = s['-5v_0_66']*3.0
    s['ref12v_pos'] = s['12pos_0_33']*3.0
    s['ref12v_neg'] = s['12neg_0_33']*3.0
    s['5vtef_tempC'] = s['5vref_temp']/0.0021-273.15
    
    # run some calculations
    if 'nir_1020_tecb' in s.keys(): s['nir_tec_V'] = s['nir_1020_tecb'] - s['nir_1020_teca']
    v2c = 100.0 #  "Temp = 10mv/DegC"
    if 'nir_block' in s.keys():
        if (min(s['nir_block']) > 0.0) & (max(s['nir_block']) < 1.0):
            s['nir_block'] = s['nir_block']*v2c
            s['nir_board'] = s['nir_board']*v2c
            s['vis_block'] = s['vis_block']*v2c
            s['vis_board'] = s['vis_board']*v2c
            s['mezz'] = s['mezz']*v2c
    if 'nir_block_temp' in s.keys():
        if (min(s['nir_block_temp']) > 0.0) & (max(s['nir_block_temp']) < 1.0):
            s['nir_block'] = s['nir_block_temp']*v2c
            s['nir_board'] = s['nir_board_temp']*v2c
            s['vis_block'] = s['vis_block_temp']*v2c
            s['vis_board'] = s['vis_board_temp']*v2c
            s['mezz'] = s['mezz']*v2c
    
    print len(s['UTC']),len(s['ref5v'])
    
    # now run through and calculate the voltages  compared to reference voltages
    V = np.zeros((9,3,len(s['ref5v'])))+np.nan
    for k in s.keys():
        if k.startswith('CH'):
            s[k.replace('CH','V')] = -s[k]-s['ref5v']
            # make array for easier use. channel, gain stage, then measurements
            icht,istaget = [int(i) for i in k.replace('CH','').split('_')]
            V[icht-1,istaget/2,:] = s[k.replace('CH','V')]
    define_long_names(s)
    #s['daystr'] = s['YYYY'][0]+s['MM'][0]+s['DD'][0]
    return s


# In[219]:


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
        '5vtef_tempC':'Temp of 5Vref [°C]',
        '12pos_0_33':'1/3 of +12V rail [V/3]',
        '12neg_0_33':'1/3 of -12v rail [V/3]',
        '-5v_0_66':'1/3 of our -5v reference signal [V/3]',
        'ref5v':'calculated -5V reference signal [V]',
        'ref12v_pos':'calculated +12V reference signal [V]',
        'ref12v_neg':'calculated -12V reference signal [V]',
        'V':'Calculated voltage output [V], reference adjusted',
        'lat':'Latitude of measurement [°]',
        'lon':'Longitude of measurement [°]',
        'alt':'Altitude of measurement [m]',
        'sza':'Solar Zenith Angle [°], calculated',
        'azi':'Solar Azimuth Angle [°], calculated',
        'sf':'Solar-Earth factor due to distance between Earth and sun [unitless], calculated',
        'dec':'Solar declination per respect to earth [°],calculated',
        'nir_block':u'NIR Block temp [°C]',
        'nir_board':u'NIR Board temp [°C]',
        'vis_block':u'VIS Block temp [°C]',
        'vis_board':u'VIS Board temp [°C]',
        'mezz':u'Mezzanine board temp [°C]',
        'nir_tec_V':'TEC difference voltage [V]'}
    for k in names.keys():
        if k.startswith('CH'):
            names[k.replace('CH','V')] = names[k]+' [V]'
            names[k] = names[k] + ' [rawV]'
    for k in s.keys():
        s[k].name = names.get(k,'Undefined')


# In[86]:


def plot_housekeeping(s):
    'function to make the housekeeping plot'
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%H:%M:%S')
    
    
    fig = plt.figure()
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

    fig1,ax = plt.subplots(4,1,sharex=True)
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
    
    if 'mezz' in s.keys():
        fig2,ax2 = plt.subplots(1,1)
        plt.plot(s['UTC'],s['vis_block'])
        plt.plot(s['UTC'],s['nir_block'])
        plt.plot(s['UTC'],s['vis_board'])
        plt.plot(s['UTC'],s['nir_board'])
        plt.plot(s['UTC'],s['mezz'])
        plt.grid()
        plt.legend()
        plt.xlabel('UTC [Hour]')
        plt.ylabel(u'Temperature [°C]')
        fig2.tight_layout(rect=[0,0,1,0.96])
        ax2.xaxis.set_major_formatter(myFmt)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(7))
        try:
            ax2.set_title('{f.instname} - {f.daystr} {f.label}- Temperatures'.format(f=s))
        except:
            pass
    else:
        fig2 = None
        
    if 'nir_tec_V' in s.keys():
        fig3,ax = plt.subplots(2,1,sharex=True)
        ax[0].plot(s['UTC'],s['nir_1020_teca'],'.')
        ax[0].plot(s['UTC'],s['nir_1020_tecb'],'.')
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


# In[87]:


def plot_v(s,gain=0):
    'function to plot out the different gain channels'
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(2,1,sharex=True)
    
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


# In[88]:


def plot_v_expect(s,gain=0):
    'function to plot out the different gain channels, gain refers to the base gain'
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(2,1,sharex=True,figsize=(9,6))
    
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


# In[89]:


def plot_channel_multigain(s,ch=1,ns=300):
    'Plotting of a single wavelength (channel) over multiple gain stages. s:5STAR pd, ch: channel number, ns: smoothing window'
    from Sp_parameters import smooth
    import matplotlib.dates as mdates
    myFmt = mdates.DateFormatter('%H:%M:%S')
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(7,7))

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


# In[162]:


def plot_gains(s,ch,ns=300):
    'Plotting the ratio between the different gain stages for one wavelength'
    from Sp_parameters import smooth
    import plotting_utils as pu

    fig,ax = plt.subplots(1,2,figsize=(7,3))
    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    ic = (ch-1)%5
    igood = (s['V{}_0'.format(ch)] < 10.3) & (s['V{}_2'.format(ch)] < 10.3) 
    ax[0].plot(smooth(s['V{}_0'.format(ch)][igood],ns,old=True),s['V{}_2'.format(ch)][igood],'.',
               alpha=0.02,label=s['V{}_0'.format(ch)].name,c=cs[ic])
    pu.plot_lin(smooth(s['V{}_0'.format(ch)][igood],ns,old=True),s['V{}_2'.format(ch)][igood],
                lblfmt='2.5f',ax=ax[0],color=cs[ic])
    ax[0].legend()
    ax[0].set_ylabel(s['V{}_2'.format(ch)].name)
    ax[0].set_xlabel('smooth ('+s['V{}_0'.format(ch)].name+')')
    ax[0].xaxis.set_major_locator(plt.MaxNLocator(5))

    igood2 = (s['V{}_2'.format(ch)] < 10.3) & (s['V{}_4'.format(ch)] < 10.3) 
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


# In[223]:


def plot_corr_temp(s,ig=0,tmp='nir_block'):
    'Plot some correlation plots of the voltages vs. temp'
    import plotting_utils as pu
    #tmp = ['nir_block','nir_board','vis_block','vis_board','mezz'][it]
    cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
    fig,ax = plt.subplots(3,3,sharex=True)
    for ia, x in enumerate(ax.ravel()):
        goods = s['V{}_{}'.format(ia+1,ig)]<10.3
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


# In[90]:


def plot_channels(s,fp,ns=300,dpi=600):
    'Run through each channel to plot the multigains'
    for i in xrange(9):
        fig,ax = plot_channel_multigain(s,ch=i+1,ns=ns)
        u = 'V{}_0'.format(i+1)
        fig.savefig(fp+'{s.instname}_{s.daystr}_rad_{wvl}_allch.png'.format(s=s,wvl=s[u].name.split()[0]),
                    dpi=dpi,transparent=True)


# In[152]:


def plot_all_gainratios(s,fpp,ns=300,dpi=300):
    'Run through and plot each wavelenght and gain stage ratios'
    for i in xrange(9):
        fig,ax = plot_gains(s,i+1,ns=ns)
        u = 'V{}_0'.format(i+1)
        fig.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_{wvl}_gainratios.png'.format(s=s,wvl=s[u].name.split()[0]),
                    dpi=dpi,transparent=True)


# In[109]:


def clean_up_and_prep_for_plots(s,fpp):
    'to clean up the label and make the path for plotting, and other pre-plotting setups'
    if len(s.label.split('_'))>1:
        s.label = s.label.split('_')[-1]
    if not os.path.isdir(fpp+'{s.daystr}_{s.label}/'.format(s=s)):
        os.makedirs(fpp+'{s.daystr}_{s.label}/'.format(s=s))


# In[93]:


def plt_hskp(s,fpp,dpi=300):
    'plot all the housekeeping files and save them'
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


# In[107]:


def plt_gains(s,fpp,dpi=300):
    'plot all the different channels seperated by the gains'
    fig1,ax1 = plot_v(s,gain=0)
    ax1[0].autoscale(axis='y')
    ax1[1].autoscale(axis='y')
    fig1.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=dpi,transparent=True)

    fig2,ax2 = plot_v(s,gain=1)
    ax2[0].autoscale(axis='y')
    ax2[1].autoscale(axis='y')
    fig2.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=dpi,transparent=True)
    
    fig3,ax3 = plot_v(s,gain=2)
    ax3[0].autoscale(axis='y')
    ax3[1].autoscale(axis='y')
    fig3.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=dpi,transparent=True)
    
    fig4,ax4 = plot_v_expect(s,gain=0)
    #plt.title('5STARG radiometers 1x 2020-11-17')
    ax4[0].autoscale(axis='y')
    ax4[1].autoscale(axis='y')
    fig4.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=dpi,transparent=True)
    
    fig5,ax5 = plot_v_expect(s,gain=1)
    #plt.title('5STARG radiometers 1x 2020-11-17')
    ax5[0].autoscale(axis='y')
    ax5[1].autoscale(axis='y')
    fig5.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=dpi,transparent=True)
    
    #return [fig1,fig2,fig3,fig4,fig5],[ax1,ax2,ax3,ax4,ax5]


# In[216]:


def plt_corrs(s,fpp,dpi=300):
    'Run through the different gains and temperatures'
    for ig in [0,2,4]:
        for it in ['nir_block','vis_block','nir_tec_V']:
            fig = plot_corr_temp(s,ig,it)
            fig.savefig(fpp+'{s.daystr}_{s.label}/{s.instname}_{s.daystr}_corr_gain{}_temp{}.png'.format(ig,it,s=s),
                         dpi=dpi,transparent=True)


# # Load files

# ## Test 2020-10-16

# In[4]:


s = pd.read_csv(fp+'data/Data_2020_10_16_22_55_43_v2.txt',header=1)


# In[5]:


s['hour_UTC'][0]


# In[6]:


s['UTC'] = pd.to_datetime(s['hour_UTC'],format='%H:%M:%S.%f')
s['UTC'] = s['UTC'] + pd.offsets.DateOffset(years=s['Year'][0]-s['UTC'][0].year,
                                            months=s['month'][0]-s['UTC'][0].month,
                                            days=s['day'][0]-s['UTC'][0].day)


# ### Calculate the sun position

# In[14]:


import Sun_utils as su
import map_utils as mu


# In[32]:


lat = 36.970087
lon = -122.052537
alt = 74.0


# In[ ]:


sza,azi,sf,dec = mu.get_sza_azi(lat,lon,s['UTC'],alt=alt,return_sunf_and_dec=True)


# ## Test 2020-11-04

# In[75]:


s1 = pd.read_csv(fp+'data/5STARG_20201104_205500_RADIOMETERS_SweepTest.dat',header=5)
s2 = pd.read_csv(fp+'data/5STARG_20201104_213654_RADIOMETERS_MoveOffTest.dat',header=5)
s3 = pd.read_csv(fp+'data/5STARG_20201104_223408_RADIOMETERS_OutdoorDarkTest.dat',header=5)
s4 = pd.read_csv(fp+'data/5STARG_20201104_230010_RADIOMETERS_IndoorDarkTest.dat',header=5)


# In[76]:


s1.keys()


# In[78]:


s1['CH2_2'].name = 'something new'


# In[30]:


s1['UTC'] = pd.to_datetime(s1['hh:mm:ss'],format='%H:%M:%S.%f')
s2['UTC'] = pd.to_datetime(s2['hh:mm:ss'],format='%H:%M:%S.%f')
s3['UTC'] = pd.to_datetime(s3['hh:mm:ss'],format='%H:%M:%S.%f')
s4['UTC'] = pd.to_datetime(s4['hh:mm:ss'],format='%H:%M:%S.%f')


# In[31]:


s1['UTC'] = s1['UTC'] + pd.offsets.DateOffset(years=s1['YYYY'][0]-s1['UTC'][0].year,
                                            months=s1['MM'][0]-s1['UTC'][0].month,
                                            days=s1['DD'][0]-s1['UTC'][0].day)
s2['UTC'] = s2['UTC'] + pd.offsets.DateOffset(years=s2['YYYY'][0]-s2['UTC'][0].year,
                                            months=s2['MM'][0]-s2['UTC'][0].month,
                                            days=s2['DD'][0]-s2['UTC'][0].day)
s3['UTC'] = s3['UTC'] + pd.offsets.DateOffset(years=s3['YYYY'][0]-s3['UTC'][0].year,
                                            months=s3['MM'][0]-s3['UTC'][0].month,
                                            days=s3['DD'][0]-s3['UTC'][0].day)
s4['UTC'] = s4['UTC'] + pd.offsets.DateOffset(years=s4['YYYY'][0]-s4['UTC'][0].year,
                                            months=s4['MM'][0]-s4['UTC'][0].month,
                                            days=s4['DD'][0]-s4['UTC'][0].day)


# In[38]:


Titles = ['Sweep Test', 'Move Off Test','Dark Test Outside', 'Dark Test Inside']


# In[40]:


ss = [s1,s2,s3,s4]


# ### Calculate the sun position

# In[33]:


sza1,azi1,sf,dec = mu.get_sza_azi(lat,lon,s1['UTC'],alt=alt,return_sunf_and_dec=True)
sza2,azi2,sf,dec = mu.get_sza_azi(lat,lon,s2['UTC'],alt=alt,return_sunf_and_dec=True)
sza3,azi3,sf,dec = mu.get_sza_azi(lat,lon,s3['UTC'],alt=alt,return_sunf_and_dec=True)
sza4,azi4,sf,dec = mu.get_sza_azi(lat,lon,s4['UTC'],alt=alt,return_sunf_and_dec=True)


# In[42]:


ss[0]['sza'] = sza1
ss[1]['sza'] = sza2
ss[2]['sza'] = sza3
ss[3]['sza'] = sza4


# In[43]:


ss[0]['azi'] = azi1
ss[1]['azi'] = azi2
ss[2]['azi'] = azi3
ss[3]['azi'] = azi4


# ### Quantify the ref values

# In[44]:


ss[0]['ref'] = ss[0]['-5v_0_66']*3
ss[1]['ref'] = ss[1]['-5v_0_66']*3
ss[2]['ref'] = ss[2]['-5v_0_66']*3
ss[3]['ref'] = ss[3]['-5v_0_66']*3


# ## Test 2020-11-17

# In[97]:


fp


# In[108]:


s = reader_RADIOMETER(fp+'data/5STARG_20201117_001420_RADIOMETERS.dat')


# In[ ]:





# ## Test 2021-01-14

# In[55]:


fp


# In[56]:


os.listdir(fp+'data/')


# In[120]:


s = reader_RADIOMETER(fp+'data/5STARG_20210114_125517_RADIOMETERS.dat')


# ## Test 2021-01-25

# In[11]:


fp


# In[12]:


os.listdir(fp+'data/')


# In[64]:


s = reader_RADIOMETER(fp+'data/5STARG_20210125_154600_RADIOMETERS.dat')


# In[58]:


s.keys()


# ## Test 2021-01-26

# In[96]:


os.listdir(fp+'data/')


# In[87]:


s = reader_RADIOMETER(fp+'data/5STARG_20210126_103658_RADIOMETERS.dat')


# In[88]:


s.keys()


# In[89]:


s.daystr


# In[90]:


s.label


# In[97]:


s2 = reader_RADIOMETER(fp+'data/5STARG_20210126_143946_RADIOMETERS_HeatTest.dat')
s3 = reader_RADIOMETER(fp+'data/5STARG_20210126_154223_RADIOMETERS_IncandescentTest.dat')


# In[98]:


s2.label, s3.label


# ## Test 2021-01-27

# In[139]:


os.listdir(fp+'data/')


# In[130]:


s = reader_RADIOMETER(fp+'data/5STARG_20210127_141658_RADIOMETERS_SpeakerWire.dat')


# ## Test 2012-01-28

# In[141]:


s = reader_RADIOMETER(fp+'data/5STARG_20210128_155528_RADIOMETERS_TECDriverVoltages.dat')


# In[142]:


s.keys()


# In[144]:


s['nir_block'] = s['nir_block']*100.0
s['nir_board'] = s['nir_board']*100.0
s['vis_block'] = s['vis_block']*100.0
s['vis_board'] = s['vis_board']*100.0
s['mezz'] = s['mezz']*100.0


# In[145]:


s['nir_tec_V'] = s['nir_1020_tecb'] - s['nir_1020_teca']


# ## Test 2021-01-29 - on/off cover

# In[11]:


os.listdir(fp+'data/')


# In[13]:


s = reader_RADIOMETER(fp+'data/5STARG_20210129_163011_RADIOMETERS_CoverONOFF.dat')


# In[14]:


s.keys()


# In[15]:


s['nir_tec_V'] = s['nir_1020_tecb'] - s['nir_1020_teca']


# In[18]:


s['nir_block'] = s['nir_block_temp']*100.0
s['nir_board'] = s['nir_board_temp']*100.0
s['vis_block'] = s['vis_block_temp']*100.0
s['vis_board'] = s['vis_board_temp']*100.0
s['mezz'] = s['mezz']*100.0


# ## Test 2021-02-01 - 6 hour test

# In[99]:


os.listdir(fp+'data/')


# In[100]:


s = reader_RADIOMETER(fp+'data/5STARG_20210201_120748_RADIOMETERS_6HOUR.dat')


# In[101]:


s.keys()


# # Plot out data

# In[12]:


import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%H:%M:%S')


# In[13]:


s.keys()


# ## Housekeeping

# ### TEst 2020-10-16

# In[23]:


plt.figure()
plt.plot(s['UTC'],s[' Spare'])
plt.plot(s['UTC'],s[' 5vRef_Temp'])
plt.plot(s['UTC'],s[' -5v_0_66']*3,label='( -5v_0_66)*3')
plt.plot(s['UTC'],s[' 12POS_0_33'])
plt.plot(s['UTC'],s[' 12NEG_0_33'])
plt.legend()
plt.xlabel('UTC [Hour]')
plt.ylabel('Voltage [V]')
plt.gca().xaxis.set_major_formatter(myFmt)


# In[18]:


ref5v = s[' -5v_0_66']*3


# In[19]:


np.nanmean(ref5v)


# ### Test 2020-11-04

# In[36]:


s1.keys()


# In[50]:


for i,si in enumerate(ss):
    
    plt.figure()
    plt.plot(si['UTC'],si['Spare'])
    plt.plot(si['UTC'],si['5vRef_Temp'])
    plt.plot(si['UTC'],si['ref'],label='(-5v_0_66)*3:\n {:02.5f}+/-{:02.5f}'.format(np.nanmean(si['ref']),np.nanstd(si['ref'])))
    plt.plot(si['UTC'],si['12POS_0_33'])
    plt.plot(si['UTC'],si['12NEG_0_33'])
    plt.legend()
    plt.xlabel('UTC [Hour]')
    plt.ylabel('Voltage [V]')
    plt.gca().xaxis.set_major_formatter(myFmt)
    plt.title(Titles[i])


# ### Test 2021-01-14

# In[ ]:





# ## Against reference signal

# ### Test 2020-10-16

# #### 1x

# In[27]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],-s['CH1_0']-ref5v,label='340 nm')
ax[0].plot(s['UTC'],-s[' CH2_0']-ref5v,label='440 nm')
ax[0].plot(s['UTC'],-s[' CH3_0']-ref5v,label='675 nm')
ax[0].plot(s['UTC'],-s[' CH4_0']-ref5v,label='870 nm')
ax[0].plot(s['UTC'],-s[' CH5_0']-ref5v,label='1020 nm')
ax[0].set_title('VIS rad 1x Channels')
#ax[0].set_xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
ax[0].grid()
ax[0].legend()

ax[1].plot(s['UTC'],-s[' CH6_0']-ref5v,label='1020 nm')
ax[1].plot(s['UTC'],-s[' CH7_0']-ref5v,label='1240 nm')
ax[1].plot(s['UTC'],-s[' CH8_0']-ref5v,label='1640 nm')
ax[1].plot(s['UTC'],-s[' CH9_0']-ref5v,label='2200 nm')
ax[1].set_title('NIR rad 1x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_rad_test_1x_ref.png',dpi=600,transparent=True)


# In[22]:


fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(sza,-s['CH1_0']-ref5v,label='340 nm')
ax[0].plot(sza,-s[' CH2_0']-ref5v,label='440 nm')
ax[0].plot(sza,-s[' CH3_0']-ref5v,label='675 nm')
ax[0].plot(sza,-s[' CH4_0']-ref5v,label='870 nm')
ax[0].plot(sza,-s[' CH5_0']-ref5v,label='1020 nm')
ax[0].set_title('VIS rad 1x Channels')
#ax[0].set_xlabel('UTC [Hour]')
#ax[0].xaxis.set_major_formatter(myFmt)
#ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
ax[0].grid()
ax[0].legend()

ax[1].plot(sza,-s[' CH6_0']-ref5v,label='1020 nm')
ax[1].plot(sza,-s[' CH7_0']-ref5v,label='1240 nm')
ax[1].plot(sza,-s[' CH8_0']-ref5v,label='1640 nm')
ax[1].plot(sza,-s[' CH9_0']-ref5v,label='2200 nm')
ax[1].set_title('NIR rad 1x Channels')
ax[1].set_xlabel('Solar Zenith Angle')
#ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_SZA_rad_test_1x_ref.png',dpi=600,transparent=True)


# #### 100x

# In[29]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],-s[' CH1_2']-ref5v,label='340 nm')
ax[0].plot(s['UTC'],-s[' CH2_2']-ref5v,label='440 nm')
ax[0].plot(s['UTC'],-s[' CH3_2']-ref5v,label='675 nm')
ax[0].plot(s['UTC'],-s[' CH4_2']-ref5v,label='870 nm')
ax[0].plot(s['UTC'],-s[' CH5_2']-ref5v,label='1020 nm')
ax[0].set_title('VIS rad 100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V]')
ax[0].grid()
ax[0].legend()

ax[1].plot(s['UTC'],-s[' CH6_2']-ref5v,label='1020 nm')
ax[1].plot(s['UTC'],-s[' CH7_2']-ref5v,label='1240 nm')
ax[1].plot(s['UTC'],-s[' CH8_2']-ref5v,label='1640 nm')
ax[1].plot(s['UTC'],-s[' CH9_2']-ref5v,label='2200 nm')
ax[1].set_title('NIR rad 100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_rad_test_100x_ref.png',dpi=600,transparent=True)


# In[38]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],-s[' CH1_2']-ref5v,label='340 nm')
ax[0].plot(s['UTC'],-s[' CH2_2']-ref5v,label='440 nm')
ax[0].plot(s['UTC'],-s[' CH3_2']-ref5v,label='675 nm')
ax[0].plot(s['UTC'],-s[' CH4_2']-ref5v,label='870 nm')
ax[0].plot(s['UTC'],-s[' CH5_2']-ref5v,label='1020 nm')

ax[0].plot(s['UTC'],(-s['CH1_0']-ref5v)*100.0,label='Expected 340 nm',ls=':',alpha=0.4,c='tab:blue')
ax[0].plot(s['UTC'],(-s[' CH2_0']-ref5v)*100.0,label='Expected 440 nm',ls=':',alpha=0.4,c='tab:orange')
ax[0].plot(s['UTC'],(-s[' CH3_0']-ref5v)*100.0,label='Expected 675 nm',ls=':',alpha=0.4,c='tab:green')
ax[0].plot(s['UTC'],(-s[' CH4_0']-ref5v)*100.0,label='Expected 870 nm',ls=':',alpha=0.4,c='tab:red')
ax[0].plot(s['UTC'],(-s[' CH5_0']-ref5v)*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:purple')


ax[0].set_title('VIS rad 100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V]')
ax[0].grid()
ax[0].legend()
ax[0].set_ylim(-0.5,10.5)

ax[1].plot(s['UTC'],-s[' CH6_2']-ref5v,label='1020 nm')
ax[1].plot(s['UTC'],-s[' CH7_2']-ref5v,label='1240 nm')
ax[1].plot(s['UTC'],-s[' CH8_2']-ref5v,label='1640 nm')
ax[1].plot(s['UTC'],-s[' CH9_2']-ref5v,label='2200 nm')

ax[1].plot(s['UTC'],(-s[' CH6_0']-ref5v)*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:blue')
ax[1].plot(s['UTC'],(-s[' CH7_0']-ref5v)*100.0,label='Expected 1240 nm',ls=':',alpha=0.4,c='tab:orange')
ax[1].plot(s['UTC'],(-s[' CH8_0']-ref5v)*100.0,label='Expected 1640 nm',ls=':',alpha=0.4,c='tab:green')
ax[1].plot(s['UTC'],(-s[' CH9_0']-ref5v)*100.0,label='Expected 2200 nm',ls=':',alpha=0.4,c='tab:red')
ax[1].set_title('NIR rad 100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()
ax[1].set_ylim(-0.5,10.5)

plt.savefig(fp+'2020101622_rad_test_100x_ref_expect.png',dpi=600,transparent=True)


# In[64]:


fig,ax = plt.subplots(2,1)
ig1 = (-s[' CH1_2']-ref5v)<10.1
ig2 = (-s[' CH2_2']-ref5v)<10.1
ig3 = (-s[' CH3_2']-ref5v)<10.1
ig4 = (-s[' CH4_2']-ref5v)<10.1
ig5 = (-s[' CH5_2']-ref5v)<10.1

g1 = np.nanmean(((-s[' CH1_2']-ref5v)/(-s['CH1_0']-ref5v))[ig1])
g2 = np.nanmean(((-s[' CH2_2']-ref5v)/(-s[' CH2_0']-ref5v))[ig2])
g3 = np.nanmean(((-s[' CH3_2']-ref5v)/(-s[' CH3_0']-ref5v))[ig3])
g4 = np.nanmean(((-s[' CH4_2']-ref5v)/(-s[' CH4_0']-ref5v))[ig4])
g5 = np.nanmean(((-s[' CH5_2']-ref5v)/(-s[' CH5_0']-ref5v))[ig5])

eg1 = np.nanstd(((-s[' CH1_2']-ref5v)/(-s['CH1_0']-ref5v))[ig1])
eg2 = np.nanstd(((-s[' CH2_2']-ref5v)/(-s[' CH2_0']-ref5v))[ig2])
eg3 = np.nanstd(((-s[' CH3_2']-ref5v)/(-s[' CH3_0']-ref5v))[ig3])
eg4 = np.nanstd(((-s[' CH4_2']-ref5v)/(-s[' CH4_0']-ref5v))[ig4])
eg5 = np.nanstd(((-s[' CH5_2']-ref5v)/(-s[' CH5_0']-ref5v))[ig5])

ax[0].plot(s['UTC'][ig1],((-s[' CH1_2']-ref5v)/(-s['CH1_0']-ref5v))[ig1],'.',
           label='340 nm :{:3.1f}+/-{:3.1f}'.format(g1,eg1),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig2],((-s[' CH2_2']-ref5v)/(-s[' CH2_0']-ref5v))[ig2],'.',
           label='440 nm :{:3.1f}+/-{:3.1f}'.format(g2,eg2),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig3],((-s[' CH3_2']-ref5v)/(-s[' CH3_0']-ref5v))[ig3],'.',
           label='675 nm :{:3.1f}+/-{:3.1f}'.format(g3,eg3),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig4],((-s[' CH4_2']-ref5v)/(-s[' CH4_0']-ref5v))[ig4],'.',
           label='870 nm :{:3.1f}+/-{:3.1f}'.format(g4,eg4),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig5],((-s[' CH5_2']-ref5v)/(-s[' CH5_0']-ref5v))[ig5],'.',
           label='1020 nm :{:3.1f}+/-{:3.1f}'.format(g5,eg5),alpha=0.4,mec='None')

ax[0].set_title('Gain of VIS 100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Gain factor')
ax[0].grid()
ax[0].legend()
ax[0].set_ylim(0,500)

ign1 = ((-s[' CH6_2']-ref5v)<10.1) & ((-s[' CH6_2']-ref5v)>-0.1)
ign2 = ((-s[' CH7_2']-ref5v)<10.1) & ((-s[' CH7_2']-ref5v)>-0.1)
ign3 = ((-s[' CH8_2']-ref5v)<10.1) & ((-s[' CH8_2']-ref5v)>-0.1)
ign4 = ((-s[' CH9_2']-ref5v)<10.1) & ((-s[' CH9_2']-ref5v)>-0.1)

gn1 = np.nanmean(((-s[' CH6_2']-ref5v)/(-s[' CH6_0']-ref5v))[ign1])
gn2 = np.nanmean(((-s[' CH7_2']-ref5v)/(-s[' CH7_0']-ref5v))[ign2])
gn3 = np.nanmean(((-s[' CH8_2']-ref5v)/(-s[' CH8_0']-ref5v))[ign3])
gn4 = np.nanmean(((-s[' CH9_2']-ref5v)/(-s[' CH9_0']-ref5v))[ign4])

egn1 = np.nanstd(((-s[' CH6_2']-ref5v)/(-s[' CH6_0']-ref5v))[ign1])
egn2 = np.nanstd(((-s[' CH7_2']-ref5v)/(-s[' CH7_0']-ref5v))[ign2])
egn3 = np.nanstd(((-s[' CH8_2']-ref5v)/(-s[' CH8_0']-ref5v))[ign3])
egn4 = np.nanstd(((-s[' CH9_2']-ref5v)/(-s[' CH9_0']-ref5v))[ign4])

ax[1].plot(s['UTC'][ign1],((-s[' CH6_2']-ref5v)/(-s[' CH6_0']-ref5v))[ign1],'.',
           label='1020 nm:{:3.1f}+/-{:3.1f}'.format(gn1,egn1),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign2],((-s[' CH7_2']-ref5v)/(-s[' CH7_0']-ref5v))[ign2],'.',
           label='1240 nm:{:3.1f}+/-{:3.1f}'.format(gn2,egn2),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign3],((-s[' CH8_2']-ref5v)/(-s[' CH8_0']-ref5v))[ign3],'.',
           label='1640 nm:{:3.1f}+/-{:3.1f}'.format(gn3,egn3),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign4],((-s[' CH9_2']-ref5v)/(-s[' CH9_0']-ref5v))[ign4],'.',
           label='2200 nm:{:3.1f}+/-{:3.1f}'.format(gn4,egn4),alpha=0.4,mec='None')

ax[1].set_title('Gain of NIR 100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Gain factor')
ax[1].grid()
ax[1].legend()
ax[1].set_ylim(0,500)

plt.savefig(fp+'2020101622_rad_test_gainof100x.png',dpi=600,transparent=True)


# #### 100x - 100x

# In[66]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],-s[' CH1_4']-ref5v,label='340 nm')
ax[0].plot(s['UTC'],-s[' CH2_4']-ref5v,label='440 nm')
ax[0].plot(s['UTC'],-s[' CH3_4']-ref5v,label='675 nm')
ax[0].plot(s['UTC'],-s[' CH4_4']-ref5v,label='870 nm')
ax[0].plot(s['UTC'],-s[' CH5_4']-ref5v,label='1020 nm')
ax[0].set_title('VIS irrad 100x-100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
ax[0].grid()
ax[0].legend()
ax[0].set_ylim(-0.5,10.5)

ax[1].plot(s['UTC'],-s[' CH6_4']-ref5v,label='1020 nm')
ax[1].plot(s['UTC'],-s[' CH7_4']-ref5v,label='1240 nm')
ax[1].plot(s['UTC'],-s[' CH8_4']-ref5v,label='1640 nm')
ax[1].plot(s['UTC'],-s[' CH9_4']-ref5v,label='2200 nm')
ax[1].set_title('NIR irrad 100x-100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()
ax[1].set_ylim(-0.5,10.5)

plt.savefig(fp+'2020101622_rad_test_100x_100x_ref.png',dpi=600,transparent=True)


# In[63]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],-s[' CH1_4']-ref5v,label='340 nm')
ax[0].plot(s['UTC'],-s[' CH2_4']-ref5v,label='440 nm')
ax[0].plot(s['UTC'],-s[' CH3_4']-ref5v,label='675 nm')
ax[0].plot(s['UTC'],-s[' CH4_4']-ref5v,label='870 nm')
ax[0].plot(s['UTC'],-s[' CH5_4']-ref5v,label='1020 nm')

ax[0].plot(s['UTC'],(-s[' CH1_2']-ref5v)*100.0,label='Expected 340 nm',ls=':',alpha=0.4,c='tab:blue')
ax[0].plot(s['UTC'],(-s[' CH2_2']-ref5v)*100.0,label='Expected 440 nm',ls=':',alpha=0.4,c='tab:orange')
ax[0].plot(s['UTC'],(-s[' CH3_2']-ref5v)*100.0,label='Expected 675 nm',ls=':',alpha=0.4,c='tab:green')
ax[0].plot(s['UTC'],(-s[' CH4_2']-ref5v)*100.0,label='Expected 870 nm',ls=':',alpha=0.4,c='tab:red')
ax[0].plot(s['UTC'],(-s[' CH5_2']-ref5v)*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:purple')

ax[0].set_title('VIS irrad 100x-100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
ax[0].grid()
ax[0].legend(loc=2)
ax[0].set_ylim(-0.5,10.5)

ax[1].plot(s['UTC'],-s[' CH6_4']-ref5v,label='1020 nm')
ax[1].plot(s['UTC'],-s[' CH7_4']-ref5v,label='1240 nm')
ax[1].plot(s['UTC'],-s[' CH8_4']-ref5v,label='1640 nm')
ax[1].plot(s['UTC'],-s[' CH9_4']-ref5v,label='2200 nm')

ax[1].plot(s['UTC'],(-s[' CH6_2']-ref5v)*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:blue')
ax[1].plot(s['UTC'],(-s[' CH7_2']-ref5v)*100.0,label='Expected 1240 nm',ls=':',alpha=0.4,c='tab:orange')
ax[1].plot(s['UTC'],(-s[' CH8_2']-ref5v)*100.0,label='Expected 1640 nm',ls=':',alpha=0.4,c='tab:green')
ax[1].plot(s['UTC'],(-s[' CH9_2']-ref5v)*100.0,label='Expected 2200 nm',ls=':',alpha=0.4,c='tab:red')

ax[1].set_title('NIR irrad 100x-100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Irradiance [V]')
ax[1].grid()
ax[1].legend()
ax[1].set_ylim(-0.5,10.5)

plt.savefig(fp+'2020101622_rad_test_100x_100x_ref_expect.png',dpi=600,transparent=True)


# In[67]:


fig,ax = plt.subplots(2,1)
ig1 = (-s[' CH1_4']-ref5v)<10.3
ig2 = (-s[' CH2_4']-ref5v)<10.3
ig3 = (-s[' CH3_4']-ref5v)<10.3
ig4 = (-s[' CH4_4']-ref5v)<10.3
ig5 = (-s[' CH5_4']-ref5v)<10.3

g1 = np.nanmean(((-s[' CH1_4']-ref5v)/(-s[' CH1_2']-ref5v))[ig1])
g2 = np.nanmean(((-s[' CH2_4']-ref5v)/(-s[' CH2_2']-ref5v))[ig2])
g3 = np.nanmean(((-s[' CH3_4']-ref5v)/(-s[' CH3_2']-ref5v))[ig3])
g4 = np.nanmean(((-s[' CH4_4']-ref5v)/(-s[' CH4_2']-ref5v))[ig4])
g5 = np.nanmean(((-s[' CH5_4']-ref5v)/(-s[' CH5_2']-ref5v))[ig5])

eg1 = np.nanstd(((-s[' CH1_4']-ref5v)/(-s[' CH1_2']-ref5v))[ig1])
eg2 = np.nanstd(((-s[' CH2_4']-ref5v)/(-s[' CH2_2']-ref5v))[ig2])
eg3 = np.nanstd(((-s[' CH3_4']-ref5v)/(-s[' CH3_2']-ref5v))[ig3])
eg4 = np.nanstd(((-s[' CH4_4']-ref5v)/(-s[' CH4_2']-ref5v))[ig4])
eg5 = np.nanstd(((-s[' CH5_4']-ref5v)/(-s[' CH5_2']-ref5v))[ig5])

ax[0].plot(s['UTC'][ig1],((-s[' CH1_4']-ref5v)/(-s[' CH1_2']-ref5v))[ig1],'.',
           label='340 nm :{:3.1f}+/-{:3.1f}'.format(g1,eg1),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig2],((-s[' CH2_4']-ref5v)/(-s[' CH2_2']-ref5v))[ig2],'.',
           label='440 nm :{:3.1f}+/-{:3.1f}'.format(g2,eg2),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig3],((-s[' CH3_4']-ref5v)/(-s[' CH3_2']-ref5v))[ig3],'.',
           label='675 nm :{:3.1f}+/-{:3.1f}'.format(g3,eg3),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig4],((-s[' CH4_4']-ref5v)/(-s[' CH4_2']-ref5v))[ig4],'.',
           label='870 nm :{:3.1f}+/-{:3.1f}'.format(g4,eg4),alpha=0.4,mec='None')
ax[0].plot(s['UTC'][ig5],((-s[' CH5_4']-ref5v)/(-s[' CH5_2']-ref5v))[ig5],'.',
           label='1020 nm :{:3.1f}+/-{:3.1f}'.format(g5,eg5),alpha=0.4,mec='None')

ax[0].set_title('Gain of VIS 100x-100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Gain factor')
ax[0].grid()
ax[0].legend()
ax[0].set_ylim(0,500)

ign1 = ((-s[' CH6_4']-ref5v)<10.3) & ((-s[' CH6_4']-ref5v)>-0.1)
ign2 = ((-s[' CH7_4']-ref5v)<10.3) & ((-s[' CH7_4']-ref5v)>-0.1)
ign3 = ((-s[' CH8_4']-ref5v)<10.3) & ((-s[' CH8_4']-ref5v)>-0.1)
ign4 = ((-s[' CH9_4']-ref5v)<10.3) & ((-s[' CH9_4']-ref5v)>-0.1)

gn1 = np.nanmean(((-s[' CH6_4']-ref5v)/(-s[' CH6_2']-ref5v))[ign1])
gn2 = np.nanmean(((-s[' CH7_4']-ref5v)/(-s[' CH7_2']-ref5v))[ign2])
gn3 = np.nanmean(((-s[' CH8_4']-ref5v)/(-s[' CH8_2']-ref5v))[ign3])
gn4 = np.nanmean(((-s[' CH9_4']-ref5v)/(-s[' CH9_2']-ref5v))[ign4])

egn1 = np.nanstd(((-s[' CH6_4']-ref5v)/(-s[' CH6_2']-ref5v))[ign1])
egn2 = np.nanstd(((-s[' CH7_4']-ref5v)/(-s[' CH7_2']-ref5v))[ign2])
egn3 = np.nanstd(((-s[' CH8_4']-ref5v)/(-s[' CH8_2']-ref5v))[ign3])
egn4 = np.nanstd(((-s[' CH9_4']-ref5v)/(-s[' CH9_2']-ref5v))[ign4])

ax[1].plot(s['UTC'][ign1],((-s[' CH6_4']-ref5v)/(-s[' CH6_2']-ref5v))[ign1],'.',
           label='1020 nm:{:3.1f}+/-{:3.1f}'.format(gn1,egn1),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign2],((-s[' CH7_4']-ref5v)/(-s[' CH7_2']-ref5v))[ign2],'.',
           label='1240 nm:{:3.1f}+/-{:3.1f}'.format(gn2,egn2),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign3],((-s[' CH8_4']-ref5v)/(-s[' CH8_2']-ref5v))[ign3],'.',
           label='1640 nm:{:3.1f}+/-{:3.1f}'.format(gn3,egn3),alpha=0.4,mec='None')
ax[1].plot(s['UTC'][ign4],((-s[' CH9_4']-ref5v)/(-s[' CH9_2']-ref5v))[ign4],'.',
           label='2200 nm:{:3.1f}+/-{:3.1f}'.format(gn4,egn4),alpha=0.4,mec='None')

ax[1].set_title('Gain of NIR 100x-100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Gain factor')
ax[1].grid()
ax[1].legend()
ax[1].set_ylim(0,500)

plt.savefig(fp+'2020101622_rad_test_gainof100x_100x.png',dpi=600,transparent=True)


# ### Test 2020-11-04

# #### 1x

# In[53]:



for i,si in enumerate(ss):
    fig,ax = plt.subplots(2,1)
    ax[0].plot(si['UTC'],-si['CH1_0']-si['ref'],label='340 nm')
    ax[0].plot(si['UTC'],-si['CH2_0']-si['ref'],label='440 nm')
    ax[0].plot(si['UTC'],-si['CH3_0']-si['ref'],label='675 nm')
    ax[0].plot(si['UTC'],-si['CH4_0']-si['ref'],label='870 nm')
    ax[0].plot(si['UTC'],-si['CH5_0']-si['ref'],label='1020 nm')
    ax[0].set_title('VIS rad 1x Channels')
    #ax[0].set_xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(si['UTC'],-si['CH6_0']-si['ref'],label='1020 nm')
    ax[1].plot(si['UTC'],-si['CH7_0']-si['ref'],label='1240 nm')
    ax[1].plot(si['UTC'],-si['CH8_0']-si['ref'],label='1640 nm')
    ax[1].plot(si['UTC'],-si['CH9_0']-si['ref'],label='2200 nm')
    ax[1].set_title('NIR rad 1x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_1x_ref'+Titles[i]+'.png',dpi=600,transparent=True)


# #### 100x

# In[54]:


for i,si in enumerate(ss):
    fig,ax = plt.subplots(2,1)
    ax[0].plot(si['UTC'],-si['CH1_2']-si['ref'],label='340 nm')
    ax[0].plot(si['UTC'],-si['CH2_2']-si['ref'],label='440 nm')
    ax[0].plot(si['UTC'],-si['CH3_2']-si['ref'],label='675 nm')
    ax[0].plot(si['UTC'],-si['CH4_2']-si['ref'],label='870 nm')
    ax[0].plot(si['UTC'],-si['CH5_2']-si['ref'],label='1020 nm')
    ax[0].set_title('VIS rad 100x Channels')
    #ax[0].xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V]')
    ax[0].grid()
    ax[0].legend()

    ax[1].plot(si['UTC'],-si['CH6_2']-si['ref'],label='1020 nm')
    ax[1].plot(si['UTC'],-si['CH7_2']-si['ref'],label='1240 nm')
    ax[1].plot(si['UTC'],-si['CH8_2']-si['ref'],label='1640 nm')
    ax[1].plot(si['UTC'],-si['CH9_2']-si['ref'],label='2200 nm')
    ax[1].set_title('NIR rad 100x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_100x_ref'+Titles[i]+'.png',dpi=600,transparent=True)


# In[57]:


for i,si, in enumerate(ss):
    fig,ax = plt.subplots(2,1)
    ax[0].plot(si['UTC'],-si['CH1_2']-si['ref'],label='340 nm')
    ax[0].plot(si['UTC'],-si['CH2_2']-si['ref'],label='440 nm')
    ax[0].plot(si['UTC'],-si['CH3_2']-si['ref'],label='675 nm')
    ax[0].plot(si['UTC'],-si['CH4_2']-si['ref'],label='870 nm')
    ax[0].plot(si['UTC'],-si['CH5_2']-si['ref'],label='1020 nm')

    ax[0].plot(si['UTC'],(-si['CH1_0']-si['ref'])*100.0,label='Expected 340 nm',ls=':',alpha=0.4,c='tab:blue')
    ax[0].plot(si['UTC'],(-si['CH2_0']-si['ref'])*100.0,label='Expected 440 nm',ls=':',alpha=0.4,c='tab:orange')
    ax[0].plot(si['UTC'],(-si['CH3_0']-si['ref'])*100.0,label='Expected 675 nm',ls=':',alpha=0.4,c='tab:green')
    ax[0].plot(si['UTC'],(-si['CH4_0']-si['ref'])*100.0,label='Expected 870 nm',ls=':',alpha=0.4,c='tab:red')
    ax[0].plot(si['UTC'],(-si['CH5_0']-si['ref'])*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:purple')


    ax[0].set_title('VIS rad 100x Channels')
    #ax[0].xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V]')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_ylim(-0.5,10.5)

    ax[1].plot(si['UTC'],-si['CH6_2']-si['ref'],label='1020 nm')
    ax[1].plot(si['UTC'],-si['CH7_2']-si['ref'],label='1240 nm')
    ax[1].plot(si['UTC'],-si['CH8_2']-si['ref'],label='1640 nm')
    ax[1].plot(si['UTC'],-si['CH9_2']-si['ref'],label='2200 nm')

    ax[1].plot(si['UTC'],(-si['CH6_0']-si['ref'])*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:blue')
    ax[1].plot(si['UTC'],(-si['CH7_0']-si['ref'])*100.0,label='Expected 1240 nm',ls=':',alpha=0.4,c='tab:orange')
    ax[1].plot(si['UTC'],(-si['CH8_0']-si['ref'])*100.0,label='Expected 1640 nm',ls=':',alpha=0.4,c='tab:green')
    ax[1].plot(si['UTC'],(-si['CH9_0']-si['ref'])*100.0,label='Expected 2200 nm',ls=':',alpha=0.4,c='tab:red')
    ax[1].set_title('NIR rad 100x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylim(-0.5,10.5)

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_100x_ref_expect'+Titles[i]+'.png',dpi=600,transparent=True)


# #### 100x - 100x

# In[ ]:


for i, si in enumerate(ss):
    fig,ax = plt.subplots(2,1)
    ax[0].plot(si['UTC'],-si['CH1_4']-si['ref'],label='340 nm')
    ax[0].plot(si['UTC'],-si['CH2_4']-si['ref'],label='440 nm')
    ax[0].plot(si['UTC'],-si['CH3_4']-si['ref'],label='675 nm')
    ax[0].plot(si['UTC'],-si['CH4_4']-si['ref'],label='870 nm')
    ax[0].plot(si['UTC'],-si['CH5_4']-si['ref'],label='1020 nm')
    ax[0].set_title('VIS irrad 100x-100x Channels')
    #ax[0].xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_ylim(-0.5,10.5)

    ax[1].plot(si['UTC'],-si['CH6_4']-si['ref'],label='1020 nm')
    ax[1].plot(si['UTC'],-si['CH7_4']-si['ref'],label='1240 nm')
    ax[1].plot(si['UTC'],-si['CH8_4']-si['ref'],label='1640 nm')
    ax[1].plot(si['UTC'],-si['CH9_4']-si['ref'],label='2200 nm')
    ax[1].set_title('NIR irrad 100x-100x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylim(-0.5,10.5)

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_100x_100x_'+Titles[i]+'ref.png',dpi=600,transparent=True)


# In[ ]:


for i,si in enumerate(ss):
    fig,ax = plt.subplots(2,1)
    ax[0].plot(si['UTC'],-si['CH1_4']-si['ref'],label='340 nm')
    ax[0].plot(si['UTC'],-si['CH2_4']-si['ref'],label='440 nm')
    ax[0].plot(si['UTC'],-si['CH3_4']-si['ref'],label='675 nm')
    ax[0].plot(si['UTC'],-si['CH4_4']-si['ref'],label='870 nm')
    ax[0].plot(si['UTC'],-si['CH5_4']-si['ref'],label='1020 nm')

    ax[0].plot(si['UTC'],(-si['CH1_2']-si['ref'])*100.0,label='Expected 340 nm',ls=':',alpha=0.4,c='tab:blue')
    ax[0].plot(si['UTC'],(-si['CH2_2']-si['ref'])*100.0,label='Expected 440 nm',ls=':',alpha=0.4,c='tab:orange')
    ax[0].plot(si['UTC'],(-si['CH3_2']-si['ref'])*100.0,label='Expected 675 nm',ls=':',alpha=0.4,c='tab:green')
    ax[0].plot(si['UTC'],(-si['CH4_2']-si['ref'])*100.0,label='Expected 870 nm',ls=':',alpha=0.4,c='tab:red')
    ax[0].plot(si['UTC'],(-si['CH5_2']-si['ref'])*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:purple')

    ax[0].set_title('VIS irrad 100x-100x Channels')
    #ax[0].xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Irradiance [V] (X+5V_ref)*-1')
    ax[0].grid()
    ax[0].legend(loc=2)
    ax[0].set_ylim(-0.5,10.5)

    ax[1].plot(si['UTC'],-si['CH6_4']-si['ref'],label='1020 nm')
    ax[1].plot(si['UTC'],-si['CH7_4']-si['ref'],label='1240 nm')
    ax[1].plot(si['UTC'],-si['CH8_4']-si['ref'],label='1640 nm')
    ax[1].plot(si['UTC'],-si['CH9_4']-si['ref'],label='2200 nm')

    ax[1].plot(si['UTC'],(-si['CH6_2']-si['ref'])*100.0,label='Expected 1020 nm',ls=':',alpha=0.4,c='tab:blue')
    ax[1].plot(si['UTC'],(-si['CH7_2']-si['ref'])*100.0,label='Expected 1240 nm',ls=':',alpha=0.4,c='tab:orange')
    ax[1].plot(si['UTC'],(-si['CH8_2']-si['ref'])*100.0,label='Expected 1640 nm',ls=':',alpha=0.4,c='tab:green')
    ax[1].plot(si['UTC'],(-si['CH9_2']-si['ref'])*100.0,label='Expected 2200 nm',ls=':',alpha=0.4,c='tab:red')

    ax[1].set_title('NIR irrad 100x-100x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Irradiance [V]')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylim(-0.5,10.5)

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_100x_100x_ref_expect'+Titles[i]+'.png',dpi=600,transparent=True)


# In[ ]:


for i,si in enumerate(ss):

    fig,ax = plt.subplots(2,1)
    ig1 = (-si['CH1_4']-si['ref'])<10.3
    ig2 = (-si['CH2_4']-si['ref'])<10.3
    ig3 = (-si['CH3_4']-si['ref'])<10.3
    ig4 = (-si['CH4_4']-si['ref'])<10.3
    ig5 = (-si['CH5_4']-si['ref'])<10.3

    g1 = np.nanmean(((-si['CH1_4']-si['ref'])/(-si['CH1_2']-si['ref']))[ig1])
    g2 = np.nanmean(((-si['CH2_4']-si['ref'])/(-si['CH2_2']-si['ref']))[ig2])
    g3 = np.nanmean(((-si['CH3_4']-si['ref'])/(-si['CH3_2']-si['ref']))[ig3])
    g4 = np.nanmean(((-si['CH4_4']-si['ref'])/(-si['CH4_2']-si['ref']))[ig4])
    g5 = np.nanmean(((-si['CH5_4']-si['ref'])/(-si['CH5_2']-si['ref']))[ig5])

    eg1 = np.nanstd(((-si['CH1_4']-si['ref'])/(-si['CH1_2']-si['ref']))[ig1])
    eg2 = np.nanstd(((-si['CH2_4']-si['ref'])/(-si['CH2_2']-si['ref']))[ig2])
    eg3 = np.nanstd(((-si['CH3_4']-si['ref'])/(-si['CH3_2']-si['ref']))[ig3])
    eg4 = np.nanstd(((-si['CH4_4']-si['ref'])/(-si['CH4_2']-si['ref']))[ig4])
    eg5 = np.nanstd(((-si['CH5_4']-si['ref'])/(-si['CH5_2']-si['ref']))[ig5])

    ax[0].plot(si['UTC'][ig1],((-si['CH1_4']-si['ref'])/(-si['CH1_2']-si['ref']))[ig1],'.',
               label='340 nm :{:3.1f}+/-{:3.1f}'.format(g1,eg1),alpha=0.4,mec='None')
    ax[0].plot(si['UTC'][ig2],((-si['CH2_4']-si['ref'])/(-si['CH2_2']-si['ref']))[ig2],'.',
               label='440 nm :{:3.1f}+/-{:3.1f}'.format(g2,eg2),alpha=0.4,mec='None')
    ax[0].plot(si['UTC'][ig3],((-si['CH3_4']-si['ref'])/(-si['CH3_2']-si['ref']))[ig3],'.',
               label='675 nm :{:3.1f}+/-{:3.1f}'.format(g3,eg3),alpha=0.4,mec='None')
    ax[0].plot(si['UTC'][ig4],((-si['CH4_4']-si['ref'])/(-si['CH4_2']-si['ref']))[ig4],'.',
               label='870 nm :{:3.1f}+/-{:3.1f}'.format(g4,eg4),alpha=0.4,mec='None')
    ax[0].plot(si['UTC'][ig5],((-si['CH5_4']-si['ref'])/(-si['CH5_2']-si['ref']))[ig5],'.',
               label='1020 nm :{:3.1f}+/-{:3.1f}'.format(g5,eg5),alpha=0.4,mec='None')

    ax[0].set_title('Gain of VIS 100x-100x Channels')
    #ax[0].xlabel('UTC [Hour]')
    ax[0].xaxis.set_major_formatter(myFmt)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Gain factor')
    ax[0].grid()
    ax[0].legend()
    ax[0].set_ylim(0,500)

    ign1 = ((-si['CH6_4']-si['ref'])<10.3) & ((-si['CH6_4']-si['ref'])>-0.1)
    ign2 = ((-si['CH7_4']-si['ref'])<10.3) & ((-si['CH7_4']-si['ref'])>-0.1)
    ign3 = ((-si['CH8_4']-si['ref'])<10.3) & ((-si['CH8_4']-si['ref'])>-0.1)
    ign4 = ((-si['CH9_4']-si['ref'])<10.3) & ((-si['CH9_4']-si['ref'])>-0.1)

    gn1 = np.nanmean(((-si['CH6_4']-si['ref'])/(-si['CH6_2']-si['ref']))[ign1])
    gn2 = np.nanmean(((-si['CH7_4']-si['ref'])/(-si['CH7_2']-si['ref']))[ign2])
    gn3 = np.nanmean(((-si['CH8_4']-si['ref'])/(-si['CH8_2']-si['ref']))[ign3])
    gn4 = np.nanmean(((-si['CH9_4']-si['ref'])/(-si['CH9_2']-si['ref']))[ign4])

    egn1 = np.nanstd(((-si['CH6_4']-si['ref'])/(-si['CH6_2']-si['ref']))[ign1])
    egn2 = np.nanstd(((-si['CH7_4']-si['ref'])/(-si['CH7_2']-si['ref']))[ign2])
    egn3 = np.nanstd(((-si['CH8_4']-si['ref'])/(-si['CH8_2']-si['ref']))[ign3])
    egn4 = np.nanstd(((-si['CH9_4']-si['ref'])/(-si['CH9_2']-si['ref']))[ign4])

    ax[1].plot(si['UTC'][ign1],((-si['CH6_4']-si['ref'])/(-si['CH6_2']-si['ref']))[ign1],'.',
               label='1020 nm:{:3.1f}+/-{:3.1f}'.format(gn1,egn1),alpha=0.4,mec='None')
    ax[1].plot(si['UTC'][ign2],((-si['CH7_4']-si['ref'])/(-si['CH7_2']-si['ref']))[ign2],'.',
               label='1240 nm:{:3.1f}+/-{:3.1f}'.format(gn2,egn2),alpha=0.4,mec='None')
    ax[1].plot(si['UTC'][ign3],((-si['CH8_4']-si['ref'])/(-si['CH8_2']-si['ref']))[ign3],'.',
               label='1640 nm:{:3.1f}+/-{:3.1f}'.format(gn3,egn3),alpha=0.4,mec='None')
    ax[1].plot(si['UTC'][ign4],((-si['CH9_4']-si['ref'])/(-si['CH9_2']-si['ref']))[ign4],'.',
               label='2200 nm:{:3.1f}+/-{:3.1f}'.format(gn4,egn4),alpha=0.4,mec='None')

    ax[1].set_title('Gain of NIR 100x-100x Channels')
    ax[1].set_xlabel('UTC [Hour]')
    ax[1].xaxis.set_major_formatter(myFmt)
    ax[1].set_ylabel('Gain factor')
    ax[1].grid()
    ax[1].legend()
    ax[1].set_ylim(0,500)

    plt.suptitle(Titles[i])
    plt.savefig(fp+'20201104_rad_test_gainof100x_100x'+Titles[i]+'.png',dpi=600,transparent=True)


# ## Just voltages

# In[38]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],s['CH1_0']*-1.0,label='340 nm')
ax[0].plot(s['UTC'],s[' CH2_0']*-1.0,label='440 nm')
ax[0].plot(s['UTC'],s[' CH3_0']*-1.0,label='675 nm')
ax[0].plot(s['UTC'],s[' CH4_0']*-1.0,label='870 nm')
ax[0].plot(s['UTC'],s[' CH5_0']*-1.0,label='1020 nm')
ax[0].set_title('VIS rad 1x Channels')
#ax[0].set_xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Radiance [V] *-1')
ax[0].grid()
ax[0].legend()

ax[1].plot(s['UTC'],s[' CH6_0']*-1.0,label='1020 nm')
ax[1].plot(s['UTC'],s[' CH7_0']*-1.0,label='1240 nm')
ax[1].plot(s['UTC'],s[' CH8_0']*-1.0,label='1640 nm')
ax[1].plot(s['UTC'],s[' CH9_0']*-1.0,label='2200 nm')
ax[1].set_title('NIR rad 1x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Radiance [V] *-1')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_rad_test_1x.png',dpi=600,transparent=True)


# In[36]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],s[' CH1_2']*-1.0,label='340 nm')
ax[0].plot(s['UTC'],s[' CH2_2']*-1.0,label='440 nm')
ax[0].plot(s['UTC'],s[' CH3_2']*-1.0,label='675 nm')
ax[0].plot(s['UTC'],s[' CH4_2']*-1.0,label='870 nm')
ax[0].plot(s['UTC'],s[' CH5_2']*-1.0,label='1020 nm')
ax[0].set_title('VIS rad 100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Radiance [V] *-1')
ax[0].grid()
ax[0].legend()

ax[1].plot(s['UTC'],s[' CH6_2']*-1.0,label='1020 nm')
ax[1].plot(s['UTC'],s[' CH7_2']*-1.0,label='1240 nm')
ax[1].plot(s['UTC'],s[' CH8_2']*-1.0,label='1640 nm')
ax[1].plot(s['UTC'],s[' CH9_2']*-1.0,label='2200 nm')
ax[1].set_title('NIR rad 100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Radiance [V] *-1')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_rad_test_100x.png',dpi=600,transparent=True)


# In[37]:


fig,ax = plt.subplots(2,1)
ax[0].plot(s['UTC'],s[' CH1_4']*-1.0,label='340 nm')
ax[0].plot(s['UTC'],s[' CH2_4']*-1.0,label='440 nm')
ax[0].plot(s['UTC'],s[' CH3_4']*-1.0,label='675 nm')
ax[0].plot(s['UTC'],s[' CH4_4']*-1.0,label='870 nm')
ax[0].plot(s['UTC'],s[' CH5_4']*-1.0,label='1020 nm')
ax[0].set_title('VIS rad 100x-100x Channels')
#ax[0].xlabel('UTC [Hour]')
ax[0].xaxis.set_major_formatter(myFmt)
ax[0].set_xticklabels([])
ax[0].set_ylabel('Radiance [V] *-1')
ax[0].grid()
ax[0].legend()

ax[1].plot(s['UTC'],s[' CH6_4']*-1.0,label='1020 nm')
ax[1].plot(s['UTC'],s[' CH7_4']*-1.0,label='1240 nm')
ax[1].plot(s['UTC'],s[' CH8_4']*-1.0,label='1640 nm')
ax[1].plot(s['UTC'],s[' CH9_4']*-1.0,label='2200 nm')
ax[1].set_title('NIR rad 100x-100x Channels')
ax[1].set_xlabel('UTC [Hour]')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].set_ylabel('Radiance [V] *-1')
ax[1].grid()
ax[1].legend()

plt.savefig(fp+'2020101622_rad_test_100x_100x.png',dpi=600,transparent=True)


# ## Plot out test 2020-11-17, with new functions

# In[176]:


fig = plot_housekeeping(s)
fig[0].gca().set_title('5STARG test: 2020-11-17')
fig[0].savefig(fp+'5STARG_20201117_Housekeeping.png',dpi=600,transparent=True)


# In[121]:


fig,ax = plot_v(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
fig.savefig(fp+'5STARG_20201117_rad_1x.png',dpi=600,transparent=True)


# In[122]:


fig,ax = plot_v(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
fig.savefig(fp+'5STARG_20201117_rad_100x.png',dpi=600,transparent=True)


# In[123]:


fig,ax = plot_v(s,gain=2)
#plt.title('5STARG radiometers 1x 2020-11-17')
fig.savefig(fp+'5STARG_20201117_rad_100x_100x.png',dpi=600,transparent=True)


# In[135]:


fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
fig.savefig(fp+'5STARG_20201117_rad_100x_expect.png',dpi=600,transparent=True)


# In[158]:


fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[1].set_ylim(-100,25)
#fig.savefig(fp+'5STARG_20201117_rad_100x_100x_expect.png',dpi=600,transparent=True)


# ## Plot out Test 2021-01-14

# In[125]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)


# In[131]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.015,0.02)
fig.savefig(fp+'{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)


# In[136]:


fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.015,0.8)
fig.savefig(fp+'{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)


# In[141]:


fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(-0.5,0.5)
ax[1].set_ylim(-0.015,10.5)
fig.savefig(fp+'{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)


# In[143]:


fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-1.0,2.0)
fig.savefig(fp+'{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[145]:


fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.5,1.5)
ax[1].set_ylim(0,11.0)
fig.savefig(fp+'{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# ### Checkout the 440 nm

# In[181]:


from Sp_parameters import smooth
myFmt = mdates.DateFormatter('%H:%M:%S')
fig,ax = plt.subplots(3,1,sharex=True,figsize=(7,7))

cs = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple']

u,u1,u2 = 'V2_0','V2_2','V2_4'
ns = 300
ax[0].plot(s['UTC'],s[u],'.',label='raw',c=cs[iu],alpha=0.02)
ax[0].plot(s['UTC'],smooth(s[u],ns,old=True),'-',label='smooth n={}'.format(ns),c='k')
ax[0].legend()
ax[1].plot(s['UTC'],s[u1],'.',label=s[u1].name,c=cs[iu],alpha=0.02)
ax[1].plot(s['UTC'],smooth(s[u1],ns,old=True),'-',label='smooth n={}'.format(ns),c='k')
ax[2].plot(s['UTC'],s[u2],'.',label=s[u1].name,c=cs[iu],alpha=0.02)
ax[2].plot(s['UTC'],smooth(s[u2],ns,old=True),'-',label='smooth n={}'.format(ns),c='k')

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

fig.savefig(fp+'{s.instname}_{s.daystr}_rad_440nm_allch.png'.format(s=s),dpi=600,transparent=True)


# In[ ]:


plot_channels(s,fp,dpi=200)


# ## Plot out test 2021-01-25

# In[30]:


s.daystr


# In[81]:


os.makedirs(fp+'plots/{s.daystr}/'.format(s=s))


# In[85]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=600,transparent=True)


# In[34]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.015,0.02)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)


# In[36]:


fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.015,0.02)
ax[1].set_ylim(-0.015,0.1)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)


# In[44]:


fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(0.4,1.1)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)


# In[45]:


fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-1.0,2.0)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[47]:


fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.5,2.5)
ax[1].set_ylim(-0.5,11.0)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[ ]:


plot_channels(s,fp+'plots/{s.daystr}/'.format(s=s),dpi=200)


# ## Plot out test 2021-01-26

# ### Battery test

# In[91]:


os.makedirs(fp+'plots/{s.daystr}/'.format(s=s))


# In[92]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=600,transparent=True)


# In[94]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.025,0.01)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.025,0.025)
ax[1].set_ylim(-0.25,0.1)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(-0.5,1.5)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-2.5,2.0)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-2.5,2.5)
ax[1].set_ylim(-1.5,11.0)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[95]:


plot_channels(s,fp+'plots/{s.daystr}/'.format(s=s),dpi=200)


# ### Heating test

# In[99]:


os.makedirs(fp+'plots/{s.daystr}_{s.label}/'.format(s=s2))


# In[127]:


fig = plot_housekeeping(s2)
fig[0].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_Housekeeping.png'.format(s=s2),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_Housekeeping_inputV.png'.format(s=s2),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_Housekeeping_Temps.png'.format(s=s2),dpi=600,transparent=True)


# In[128]:


fig,ax = plot_v(s2,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.015,0.015)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_rad_1x.png'.format(s=s2),dpi=600,transparent=True)
fig,ax = plot_v(s2,gain=1)
ax[0].set_ylim(-0.025,0.025)
ax[1].set_ylim(-0.0,0.55)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_rad_100x.png'.format(s=s2),dpi=600,transparent=True)
fig,ax = plot_v(s2,gain=2)
ax[0].set_ylim(-0.2,1.2)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_rad_100x_100x.png'.format(s=s2),dpi=600,transparent=True)
fig,ax = plot_v_expect(s2,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-1.0,1.5)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_rad_100x_expect.png'.format(s=s2),dpi=600,transparent=True)
fig,ax = plot_v_expect(s2,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.0,2.5)
ax[1].set_ylim(-0.5,11.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_{s.label}_rad_100x_100x_expect.png'.format(s=s2),dpi=600,transparent=True)


# In[106]:


plot_channels(s2,fp+'plots/{s.daystr}_{s.label}/'.format(s=s2),dpi=200)


# ### Incandescent light test

# In[118]:


s3.label = s3.label.split('_')[-1]


# In[119]:


os.makedirs(fp+'plots/{s.daystr}_{s.label}/'.format(s=s3))


# In[120]:


fig = plot_housekeeping(s3)
fig[0].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_Housekeeping.png'.format(s=s3),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_Housekeeping_inputV.png'.format(s=s3),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_Housekeeping_Temps.png'.format(s=s3),dpi=600,transparent=True)


# In[125]:


fig,ax = plot_v(s3,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.02,0.2)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_rad_1x.png'.format(s=s3),dpi=600,transparent=True)
fig,ax = plot_v(s3,gain=1)
ax[0].set_ylim(-0.015,2.0)
ax[1].set_ylim(-0.0,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_rad_100x.png'.format(s=s3),dpi=600,transparent=True)
fig,ax = plot_v(s3,gain=2)
ax[0].set_ylim(-0.5,10.5)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_rad_100x_100x.png'.format(s=s3),dpi=600,transparent=True)
fig,ax = plot_v_expect(s3,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-0.5,1.5)
ax[1].set_ylim(-2.5,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_rad_100x_expect.png'.format(s=s3),dpi=600,transparent=True)
fig,ax = plot_v_expect(s3,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.5,10.5)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}/{s.instname}_{s.daystr}_{s.label}_rad_100x_100x_expect.png'.format(s=s3),dpi=600,transparent=True)


# In[126]:


plot_channels(s3,fp+'plots/{s.daystr}_{s.label}/'.format(s=s3),dpi=200)


# ## Plot out test 2021-01-27

# In[133]:


s.label = s.label.split('_')[-1]


# In[134]:


os.makedirs(fp+'plots/{s.daystr}_{s.label}/'.format(s=s))


# In[135]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=600,transparent=True)


# In[ ]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.01,0.01)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.01,0.025)
ax[1].set_ylim(0.0,0.12)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(0.4,1.6)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-2.5,2.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-2.5,2.5)
ax[1].set_ylim(-1.5,11.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[138]:


plot_channels(s,fp+'plots/{s.daystr}_{s.label}/'.format(s=s),dpi=200)


# ## Plot out test 2021-01-28 - TEC voltages

# In[146]:


s.label = s.label.split('_')[-1]


# In[147]:


os.makedirs(fp+'plots/{s.daystr}_{s.label}/'.format(s=s))


# In[148]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=600,transparent=True)


# In[165]:


fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(s['UTC'],s['nir_1020_teca'],'.')
ax[0].plot(s['UTC'],s['nir_1020_tecb'],'.')
ax[1].plot(s['UTC'],s['nir_tec_V'],'-',label='tecb - teca',lw=0.2)
ax[0].legend()
ax[1].legend()
ax[1].set_xlabel('UTC [Hour]')
ax[0].set_ylabel('TEC Voltages [V]')
ax[1].set_ylabel('TEC Difference [V]')
ax[0].set_title('{s.instname} - {s.daystr} {s.label}- NIR 1020 nm TEC Voltages'.format(s=s))
ax[1].xaxis.set_major_locator(plt.MaxNLocator(7))
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%H:%M:%S')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].grid()
ax[0].grid()
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_NIR_TECvoltages.png'.format(s=s),dpi=600,transparent=True)


# In[161]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.01,0.015)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.01,0.025)
ax[1].set_ylim(0.0,0.22)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(0.4,1.6)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-2.5,2.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-0.5,2.5)
ax[1].set_ylim(-1.5,11.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[162]:


plot_channels(s,fp+'plots/{s.daystr}_{s.label}/'.format(s=s),dpi=200)


# ## Plot out test 2021-01-29 - on/off cover

# In[19]:


s.label = s.label.split('_')[-1]


# In[20]:


os.makedirs(fp+'plots/{s.daystr}_{s.label}/'.format(s=s))


# In[21]:


fig = plot_housekeeping(s)
fig[0].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping.png'.format(s=s),dpi=600,transparent=True)
fig[1].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_inputV.png'.format(s=s),dpi=600,transparent=True)
fig[2].savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_Housekeeping_Temps.png'.format(s=s),dpi=600,transparent=True)


# In[22]:


fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(s['UTC'],s['nir_1020_teca'],'.')
ax[0].plot(s['UTC'],s['nir_1020_tecb'],'.')
ax[1].plot(s['UTC'],s['nir_tec_V'],'-',label='tecb - teca',lw=0.2)
ax[0].legend()
ax[1].legend()
ax[1].set_xlabel('UTC [Hour]')
ax[0].set_ylabel('TEC Voltages [V]')
ax[1].set_ylabel('TEC Difference [V]')
ax[0].set_title('{s.instname} - {s.daystr} {s.label}- NIR 1020 nm TEC Voltages'.format(s=s))
ax[1].xaxis.set_major_locator(plt.MaxNLocator(7))
import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%H:%M:%S')
ax[1].xaxis.set_major_formatter(myFmt)
ax[1].grid()
ax[0].grid()
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_NIR_TECvoltages.png'.format(s=s),dpi=600,transparent=True)


# In[23]:


fig,ax = plot_v(s,gain=0)
ax[0].set_ylim(-0.015,0.015)
ax[1].set_ylim(-0.01,0.015)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_1x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=1)
ax[0].set_ylim(-0.01,0.025)
ax[1].set_ylim(0.0,0.22)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v(s,gain=2)
ax[0].set_ylim(0.4,1.6)
ax[1].set_ylim(-1.5,10.5)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=0)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-1.2,1.2)
ax[1].set_ylim(-2.5,2.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_expect.png'.format(s=s),dpi=600,transparent=True)
fig,ax = plot_v_expect(s,gain=1)
#plt.title('5STARG radiometers 1x 2020-11-17')
ax[0].set_ylim(-0.5,2.5)
ax[1].set_ylim(-1.5,11.0)
fig.savefig(fp+'plots/{s.daystr}_{s.label}/{s.instname}_{s.daystr}_rad_100x_100x_expect.png'.format(s=s),dpi=600,transparent=True)


# In[24]:


plot_channels(s,fp+'plots/{s.daystr}_{s.label}/'.format(s=s),dpi=200)


# ## Plot out test 2021-02-01 - 6 hour drift

# In[105]:


clean_up_and_prep_for_plots(s,fp+'plots/')


# In[106]:


plt_hskp(s,fp+'plots/')


# In[108]:


plt_gains(s,fp+'plots/')


# In[163]:


plot_all_gainratios(s,fp+'plots/')


# In[54]:


plot_channels(s,fp+'plots/{s.daystr}_{s.label}/'.format(s=s),dpi=200)


# In[224]:


plt_corrs(s,fp+'plots/')

