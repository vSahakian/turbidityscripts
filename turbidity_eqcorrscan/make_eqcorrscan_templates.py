#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 16:44:36 2022

@author: vjs
"""

## Make templates within eqcorrscan for 2019

import logging

from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.core.template_gen import template_gen

#from eqcorrscan.core import day_gen

import obspy as obs
from obspy.io.sac import SACTrace

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from obspy.core.event import Event, Origin, WaveformStreamID, Pick, Catalog


from glob import glob
import os

# %% 


## Input template mseed directory from csv:
input_mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/csv2mseed/'


## Template output mseed directory:
template_mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/template_mseed/'

## Template figure directory
template_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/template_mseed_plots/'

## Template processing parameters:
        ## Note: data are 1 min sampled, so 1/60 = 0.01667 Hz. highcut = fN = 0.008333.
        ## Want to keep things with at least an hour periodicity = 1/(60*60) Hz = 0.000277
        ##   so set lowcut to 0.0002
        ## Want to keep 40 min for waveform, so set length to be 40 min * 60 sec/min
filter_lowcut=None
filter_highcut=0.008333
template_samp_rate=1/60
template_filt_order=4  ## originally 4 in example...
template_length=60*60
template_swin='P'
event_prepick=0.05
data_pre_processed=True
data_all_horiz=False


# %%
# Get all the SAC-files associated with one event.
#sac_files = glob.glob(TEST_PATH + '/SAC/2014p611252/*')
template_mseed_list = glob(input_mseed_directory + '*.mseed')

for i_file in range(len(template_mseed_list)):
    ## Define file path:
    i_mseed_path = template_mseed_list[i_file]
    
    print('working on %s \n' % i_mseed_path)
    ## Basename:
    i_mseed_basename = os.path.basename(i_mseed_path).split('.mseed')[0]
    
    ## Read in stream object:
    i_stream = obs.core.read(i_mseed_path)
    
    ## Plot the data:
    i_fig = plt.figure()

    
    plt.plot(i_stream[0].times("utcdatetime"),i_stream[0].data,color='orange',label='Unfiltered')
    
    # Now we need to create an event!
    i_event = Event()
    i_event.origins.append(Origin())
    # print(st[0].stats.sac.keys())
    i_reference_time = obs.UTCDateTime(
        year=i_stream[0].stats.starttime.year, julday=i_stream[0].stats.starttime.julday,
        hour=i_stream[0].stats.starttime.hour, minute=i_stream[0].stats.starttime.minute,
        second=i_stream[0].stats.starttime.second,
        microsecond=i_stream[0].stats.starttime.microsecond * 1000)

    i_event.origins[0].time = i_reference_time

    ## No location info:
    i_event.origins[0].latitude = None
    i_event.origins[0].longitude = None
    i_event.origins[0].depth = None

    ## Add picks to it:
    i_pick_time = i_reference_time + 0.10 ## make the pick time be a little after hte start of the file
    i_phase_hint = 'P'
    
    i_waveform_id = WaveformStreamID(station_code=i_stream[0].stats.station,
                network_code=i_stream[0].stats.network,
                channel_code=i_stream[0].stats.channel)
    
    i_pick = Pick(waveform_id=i_waveform_id,
                phase_hint=i_phase_hint,
                time=i_pick_time)
    
    i_event.picks.append(i_pick)
    
    ## Make a catalog object:
    i_catalog = Catalog([i_event])
        
    ## Generate template:
    i_templates = template_gen(
        method='from_meta_file', meta_file=i_catalog, st=i_stream, lowcut=filter_lowcut, highcut=filter_highcut,
        samp_rate=template_samp_rate, filt_order=template_filt_order, length=template_length, 
        swin=template_swin, prepick=event_prepick, delayed=False,
        pre_processed=data_pre_processed,all_horiz=data_all_horiz)
    
    plt.plot(i_templates[0][0].times("utcdatetime"),i_templates[0][0].data,color='green',label='Template')

    ## Add legend to figure:
    plt.legend()
    ## And title
    i_title_statement = '%s'
    if filter_lowcut == None and filter_highcut != None:
        plt.title('%s, \n highcut %.3f Hz' % (i_mseed_basename,filter_highcut))
    elif filter_lowcut == None and filter_highcut == None:
        plt.title('%s' % (i_mseed_basename))
    else:
        plt.title('%s, \n lowcut %.3f Hz, highcut %.3f Hz' % (i_mseed_basename,filter_lowcut,filter_highcut))

    ## Set axis formats:
    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    
    ## SAve:
    i_fig.savefig(template_figure_directory + i_mseed_basename + '.png')
    
    plt.close(i_fig)
    
    ## SAve templates to file as mseed:
    i_template_path = template_mseed_directory + 'template_' + i_mseed_basename + '.mseed'
    i_templates[0].write(i_template_path, format="MSEED")
    
    
    
    
    
    