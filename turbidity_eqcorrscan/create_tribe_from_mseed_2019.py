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

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.core.match_filter import Tribe
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


## Input template mseed directory from csv:
input_mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed/'


## Tribe output mseed directory:
#template_mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/template_mseed/'

## Template figure directory
template_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/template_mseed_plots/'

## daily mseed directory:
daily_mseed_directory = '/Users/vjs/turbidites/observational/data/BACAX_ntu_2019/mseed_unmerged/'

## Plot directory:
plot_dir_tripleplots = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/filtertest_plots/'

## Template processing parameters:
        ## Note: data are 1 min sampled, so 1/60 = 0.01667 Hz. highcut = fN = 0.008333.
        ## Want to keep things with at least an hour periodicity = 1/(60*60) Hz = 0.000277
        ##   so set lowcut to 0.0002
        ## Want to keep 40 min for waveform, so set length to be 40 min * 60 sec/min
filter_lowcut=None
filter_highcut=0.008333
template_samp_rate=1/60
template_filt_order=4  ## originally 4 in example...
template_length=60*30
template_swin='P'
event_prepick=0.05
data_pre_processed=True
data_all_horiz=False


# %%
# Get all the SAC-files associated with one event.
#sac_files = glob.glob(TEST_PATH + '/SAC/2014p611252/*')
#template_mseed_list = glob(input_mseed_directory + '*.mseed')
template_mseed_list = [input_mseed_directory+'reformated_BACAX_ntu_2019_template_notfiltered_2_15_1.1254_1_0_0.5_20191224T035330_signal_only.mseed']


## Make empty list for catalog to append to:
tribe_catalog_list = []

## Make empty stream to append to:
tribe_catalog_stream = obs.Stream()

for i_file in range(len(template_mseed_list)):
    ## Define file path:
    i_mseed_path = template_mseed_list[i_file]
    
    print('working on %s \n' % i_mseed_path)
    ## Basename:
    i_mseed_basename = os.path.basename(i_mseed_path).split('.mseed')[0]
    
    ## Read in stream object:
    i_stream = obs.core.read(i_mseed_path)
    
    ## Append merged stream to tribe stream list:
    tribe_catalog_stream += i_stream.merge()
    
    # ## Plot the data:
    # i_fig = plt.figure()

    
    # plt.plot(i_stream[0].times("utcdatetime"),i_stream[0].data,color='orange',label='Unfiltered')
    
    ### Create Event for catalog:
    # Now we need to create an event!
    i_event = Event()
    i_event.origins.append(Origin())
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
    
    ## Add the pick
    i_event.picks.append(i_pick)
    
    ## Append to catalog list:
    tribe_catalog_list.append(i_event)
        
    # ## Generate template:
    # i_templates = template_gen(
    #     method='from_meta_file', meta_file=i_catalog, st=i_stream, lowcut=filter_lowcut, highcut=filter_highcut,
    #     samp_rate=template_samp_rate, filt_order=template_filt_order, length=template_length, 
    #     swin=template_swin, prepick=event_prepick, delayed=False,
    #     pre_processed=data_pre_processed,all_horiz=data_all_horiz)
    
    # plt.plot(i_templates[0][0].times("utcdatetime"),i_templates[0][0].data,color='green',label='Template')

    # ## Add legend to figure:
    # plt.legend()
    # ## And title
    # i_title_statement = '%s'
    # if filter_lowcut == None and filter_highcut != None:
    #     plt.title('%s, \n highcut %.3f Hz' % (i_mseed_basename,filter_highcut))
    # elif filter_lowcut == None and filter_highcut == None:
    #     plt.title('%s' % (i_mseed_basename))
    # else:
    #     plt.title('%s, \n lowcut %.3f Hz, highcut %.3f Hz' % (i_mseed_basename,filter_lowcut,filter_highcut))

    # ## Set axis formats:
    # #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    # #plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    
    # ## SAve:
    # i_fig.savefig(template_figure_directory + i_mseed_basename + '.png')
    
    # plt.close(i_fig)
    
    # ## SAve templates to file as mseed:
    # i_template_path = template_mseed_directory + 'template_' + i_mseed_basename + '.mseed'
    # i_templates[0].write(i_template_path, format="MSEED")
    
## Make a catalog object with the events:
tribe_catalog = Catalog(tribe_catalog_list)    

## Make the tribe object:
tribe2019turbidity = Tribe().construct(
    method="from_meta_file", meta_file=tribe_catalog, st=tribe_catalog_stream, 
    lowcut=filter_lowcut, highcut=filter_highcut, samp_rate=template_samp_rate, 
    length=template_length, filt_order=template_filt_order, 
    prepick=event_prepick, swin=template_swin, process_len=86400) 
    

    

# %%  MATCH FILTER

## Read in each daily file, rpocess, cross correlate
daily_mseed_list = sorted(glob(daily_mseed_directory + '*.mseed'))
# daily_mseed_list = [daily_mseed_directory + 'BarkleyCanyon_BarkleyCanyonAxis_variables_TurbidityNTU_20191224T000000Z_20191225T000000Z-clean_avg1minute.mseed']

crosscorr_stream_unmerged = obs.Stream()

for i_day_index in range(len(daily_mseed_list)):
    ## Get the path
    i_day_path = daily_mseed_list[i_day_index]
    
    ## REad it in:
    i_day_raw_st = obs.read(i_day_path)
    
    ## Merge it:
    i_day_raw_merged_st = i_day_raw_st.merge()

    
    ## Pre-process using same filters as templates:
   # i_day_st = shortproc(st=i_day_raw_merged_st, lowcut=filter_lowcut, highcut=filter_highcut, 
                            # filt_order=template_filt_order, samp_rate=template_samp_rate,
                            # parallel=False, num_cores=2)
                            
    ## Add to stream to use in cross correlation:
    crosscorr_stream_unmerged += i_day_raw_merged_st
    
# ## Cross correlate with tribe:
# party,families = tribe2019turbidity.detect(stream=crosscorr_stream, threshold=5, 
#                                            threshold_type='MAD', trig_int=60*3, 
#                                            plot=False, plotdir=None, 
#                                            daylong=True, parallel_process=False,
#                                            xcorr_func=None, concurrency=None, cores=None,
#                                            ignore_length=False, ignore_bad_data=False, 
#                                            group_size=None, overlap="calculate", 
#                                            full_peaks=False, save_progress=False,
#                                            process_cores=None, pre_processed=False)

## Try merging the stream.
crosscorr_stream = crosscorr_stream_unmerged.merge()

## Cross correlate with tribe:
party = tribe2019turbidity.detect(stream=crosscorr_stream, threshold=8.0, 
                                           threshold_type='absolute', trig_int=60*3, 
                                           plot=True, plotdir=plot_dir_tripleplots,
                                           daylong=True, parallel_process=False,
                                           pre_processed=False, ignore_length=True)

# %%
## Plot detections:
detectionstreams = party[0].extract_streams(stream=crosscorr_stream, length=10*60, prepick=event_prepick)
print(party[0].detections[0])
fig = detectionstreams[party[0].detections[7].id].plot(equal_scale=False, size=(800, 600))
    