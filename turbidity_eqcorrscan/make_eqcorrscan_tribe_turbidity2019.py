#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:25:06 2022

@author: vjs
"""
## Make templates within eqcorrscan for 2019

import logging

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.core.template_gen import template_gen
from eqcorrscan.core.match_filter import Tribe
from eqcorrscan.utils.pre_processing import shortproc

#from eqcorrscan.core import day_gen

import obspy as obs

import numpy as np
import pandas as pd

import onc_download_module as odm

from glob import glob
import os

# %% 

## Template output mseed directory:
template_mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/template_mseed/'

daily_mseed_directory = '/Users/vjs/turbidites/observational/data/BACAX_ntu_2019/mseed_unmerged/'

## 
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

## Glob together the templates into a list to open:
# template_path_list = glob(template_mseed_directory + '*.mseed')
template_path_list = [template_mseed_directory+'template_reformated_BACAX_ntu_2019_template_notfiltered_2_15_1.1254_1_0_0.5_20191224T035330_signal_only.mseed']

## Make an empty list to put templates into:
template_object_list = []

## Loop through them:
for i_template_index in range(len(template_path_list)):
    ## Get the template path:
    i_template_path = template_path_list[i_template_index]
    
    ## Read in the template:
    i_template_st = obs.read(i_template_path)
    
    ## Append to list:
    template_object_list.append(i_template_st)
    
## Make tribe object:
tribe2019turbidity = Tribe(templates=template_object_list)
    
    
## Read in each daily file, rpocess, cross correlate
# daily_mseed_list = glob(daily_mseed_directory + '*.mseed')
daily_mseed_list = [daily_mseed_directory + 'BarkleyCanyon_BarkleyCanyonAxis_variables_TurbidityNTU_20191225T000000Z_20191226T000000Z-clean_avg1minute.mseed']


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
    
    ## Cross correlate with tribe:
    party,families = tribe2019turbidity.detect(i_day_raw_merged_st, threshold=5, threshold_type='MAD', trig_int=60*3, 
                                               plot=False, plotdir=None, 
                                               daylong=True, parallel_process=True,
                                               xcorr_func=None, concurrency=None, cores=None,
                                               ignore_length=False, ignore_bad_data=False, 
                                               group_size=None, overlap="calculate", 
                                               full_peaks=False, save_progress=False,
                                               process_cores=None, pre_processed=True)

    
    