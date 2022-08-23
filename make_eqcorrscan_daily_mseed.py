#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 12:43:09 2022

@author: vjs
"""

## Convert turbidity day files into mseeds

import logging

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.utils.catalog_utils import filter_picks
#from eqcorrscan.core import day_gen

import obspy as obs

import numpy as np
import pandas as pd

import onc_download_module as odm

from glob import glob
import os

# %% PARAMTERS

# csv_directory_list = ['/Users/vjs/turbidites/observational/data/BACAX_ntu_2019/*.csv']
# ## DAta type/column name for htese templates:
# data_type_list = ['turbntu']
# ## Unit of data:
# data_unit_list = ['NTU']
# ## Station names:
# station_name_list = ['BACAX']


#### for fourmile:
## List of directories to go through to convert
csv_directory_list= ['/home/dkilb/barkley/data/BACAX/turbidityntu/*.csv','/home/dkilb/barkley/data/BACHY/turbidityntu/*.csv','/home/dkilb/barkley/data/BACUS/turbidityftu/*.csv']
## DAta type/column name for htese templates:
data_type_list = ['turbntu','turbntu','turbftu']
## Unit of data:
data_unit_list = ['NTU','NTU','FTU']
## Station names:
station_name_list = ['BACAX','BACHY','BACUS']

## dta network:
data_network = 'ONC'





# %%



# %%


## Loop through the directories to convert all files in them
for j_directory_ind in range(len(csv_directory_list)):
    ## Get the glob directory wildcard:
    j_directory = csv_directory_list[j_directory_ind]
    print('working on %s' % j_directory)

    ## Get the glob list:
    j_all_paths = glob(j_directory)
    
    ## Get the common basepath for this directory:
    j_commondir = os.path.dirname(j_all_paths[0])
    
    ## Make a miniseed directory within it, if it doesn't exist yet:
    j_mseed_unmerged_dir = j_commondir + '/mseed_unmerged/'
    if os.path.exists(j_mseed_unmerged_dir) == False:
        print('making unmerged directory %s' % j_mseed_unmerged_dir)
        os.mkdir(j_mseed_unmerged_dir)
        
    j_mseed_merged_dir = j_commondir + '/mseed_merged/'
    if os.path.exists(j_mseed_merged_dir) == False:
        print('making merged directory %s' % j_mseed_merged_dir)
        os.mkdir(j_mseed_merged_dir)
        
    ## Get the datatype:
    j_data_type = data_type_list[j_directory_ind]
    
    ## And station name:
    j_station_name = station_name_list[j_directory_ind]
    

    ## Then loop through the files to convert:
    for i_day in range(len(j_all_paths)):
        i_day_path = j_all_paths[i_day]
        
        ## Get basename:
        i_day_basename = os.path.basename(i_day_path)
        i_day_dirname = os.path.dirname(i_day_path)
    
        
        print('working on %s \n \n ' % i_day_path)
        
        i_day_df = pd.read_csv(i_day_path,names=['datetime_string',j_data_type,'QC','turbcount'],comment='#')
    
        ## Convert datetime string to datetime:
        i_day_df['datetime'] = i_day_df.datetime_string.astype('datetime64[ms]')
    
        ## convert to a stream object
        i_stream_unmerged = odm.convert_df2stream(i_day_df,j_data_type,datetime_name='datetime',network=data_network,station=j_station_name,location='',channel=j_data_type)
    
        ## Get basename without file extension:
        i_mseed_base = os.path.basename(i_day_path).split('.csv')[0]
        
        ## Save it to file - first umerged:
        i_mseed_unmerged_path = j_mseed_unmerged_dir + i_mseed_base + '.mseed'
        i_stream_unmerged.write(i_mseed_unmerged_path,format='MSEED')
        
        
        ## Do NOT merge the data, see note from eqcorrscan, section 5.2.8 "Data Gaps"
        ## "You should provide data with gaps maintained, but merged 
        ##      (e.g. run st = st.merge() before passing the data to 
        ##      those functions)."
        
        # ## Then merge it...
        # i_stream_merged = i_stream_unmerged.copy().merge()
        # 
        # ## Save it to file - thenmerged:
        # i_mseed_merged_path = j_mseed_unmerged_dir + i_mseed_base + '.mseed'
        # i_stream_merged.write(i_mseed_merged_path,format='MSEED')
        
        
        