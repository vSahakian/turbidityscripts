#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:11:19 2023

@author: vjs
"""

import numpy as np
import matplotlib.pyplot as plt
import onctools as ont
import pandas as pd    
from datetime import datetime, timedelta
import netCDF4 as nc4
import obspy as obs
import xarray as xr
import shutil



##### PARAMETERS and paths

## Year to process
# templateyear = '2018'
# templateyear = '2019'
# templateyear = '2020'
templateyear = '2021'
# templateyear = '2022'

## Main directory for year
template_main_dir = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu{}'.format(templateyear)
## Directory with all the template miniseeds
template_master_dir = '{}/csv2mseed/filtered'.format(template_main_dir)
## Csv with list of templates
template_list_path = '{}/{}_metadata_filtered.csv'.format(template_main_dir,templateyear)
## Directory to store figures
master_fig_dir = '{}/multidata_plots'.format(template_main_dir)


####
## Master time series directory
time_series_masterdir = '/Users/vjs/turbidites/observational/data/barkleycanyon/BACAX'


## Metadata paths with start times ane end times of each file:
chloro_metadata_path = '{}/chlorophyll/metadata/download_metadata_v0.csv'.format(time_series_masterdir)
temp_metadata_path = '{}/seawatertemperature/metadata/download_metadata_v0.csv'.format(time_series_masterdir)
oxy_metadata_path = '{}/oxygen/metadata/download_metadata_v0.csv'.format(time_series_masterdir)



##### OTHER PARAMS ####

#### Time deltas:
## Time to plot before event:
pre_hours = 1.5    
    
## Time to plot after event:
post_hours = 1.5


## Figure dimensions:
figure_dims = (16,13)

# %%
############################################################################## 

## Read in the template list path:
template_list = pd.read_csv(template_list_path,skiprows=6,header=0)

# ## Also the ADCP metadata path:
# adcp_metadata = pd.read_csv(adcp_metadata_path)
# ## Covert start and end times to datetime:
# adcp_metadata['starttime_dt'] = adcp_metadata['starttimes'].astype('datetime64[ns]')
# adcp_metadata['endtime_dt'] = adcp_metadata['endtimes'].astype('datetime64[ns]')


## Convert the time columns to datetime:
template_list['starttime'] = pd.to_datetime(template_list['time_min'], format='%Y%m%dT%H%M%S')
template_list['endtime'] = pd.to_datetime(template_list['time_max'], format='%Y%m%dT%H%M%S')
template_list['eventtime'] = pd.to_datetime(template_list['time_event'], format='%Y%m%dT%H%M%S')



## Also open the metadata for the multidata types:
chlorophyll_metadata = pd.read_csv(chloro_metadata_path).drop(labels=['Unnamed: 0.1'],axis=1)
temp_metadata = pd.read_csv(temp_metadata_path).drop(labels=['Unnamed: 0.1'],axis=1)
oxy_metadata = pd.read_csv(oxy_metadata_path).drop(labels=['Unnamed: 0.1'],axis=1)

## Convert the time columns to datetime:
chlorophyll_metadata['dateFromDt'] = chlorophyll_metadata['dateFrom'].astype('datetime64[ns]')
chlorophyll_metadata['dateToDt'] = chlorophyll_metadata['dateTo'].astype('datetime64[ns]')

temp_metadata['dateFromDt'] = temp_metadata['dateFrom'].astype('datetime64[ns]')
temp_metadata['dateToDt'] = temp_metadata['dateTo'].astype('datetime64[ns]')

oxy_metadata['dateFromDt'] = oxy_metadata['dateFrom'].astype('datetime64[ns]')
oxy_metadata['dateToDt'] = oxy_metadata['dateTo'].astype('datetime64[ns]')


# %%
print('iterating through templates...')
# Iterate through, get the turbidity template name, as well as the start time
for index,row in template_list.iterrows():
    print('On template {}'.format(row.fname))
    
    ## Get the file base name:
    i_basename = template_list.loc[index]['fname'].split('/')[-1].split('.csv')[0]
    
    ## Get the path for the template:
    i_template_path = template_master_dir + 'csv2mseed/filtered/' + i_basename + '.mseed'
    
    ## Get the start time of the turbidity file:
    i_eventtime = row.eventtime
  
    ## Get the time to start using and downloading in teh ADCP data:
    i_plot_start = i_eventtime - timedelta(hours=pre_hours)
    i_plot_end = i_eventtime + timedelta(hours=post_hours)

   
    ## Get the chlorophyll files to plot:
    i_chloro_basepaths = ont.get_template_files4plot(i_plot_start,i_plot_end,chlorophyll_metadata)
    i_temp_basepaths = ont.get_template_files4plot(i_plot_start,i_plot_end,temp_metadata)
    i_oxy_basepaths = ont.get_template_files4plot(i_plot_start,i_plot_end,oxy_metadata)

    ## Modify paths to have correect full path:
    for j_path in range(len(i_chloro_basepaths)):
        i_chloro_basepaths[j_path] = f'{time_series_masterdir}/chlorophyll/{i_chloro_basepaths[j_path]}'
        
    for j_path in range(len(i_temp_basepaths)):
        i_temp_basepaths[j_path] = f'{time_series_masterdir}/seawatertemperature/{i_temp_basepaths[j_path]}'
        
    for j_path in range(len(i_oxy_basepaths)):
        i_oxy_basepaths[j_path] = f'{time_series_masterdir}/oxygen/{i_oxy_basepaths[j_path]}'
    
    
    ##### Opening files...
    print ('opening multidata files')
    i_chloro_df = ont.open_merge_csv4multiplot(i_chloro_basepaths,'chlorophyll')
    i_temp_df = ont.open_merge_csv4multiplot(i_temp_basepaths,'seawatertemperature')
    i_oxy_df = ont.open_merge_csv4multiplot(i_oxy_basepaths,'oxygen')


    print('opening turbidity')
    ## And turbidity:
    ## Get the path to the turbidity data:
    i_turbidity_base = row.fname.split('/')[-1].split('.csv')[0]
    i_turbidity_path = '{}/{}.mseed'.format(template_master_dir,i_turbidity_base)
    i_turbidity = obs.read(i_turbidity_path)

    ## Turn into a dataframe for the plotting script:
    # First turn utcdatetimes of obspy stream/trace to datetimes:
    i_turb_utcdt = i_turbidity[0].times('UTCDateTime')
    i_turb_datetimes = [datetime.utcfromtimestamp(utcdatetime.timestamp) for utcdatetime in i_turb_utcdt]
    
    i_turbidity_df = pd.DataFrame({'datetime':i_turb_datetimes,'turbidityntu':i_turbidity[0].data})

    ###################
    print('plotting')
    
    ## Make the figure!!!!!!
    i_figure = ont.plot_multidata_timeseries(i_oxy_df,i_chloro_df,i_temp_df,i_turbidity_df,figure_dims,i_plot_start,i_plot_end,i_eventtime)




    ## Output figure path:
    i_template_num = row['template#']
    i_path_pdf = f'{master_fig_dir}/pdf/template{i_template_num}_multidata.pdf'
    i_path_png = f'{master_fig_dir}/png/template{i_template_num}_multidata.png'
    
    print('saving')
    # Save:
    i_figure.savefig(i_path_pdf,transparent=True)
    i_figure.savefig(i_path_png,transparent=True)

    plt.close('all')
    