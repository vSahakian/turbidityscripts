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

## Csv with list of templates
template_list_path = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019_mar2023/2019_metadata_filtered.csv'

## Directory of raw/downloaded ADCP
adcp_directory = '/Users/vjs/turbidites/observational/data/BACAX_adcp_2019/netCDF'

## Directory with all the template miniseeds
template_master_dir = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019_mar2023/csv2mseed/filtered'

## Diretory to store merged/temporary adcp files per template
temp_adcp_dir = '/Users/vjs/turbidites/observational/data/templates/BACAXadcp2019'
    
## ADCP metadata path with start times ane end times of each file:
adcp_metadata_path = '{}/{}.csv'.format(adcp_directory,'metadata_startendtimes')

## Figure directories:
master_fig_dir = '/Users/vjs/turbidites/observational/data/templates/BACAX2019_adcp_turbidity_coplot'

##### OTHER PARAMS ####

#### Time deltas:
## Time to plot before event:
pre_hours = 1.5    
    
## Time to plot after event:
post_hours = 15


# %%
############################################################################## 

## Read in the template list path:
template_list = pd.read_csv(template_list_path,skiprows=6,header=0)

## Also the ADCP metadata path:
adcp_metadata = pd.read_csv(adcp_metadata_path)
## Covert start and end times to datetime:
adcp_metadata['starttime_dt'] = adcp_metadata['starttimes'].astype('datetime64[ns]')
adcp_metadata['endtime_dt'] = adcp_metadata['endtimes'].astype('datetime64[ns]')


## Convert the time columns to datetime:
template_list['starttime'] = pd.to_datetime(template_list['time_min'], format='%Y%m%dT%H%M%S')
template_list['endtime'] = pd.to_datetime(template_list['time_max'], format='%Y%m%dT%H%M%S')
template_list['eventtime'] = pd.to_datetime(template_list['time_event'], format='%Y%m%dT%H%M%S')

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

   

    ################
    ## Get the files...
    ## Empty list for the files, if there are more than one:
    i_adcp_files_list = []
    
    print('finding start and end time for ADCP files')
    ## Find where in the starttime this snippet is (take the closest/last to the start time):
    i_startfile_index = np.where(adcp_metadata.starttime_dt <= i_plot_start)[0][-1]
    
    ## And end time (take the closest/first to the end time):
    i_endfile_index = np.where(adcp_metadata.endtime_dt >= i_plot_end)[0][0]
    
    
    ## If they are the same, append only one; otherwise, append both
    if i_startfile_index == i_endfile_index:
        ## Get the full path name:
        i_adcpstart_path = '{}/{}'.format(adcp_directory,adcp_metadata.loc[i_startfile_index]['baseFilePath'])
        i_adcp_files_list.append(i_adcpstart_path)
    else:
        i_adcpstart_path = '{}/{}'.format(adcp_directory,adcp_metadata.loc[i_startfile_index]['baseFilePath'])
        i_adcpend_path = '{}/{}'.format(adcp_directory,adcp_metadata.loc[i_endfile_index]['baseFilePath'])
        ## append
        i_adcp_files_list.append(i_adcpstart_path)
        i_adcp_files_list.append(i_adcpend_path)
    
    
    ##### Opening files...
    print('opening/merging netcdfs')
    ## If there is more than one ADCP file, open them both and concatenate into one NetCDF:
    if len(i_adcp_files_list) > 1:
               
        ## Open many with xarray:
        i_adcp_nc_merge = xr.open_mfdataset(i_adcp_files_list)
        ## save out to template adcp location:
        i_adcp_templatepath = '{}/{}.nc'.format(temp_adcp_dir,row['template#'])
        i_adcp_nc_merge.to_netcdf(i_adcp_templatepath)
        
        ## Reads back in as one file to pass to function:
        i_adcp = nc4.Dataset(i_adcp_templatepath)

            
    elif len(i_adcp_files_list) == 1:
        ## Make the path to read in the adcp data:
        i_adcp_path = i_adcp_files_list[0]
        ## Read in ADCP data:
        i_adcp = nc4.Dataset(i_adcp_path)
        
        ## Also save to the template location:
        i_adcp_templatepath = '{}/{}.nc'.format(temp_adcp_dir,row['template#'])
        shutil.copy(i_adcp_files_list[0], i_adcp_templatepath)


    print('opening turbidity')
    ## And turbidity:
    ## Get the path to the turbidity data:
    i_turbidity_base = row.fname.split('/')[-1].split('.csv')[0]
    i_turbidity_path = '{}/{}.mseed'.format(template_master_dir,i_turbidity_base)
    i_turbidity = obs.read(i_turbidity_path)

    ###################
    print('plotting')
    ## Make the figure!!!!!!
    i_adcp_figure = ont.plot_adcp_timeseries(i_adcp, i_turbidity, i_plot_start, i_plot_end,
                                             event_starttime=row.eventtime,event_plus12=True,event_plus12p4=True)  
    i_adcp_figure.show()

    ## Output figure path:
    i_adcp_path_pdf = '/Users/vjs/turbidites/observational/figs/pdf/adcpDec2019.pdf'
    i_adcp_path_png = '/Users/vjs/turbidites/observational/figs/png/adcpDec2019.png'
    
    print('saving')
    # Save:
    i_adcp_figure.savefig('{}/pdf/{}.pdf'.format(master_fig_dir,row['template#']),transparent=True)
    i_adcp_figure.savefig('{}/png/{}.png'.format(master_fig_dir,row['template#']),transparent=True)

    plt.close('all')
    