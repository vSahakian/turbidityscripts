#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 13:15:52 2022

@author: vjs
"""
import logging

from obspy.clients.fdsn import Client
from obspy.core.event import Catalog

from eqcorrscan.utils.catalog_utils import filter_picks
from eqcorrscan.core import template_gen

import obspy as obs

import numpy as np
import pandas as pd

import onc_download_module as odm

from glob import glob

import matplotlib.pyplot as plt

#################################
# %% PARAMTERS/PATHS

## Template csv directory:
#csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv/'
#csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019_mar2023/csv/filtered/'
# csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2018/csv/filtered/'
# csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv/filtered/'
# csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2020/csv/filtered/'
# csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2021/csv/filtered/'
csv_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2022/csv/filtered/'


   
## Template output mseed directory:
#mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed/'
#mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019_mar2023/csv2mseed/filtered/'
# mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2018/csv2mseed/filtered/'
# mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed/filtered/'
# mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2020/csv2mseed/filtered/'
# mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2021/csv2mseed/filtered/'
mseed_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2022/csv2mseed/filtered/'


## Template output sac directory:
#sac_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu/sac/'

## Template figure directory:
#mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed_plots/'
#mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019_mar2023/csv2mseed_plots/filtered/'
# mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2018/csv2mseed_plots/filtered/'
# mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2019/csv2mseed_plots/filtered/'
# mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2020/csv2mseed_plots/filtered/'
# mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2021/csv2mseed_plots/filtered/'
mseed_figure_directory = '/Users/vjs/turbidites/observational/data/templates/BACAXturbidityntu2022/csv2mseed_plots/filtered/'



## DAta type/column name for htese templates:
data_type = 'turbntu'

## Unit of data:
data_unit = 'NTU'

## dta network:
data_network = 'ONC'
station_name = 'BACAX'

## Figure size:
plotsize = (8,4)

    
# %%

def convert_template_2mseed(csv_directory,data_type,data_network,station_name,mseed_directory,mseed_figure_directory):
    '''
    

    Parameters
    ----------
    csv_directory : string
       Path to the csv template files
    data_type : string
        What type of data are they (i.e., 'turbntu'), according to ONC codes
    data_network : string
        The network the data are from, i.e., 'ONC'
    station_name : string
        Name of the station for the mseed
    mseed_directory : string
        Path to the directory where the mseeds are to be stored
    mseed_figure_directory : string
        Path to the directory where the figures of the mseed templates are to be stored

    Returns
    -------
    None.

    '''
    ## Loop through the files in teh template directory:
    all_csv_paths = glob(csv_directory + '*.csv')
    
    ## Loop through them:
    for i_template in range(len(all_csv_paths)):
        i_template_path = all_csv_paths[i_template]
        
        print('working on %s' % i_template_path)
        
        i_template_df = pd.read_csv(i_template_path,names=['datetime_string',data_type])
    
        ## Convert datetime string to datetime:
        i_template_df['datetime'] = i_template_df.datetime_string.astype('datetime64[ms]')
    
        ## convert to a stream object
        i_stream = odm.convert_df2stream(i_template_df,data_type,datetime_name='datetime',network=data_network,station=station_name,location='',channel=data_type)
    
        ## Save it to file:
        i_mseed_base = i_template_path.split('/')[-1].split('.csv')[0]
        i_mseed_path = mseed_directory + i_mseed_base + '.mseed'
        i_stream.write(i_mseed_path,format='MSEED')
        
        ## Plot:
        i_plot_path = mseed_figure_directory + i_mseed_base + '.png'
        
        ## INitiate plot and other commands before plotting:
        i_fig, i_ax = plt.subplots(figsize=plotsize)
        i_ax.set_ylabel(data_unit)
        i_ax.xaxis.set_ticklabels([])
        i_ax.yaxis.set_ticklabels([])
        
        i_stream.plot(fig=i_fig)
    
        
        print('saving to %s \n \n \n' % i_plot_path)
        ## Save plot:
        i_fig.savefig(i_plot_path)
    
        ## Close plot:
        plt.close('all')


convert_template_2mseed(csv_directory, data_type, data_network, station_name, mseed_directory, mseed_figure_directory)


